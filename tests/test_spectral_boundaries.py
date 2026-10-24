# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Spectral loaders and optional dependency boundaries without network access."""

import dataclasses
import sys

import numpy as np
import pytest

import optmat
from optmat.conditions import State
from optmat.datasets import load_source
from optmat.exceptions import DataError
from optmat.gas import absorption, eos, euv, virial
from optmat.gas.gases import Gas
from optmat.gas.models import fitted


@pytest.fixture
def sample(tmp_path):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[datasets]]\nid="sigma"\n'
        'material="Ar"\nquantity="absorption_cross_section"\nunit="Mb"\n'
        'abscissa_unit="nm"\npoints=[["b","500","2"],["a","400","1"]]\n',
        encoding="utf-8",
    )
    return load_source(path).datasets[0]


def test_absorption_tables_sort_wavelengths_and_convert_area(sample, monkeypatch):
    monkeypatch.setattr(absorption._catalogue, "by_key", lambda: {sample.key: sample})
    wavelength, area = absorption._read_dataset(sample.key)
    # Megabarns to square metres: 1 Mb = 1e-22 m^2.
    np.testing.assert_allclose(wavelength, [400, 500], rtol=2e-15, atol=0)
    np.testing.assert_allclose(area, [1e-22, 2e-22], rtol=2e-15, atol=0)
    with pytest.raises(DataError, match="unknown absorption dataset"):
        absorption._read_dataset("missing")
    wrong = dataclasses.replace(sample, quantity="photoionization_cross_section")
    monkeypatch.setattr(absorption._catalogue, "by_key", lambda: {sample.key: wrong})
    with pytest.raises(DataError, match="not an absorption_cross_section"):
        absorption._read_dataset(sample.key)


def test_absorption_configuration_requires_a_known_variant_and_complete_pieces(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(absorption, "ABSORPTION_DIRECTORY", tmp_path)
    absorption.measured_absorption.cache_clear()
    try:
        with pytest.raises(DataError, match="no measured absorption spectrum"):
            absorption.measured_absorption("Ar")
        path = tmp_path / "Ar.toml"
        for text, message in [
            ("pieces={}\n", "pieces must be an array"),
            ("pieces=[{}]\n", "a piece lacks"),
            ("pieces=[]\ntemperature_k=300\n", "no variant"),
        ]:
            path.write_text(text, encoding="utf-8")
            with pytest.raises(DataError, match=message):
                absorption.measured_absorption("Ar", variant="missing")
    finally:
        absorption.measured_absorption.cache_clear()


def test_reloading_a_fit_detects_changed_absorption_provenance(monkeypatch):
    spectrum = absorption.measured_absorption("O3")
    monkeypatch.setattr(fitted, "measured_absorption", lambda *args, **kwargs: spectrum)
    record = {
        "gas": "O3",
        "name": "synthetic",
        "absorption": {
            "fwhm_cm": spectrum.fwhm_cm,
            "pieces": [[p.source, p.from_nm, p.to_nm] for p in spectrum.pieces],
        },
    }
    assert fitted._absorption(record) is spectrum
    record["absorption"]["pieces"] = []
    with pytest.raises(DataError, match="pieces changed since the fit"):
        fitted._absorption(record)


def test_continua_require_photon_energies_and_supported_quantities(
    sample, monkeypatch, curated
):
    curated(sample)
    with pytest.raises(DataError, match="needs photon energies in eV"):
        euv._continuum_part({"dataset": sample.key}, Gas.Ar)
    energy = dataclasses.replace(sample, abscissa_unit="eV")
    curated(energy)
    with pytest.raises(DataError, match="not a continuum"):
        euv._continuum_part({"dataset": sample.key}, Gas.Ar)
    with pytest.raises(ValueError, match="energies must increase"):
        euv.OscillatorDistribution(
            Gas.Ar,
            "synthetic",
            np.array([1.0]),
            np.array([2.0]),
            np.array([3.0, 2.0]),
            np.array([1.0, 1.0]),
        )


def test_missing_gas_extra_has_an_actionable_error_and_ideal_density_still_works(
    monkeypatch,
):
    eos._coolprop.cache_clear()
    monkeypatch.setitem(sys.modules, "CoolProp", None)
    monkeypatch.setitem(sys.modules, "CoolProp.CoolProp", None)
    state = State(100000.0, 300.0)
    # Exercise the uncached dependency boundary so prior real-gas calls cannot hide it.
    try:
        with pytest.raises(ImportError, match=r"install optmat\[gas\]"):
            eos._coolprop()
    finally:
        monkeypatch.undo()
        eos._coolprop.cache_clear()
    assert eos.number_density("Ar", state, ideal=True) == 100000 / (1.380649e-23 * 300)


def test_unknown_virial_dataset_is_rejected(tmp_path, monkeypatch):
    path = tmp_path / "virial.toml"
    path.write_text('[Ar]\nB_R="missing"\nA_R="missing"\n', encoding="utf-8")
    monkeypatch.setattr(virial, "VIRIAL_FILE", path)
    virial.virial_sources.cache_clear()
    try:
        with pytest.raises(DataError, match="unknown dataset"):
            virial.virial_sources()
    finally:
        virial.virial_sources.cache_clear()


def test_lazy_package_submodules_are_the_imported_modules():
    assert optmat.__getattr__("datasets").load_source is load_source


def test_unsmoothed_absorption_keeps_the_measured_grid(sample, monkeypatch):
    monkeypatch.setattr(absorption._catalogue, "by_key", lambda: {sample.key: sample})
    spectrum = absorption.assemble(
        "Ar",
        (absorption.Piece(sample.key, 400, 500, "synthetic"),),
        temperature_k=300,
        fwhm_cm=0.0,
        step_cm=5000,
    )
    # Grid points are 0,5000,...,25000 cm-1. The only samples inside are 500 and 400 nm.
    np.testing.assert_allclose(
        spectrum.sigma_m2, [0, 0, 0, 0, 2e-22, 1e-22], rtol=2e-15, atol=0
    )


def test_nist_lines_skip_excited_lower_states_and_use_observed_wavelength(
    tmp_path, monkeypatch
):
    path = tmp_path / "external" / "nist_asd" / "Ar_I_lines.tsv"
    path.parent.mkdir(parents=True)
    path.write_text(
        "fik\tJ_k\tEi(cm-1)\tritz_wl_vac(nm)\tobs_wl_vac(nm)\tconf_k\n"
        "0.2\t1\t10\t100\t\texcited\n"
        "0.3\t1\t0\t\t200\tground\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(euv, "DATA_DIRECTORY", tmp_path)
    lines = euv.nist_resonance_lines(Gas.Ar)
    assert len(lines) == 1
    energy, strength, label = lines[0]
    assert energy == pytest.approx(45.56335252767 / 200, rel=2e-15, abs=0.0)
    assert strength == 0.3
    assert label == "ground"


def test_unavailable_phase_matching_partner_is_omitted_and_missing_state_is_rejected(
    sample, monkeypatch
):
    from optmat.exceptions import OptmatError  # noqa: PLC0415
    from optmat.gas import crosschecks, reduction  # noqa: PLC0415

    with pytest.raises(DataError, match="no state"):
        reduction.density_at_state(sample)
    ratio = dataclasses.replace(
        sample,
        quantity="phase_match_density_ratio",
        unit="1",
        parameters={"reference_material": "He", "harmonic": 2},
    )
    monkeypatch.setattr(crosschecks, "datasets", lambda **kwargs: (ratio,))

    def unavailable(gas):
        raise OptmatError("no recommended model")

    monkeypatch.setattr(crosschecks, "default_model", unavailable)
    assert crosschecks.phase_match_ratios() == []


def test_a_listed_euv_strength_equal_to_the_total_adds_no_pseudo_line(
    sample, monkeypatch, curated
):
    listed = dataclasses.replace(
        sample, quantity="oscillator_strength", unit="1", abscissa_unit="eV"
    )
    total = dataclasses.replace(
        listed, id="total", rows=("0",), abscissae=("1",), values=("3",), sigmas=("",)
    )
    curated(listed, total)
    terms = euv._lines({"dataset": listed.key, "total": total.key}, Gas.Ar)
    assert len(terms) == 2
    assert sum(f for _, f in terms) == 3


def test_zero_order_absorbing_dispersion_is_frequency_times_polarizability():
    from optmat.gas.models.oscillator import OscillatorModel  # noqa: PLC0415

    spectrum = absorption.MeasuredAbsorption(
        Gas.Ar,
        300.0,
        (),
        0.0,
        10000.0,
        np.array([0.0, 10000.0, 20000.0]),
        np.zeros(3),
        np.full(3, 2.0),
    )
    model = OscillatorModel(
        gas=Gas.Ar,
        name="synthetic",
        citation="synthetic",
        wavelength_range=(200e-9, 2e-6),
        kind="fitted",
        strengths=(1.0,),
        energies=(1.0,),
        labels=("pole",),
        offset_au=0.1,
        retardation_au=0.2,
        absorption=spectrum,
    )
    w = np.array([0.01, 0.02])
    # The supplied spectrum has constant Re alpha=2. The order-zero dispersion is
    # omega * [1/(1-omega^2) + 0.1 + 0.2 omega^2 + 2].
    expected = w * (1 / (1 - w * w) + 0.1 + 0.2 * w * w + 2)
    np.testing.assert_allclose(model.dispersion_au(w, 0), expected, rtol=2e-15, atol=0)

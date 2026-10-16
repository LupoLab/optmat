# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The two-photon ionisation models (OEP 6, Phase 1): derivation, values and identity.

Each model (``data/gas/models/two_photon/<gas>_<polarization>.json``) is derived by
``tools/fit/two_photon_record.py`` from ``data/gas/two_photon.toml`` and the curated
datasets. These tests check that the stored models are fresh, that a reference
reproduces its record and the numbers its paper prints, the thresholds against NIST, the
regions and their values (zero below threshold, NaN where unavailable, never an error on
a grid), the conversions of OEP 6 section 4, and the identity.
"""

import dataclasses
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from optmat.constants import BAR, E_CHARGE, HBAR, K_B, C
from optmat.datasets import datasets
from optmat.exceptions import (
    ArgumentError,
    DataError,
    ModelRangeWarning,
    UnavailableModelError,
)
from optmat.gas import _two_photon_derivation as derivation
from optmat.gas import two_photon
from optmat.gas.gases import Gas
from optmat.gas.two_photon import (
    TwoPhotonIonisationModel,
    TwoPhotonRegion,
    two_photon_ionisation_model,
    two_photon_ionisation_models,
)

MODELS = sorted((g.value, p) for g, p in two_photon_ionisation_models())
# Products and sums of a few floats: they agree to rounding, about 1e-16 relative.
RTOL = 1e-12
DASH = "\N{EN DASH}"  # the manual's ranges


def _omega(ev):
    """The angular frequency of a photon of an energy in eV, rad/s."""
    return ev * E_CHARGE / HBAR


def _nm(omega: float) -> float:
    return 2 * math.pi * C / omega * 1e9


def _quiet(model: TwoPhotonIonisationModel, omega):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return model.cross_section(omega)


def _assert_same_record(fresh, stored, path="model"):
    # Platforms differ in the last bits of the derived numbers (Linux and macOS by
    # 1 ulp, CI 2026-10-16): numbers to 1e-12 relative, everything else exactly, and a
    # stored zero stays zero.
    if isinstance(stored, float) or isinstance(fresh, float):
        assert isinstance(fresh, int | float) and isinstance(stored, int | float), path
        assert math.isclose(fresh, stored, rel_tol=1e-12, abs_tol=0.0), path
    elif isinstance(stored, dict):
        assert isinstance(fresh, dict) and fresh.keys() == stored.keys(), path
        for key in stored:
            _assert_same_record(fresh[key], stored[key], f"{path}.{key}")
    elif isinstance(stored, list):
        assert isinstance(fresh, list) and len(fresh) == len(stored), path
        for i, (now, before) in enumerate(zip(fresh, stored, strict=True)):
            _assert_same_record(now, before, f"{path}[{i}]")
    else:
        assert fresh == stored, path


@pytest.mark.parametrize(("gas", "polarization"), MODELS)
def test_the_stored_models_are_fresh(gas, polarization):
    stored = (two_photon._DIRECTORY / f"{gas}_{polarization}.json").read_text(
        encoding="utf-8"
    )
    fresh = two_photon._dumps(derivation.derive(gas, polarization))
    _assert_same_record(json.loads(fresh), json.loads(stored))


def test_the_models_are_the_approved_references():
    # TC-943: the author's rule-8 approvals.
    expected = {
        ("He", "linear"): (
            "mihelic_2021/fig3a_sigma2_L0",
            "mihelic_2021/fig3a_sigma2_L2",
        ),
        ("Ne", "linear"): ("mckenna_2004/fig1_sigma2_Ne",),
        ("Ar", "linear"): ("mckenna_2004/fig5_sigma2_Ar",),
        ("Xe", "linear"): ("lagutin_2017/fig8_cipfcp_linear",),
        ("Xe", "circular"): ("lagutin_2017/fig8_cipfcp_circular",),
    }
    got = {
        (g.value, p): m.reference.datasets
        for (g, p), m in two_photon_ionisation_models().items()
    }
    assert got == expected


@pytest.mark.parametrize(
    ("gas", "nm"),
    [
        # The thresholds of the manual's table (TC-902), to its printed 0.01 nm.
        ("He", (100.85,)),
        ("Ne", (114.99, 114.48)),
        ("Ar", (157.34, 155.59)),
        ("Xe", (204.43, 184.55)),
    ],
)
def test_the_thresholds_are_nists(gas, nm):
    model = two_photon_ionisation_model(gas)
    assert [round(_nm(t.omega), 2) for t in model.thresholds] == list(nm)
    assert all(t.photons == 2 for t in model.thresholds)
    assert model.e1.photons == 1


def test_argon_threshold_from_the_printed_nist_values():
    # 15.7596119 eV and the Ar+ 2P1/2 level 1431.5831 cm^-1 (nist_asd_2026), halved.
    t = two_photon_ionisation_model("Ar").thresholds
    ev = t[0].energy / E_CHARGE
    assert ev == pytest.approx(15.7596119, rel=RTOL)
    upper = 15.7596119 + 1431.5831 * 2 * math.pi * HBAR * C * 100 / E_CHARGE
    assert t[1].energy / E_CHARGE == pytest.approx(upper, rel=RTOL)
    assert t[1].omega == pytest.approx(t[1].energy / (2 * HBAR), rel=RTOL)


def test_krypton_has_no_model():
    with pytest.raises(UnavailableModelError, match="Q2"):
        two_photon_ionisation_model("Kr")


@pytest.mark.parametrize(
    ("gas", "polarization", "error", "match"),
    [
        ("He", "circular", UnavailableModelError, "no approved circular"),
        ("Ar", "elliptical", ArgumentError, "polarization"),
        ("N2", "linear", UnavailableModelError, "noble gases"),
    ],
)
def test_refusals(gas, polarization, error, match):
    with pytest.raises(error, match=match):
        two_photon_ionisation_model(gas, polarization)


@pytest.mark.parametrize(("gas", "polarization"), MODELS)
def test_the_reference_reproduces_its_vertices(gas, polarization):
    model = two_photon_ionisation_model(gas, polarization)
    omega, sigma = derivation.summed_curve(model.reference.datasets)
    keep = (omega >= model.thresholds[0].omega) & (omega <= model.e1.omega)
    np.testing.assert_array_equal(model.omega, omega[keep])
    np.testing.assert_allclose(
        _quiet(model, model.omega[:-1]), sigma[keep][:-1], rtol=RTOL
    )


def test_staircase_vertices_become_their_geometric_mean():
    # McKenna's Ar drawing has vertical steps (the drawing's 0.13 pt quantum): several
    # vertices at one photon energy, which the model takes as one, their geometric mean.
    (dataset,) = [d for d in datasets() if d.key == "mckenna_2004/fig5_sigma2_Ar"]
    w = dataset.angular_frequency()
    repeated = w[np.isin(w, w[np.diff(w, prepend=-1.0) == 0])][0]
    values = dataset.scaled_values()[w == repeated] * 1e-8
    omega, sigma = derivation.vertices(dataset)
    assert np.all(np.diff(omega) > 0)
    assert sigma[omega == repeated][0] == pytest.approx(
        math.exp(np.mean(np.log(values))), rel=RTOL
    )


@pytest.mark.parametrize(
    ("gas", "low", "high", "value", "margin"),
    [
        # mckenna_2004 (TC-856): Ar 'within 8% of 1.2e-50' below 11.5 eV, read 1.119
        # to 1.268e-50 from 8.0 to 11.4 eV (in cm^4 s; 1e-8 to m^4 s); the supported
        # band ends at 11.02 eV. Ne 'within 10% of 1.65e-51' from 11 to 16 eV.
        ("Ar", 8.0, 11.0, 1.2e-58, 0.08),
        ("Ne", 11.0, 16.0, 1.65e-59, 0.11),
    ],
)
def test_the_printed_plateaus(gas, low, high, value, margin):
    model = two_photon_ionisation_model(gas)
    ev = np.linspace(low, high, 101)
    sigma = _quiet(model, _omega(ev))
    assert np.all(np.abs(sigma / value - 1) <= margin)


def test_regions_and_values_on_a_propagation_grid():
    argon = two_photon_ionisation_model("Ar")
    nm = np.linspace(300, 60, 2000)
    omega = 2 * math.pi * C / (nm * 1e-9)
    with pytest.warns(ModelRangeWarning, match="unsupported-detail"):
        sigma = argon.cross_section(omega)
    labels = argon.region(omega)
    available = argon.available(omega)
    assert set(labels) == {
        "below_threshold",
        "unavailable",
        "unsupported_detail",
        "supported",
    }
    assert np.all(sigma[labels == "below_threshold"] == 0)
    assert np.all(np.isnan(sigma[~available]))
    assert np.all(np.isfinite(sigma[available]))
    assert np.all(sigma[labels == "supported"] > 0)
    # Unavailable: between the threshold and the first vertex, and at and beyond E1.
    assert not argon.available(argon.e1.omega)
    assert not argon.available(argon.thresholds[0].omega)
    assert bool(argon.available(argon.omega[0]))


def test_supported_values_do_not_warn():
    argon = two_photon_ionisation_model("Ar")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert float(argon.cross_section(_omega(9.0))) > 0


@pytest.mark.parametrize("bad", [0.0, -1.0, math.inf, math.nan])
def test_frequencies_must_be_finite_and_positive(bad):
    argon = two_photon_ionisation_model("Ar")
    with pytest.raises(ArgumentError, match="frequencies"):
        argon.cross_section([_omega(9.0), bad])


@pytest.mark.parametrize("bad", [-1.0, math.inf])
def test_intensities_and_densities_must_be_finite_and_non_negative(bad):
    argon = two_photon_ionisation_model("Ar")
    with pytest.raises(ArgumentError, match="intensity"):
        argon.rate(_omega(9.0), bad)
    with pytest.raises(ArgumentError, match="number_density"):
        argon.loss_coefficient(_omega(9.0), bad)


def test_rate_loss_and_pair_conversions():
    argon = two_photon_ionisation_model("Ar")
    w = _omega(np.array([8.5, 9.0, 10.0]))
    sigma = argon.cross_section(w)
    intensity = 1e17  # W/m^2
    flux = intensity / (HBAR * w)
    np.testing.assert_allclose(argon.rate(w, intensity), sigma * flux**2, rtol=RTOL)
    density = BAR / (K_B * 293.15)
    beta = argon.loss_coefficient(w, density)
    # Two photons per event: dI/dz = -2 hbar omega N W = -beta I^2.
    np.testing.assert_allclose(
        beta * intensity**2,
        2 * HBAR * w * density * argon.rate(w, intensity),
        rtol=RTOL,
    )
    np.testing.assert_allclose(argon.pair_cross_section(2 * w), sigma, rtol=RTOL)


def test_the_loss_law_against_an_integration():
    # dI/dz = -beta I^2 has I(z) = I0 / (1 + beta I0 z): a fourth-order Runge-Kutta
    # integration agrees to its truncation error, far below 1e-9 with this step.
    argon = two_photon_ionisation_model("Ar")
    w = _omega(8.0)
    beta = float(argon.loss_coefficient(w, BAR / (K_B * 293.15)))
    i0, length, steps = 1e17, 0.01, 1000
    h, i = length / steps, i0
    for _ in range(steps):
        k1 = -beta * i**2
        k2 = -beta * (i + h * k1 / 2) ** 2
        k3 = -beta * (i + h * k2 / 2) ** 2
        k4 = -beta * (i + h * k3) ** 2
        i += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    assert i == pytest.approx(i0 / (1 + beta * i0 * length), rel=1e-9)
    # About 5e-15 m/W at 1 bar near 155 nm (OEP 6, section 2): with McKenna's
    # 1.27e-50 cm^4 s at 8.0 eV, 2 N sigma / hbar omega = 4.9e-15; 1 %: the rounding.
    assert beta == pytest.approx(4.9e-15, rel=0.01)


def test_the_uncertainty_keeps_its_kinds_apart():
    argon = two_photon_ionisation_model("Ar").uncertainty
    assert argon.source is None
    assert argon.digitising == 0.02
    # Petrov et al.'s two gauges against McKenna over the supported band (TC-944).
    assert argon.envelope == pytest.approx((0.911, 0.992), abs=0.001)
    assert two_photon_ionisation_model("Ne").uncertainty.envelope is None
    xenon = two_photon_ionisation_model("Xe").uncertainty
    assert xenon.source == 0.077
    # The reference's velocity form is a diagnostic: Fink and Johnson alone set the
    # envelope.
    (gauge,) = [s for s in xenon.spreads if s.role == "gauge"]
    (fink,) = [s for s in xenon.spreads if s.role == "independent"]
    assert xenon.envelope == (fink.lower_ratio, fink.upper_ratio)
    assert gauge.lower_ratio > 1


def test_status_and_identity():
    neon = two_photon_ionisation_model("Ne")
    assert neon.status == "untested"
    assert "single" in neon.status_reason
    argon = two_photon_ionisation_model("Ar")
    assert argon.status == "validated"
    assert argon.identity.name == "two_photon/Ar/linear"
    assert argon.identity.kind == "two_photon_ionisation"
    assert {"mckenna_2004", "nist_asd_2026", "petrov_2016"} <= set(
        argon.identity.sources
    )
    assert argon.digest == argon.identity.digest
    changed = dataclasses.replace(argon, sigma=argon.sigma * 1.01)
    assert changed.digest != argon.digest
    assert not argon.omega.flags.writeable


def test_the_model_checks_its_grid_and_regions():
    argon = two_photon_ionisation_model("Ar")
    with pytest.raises(DataError, match="grid"):
        dataclasses.replace(argon, sigma=-argon.sigma)
    with pytest.raises(DataError, match="regions"):
        dataclasses.replace(argon, regions=argon.regions[:-1])
    bad = (*argon.regions[:-1], TwoPhotonRegion("odd", argon.e1.omega, math.inf, ""))
    with pytest.raises(DataError, match="regions"):
        dataclasses.replace(argon, regions=bad)


@pytest.fixture
def stored(monkeypatch, tmp_path):
    """Point the loader at an empty directory, and forget the cached models."""
    monkeypatch.setattr(two_photon, "_DIRECTORY", tmp_path)
    two_photon._stored.cache_clear()
    yield tmp_path
    two_photon._stored.cache_clear()


def test_a_missing_stored_model_is_an_error(stored):
    with pytest.raises(DataError, match="missing"):
        two_photon_ionisation_model("Ar")


def test_a_stored_model_of_another_gas_is_an_error(stored):
    text = two_photon._dumps(derivation.derive("Ne", "linear"))
    (stored / "Ar_linear.json").write_text(text)
    with pytest.raises(DataError, match="holds the model of Ne"):
        two_photon_ionisation_model("Ar")


def test_an_old_or_malformed_stored_model_is_an_error(stored):
    record = json.loads(two_photon._dumps(derivation.derive("Ar", "linear")))
    (stored / "Ar_linear.json").write_text(
        json.dumps(record | {"derivation_version": 0})
    )
    with pytest.raises(DataError, match="re-derive"):
        two_photon_ionisation_model("Ar")
    two_photon._stored.cache_clear()
    del record["e1"]
    (stored / "Ar_linear.json").write_text(json.dumps(record))
    with pytest.raises(DataError, match="needs exactly the fields"):
        two_photon_ionisation_model("Ar")


def test_the_stored_model_round_trips():
    for (gas, polarization), model in two_photon_ionisation_models().items():
        again = two_photon._loads(two_photon._dumps(model), "test")
        assert again.digest == model.digest, (gas, polarization)


def test_the_derivation_refuses_what_it_cannot_use(monkeypatch):
    (mccown,) = [d for d in datasets() if d.key == "mccown_1982/alpha_Xe"]
    with pytest.raises(DataError, match="not a two-photon cross section"):
        derivation.vertices(mccown)
    (argon,) = [d for d in datasets() if d.key == "mckenna_2004/fig5_sigma2_Ar"]
    negative = dataclasses.replace(argon, values=("-1", *argon.values[1:]))
    with pytest.raises(DataError, match="not positive"):
        derivation.vertices(negative)
    with pytest.raises(DataError, match="unknown datasets"):
        derivation.summed_curve(["nobody_2000/none"])
    with pytest.raises(DataError, match="do not overlap"):
        derivation.summed_curve(
            ["mckenna_2004/fig5_sigma2_Ar", "mckenna_2004/fig7_sigma2_Ar"]
        )
    with pytest.raises(DataError, match="no linear reference"):
        derivation.derive("Kr", "linear")
    selection = two_photon._selection()
    argon_entry = dict(selection["Ar"])
    argon_entry["linear"] = dict(argon_entry["linear"]) | {
        "validators": [
            {
                "label": "outside",
                "datasets": ["mckenna_2004/fig7_sigma2_Ar"],
                "form": "vector",
                "role": "independent",
            }
        ]
    }
    monkeypatch.setattr(
        derivation, "_selection", lambda: selection | {"Ar": argon_entry}
    )
    with pytest.raises(DataError, match="no overlap"):
        derivation.derive("Ar", "linear")


def test_the_registry_has_every_approved_model():
    models = two_photon_ionisation_models()
    assert sorted((g.value, p) for g, p in models) == [
        ("Ar", "linear"),
        ("He", "linear"),
        ("Ne", "linear"),
        ("Xe", "circular"),
        ("Xe", "linear"),
    ]
    assert models[Gas.Ar, "linear"] is two_photon_ionisation_model("Ar")


def test_the_manual_table_is_the_models():
    # The manual's table of the models (docs/two_photon.md) states each one's supported
    # and unsupported bands, nm, to 0.01 nm, and its envelope to 0.01.
    page = Path(__file__).parents[1] / "docs" / "two_photon.md"
    text = page.read_text(encoding="utf-8")
    text = text.split("## The model")[1].split("**The regions**")[0]
    rows = {line.split(" | ")[0].strip("| "): line for line in text.splitlines()}
    names = {("Xe", "linear"): "Xe, linear", ("Xe", "circular"): "Xe, circular"}
    for (gas, polarization), model in two_photon_ionisation_models().items():
        row = rows[names.get((gas.value, polarization), gas.value)]
        for region in model.regions:
            if region.label in {"supported", "unsupported_detail"}:
                band = f"{_nm(region.lower):.2f}{DASH}{_nm(region.upper):.2f}"
                assert band in row, (gas, polarization, band)
        assert model.status in row
        envelope = model.uncertainty.envelope
        if envelope is not None:
            assert f"{envelope[0]:.2f}{DASH}{envelope[1]:.2f}" in row


def test_regions_of_a_reference_that_spans_the_whole_window():
    # A reference from the threshold itself to beyond E1 leaves no unavailable band
    # inside the window: the regions run below_threshold, supported, detail, then
    # unavailable from E1.
    found, e1 = derivation.thresholds("He")
    choice = {"wing_ev": 0.6, "resolves_fine_structure": True}
    regions = derivation._regions(found[0].omega, 2 * e1.omega, found, e1, choice)
    assert [r.label for r in regions] == [
        "below_threshold",
        "supported",
        "unsupported_detail",
        "unavailable",
    ]
    assert regions[1].lower == found[0].omega
    assert regions[2].upper == regions[3].lower == e1.omega

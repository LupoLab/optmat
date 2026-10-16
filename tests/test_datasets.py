# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import dataclasses

import numpy as np
import pytest

from optmat.conditions import ONE_ATM_0C
from optmat.constants import AU_POLARIZABILITY
from optmat.datasets import (
    _vacuum_wavelength,
    datasets,
    load_source,
    sources,
)
from optmat.exceptions import DataError
from optmat.gas import reduction
from optmat.gas.eos import number_density
from optmat.materials import materials
from optmat.vocabulary import POLARIZABILITY_QUANTITIES


def test_every_source_loads():
    assert sources()


def test_dataset_keys_are_unique():
    keys = [d.key for d in datasets()]
    assert len(keys) == len(set(keys))


@pytest.mark.parametrize("dataset", datasets(), ids=lambda d: d.key)
def test_dataset_is_finite_and_physical(dataset):
    if dataset.abscissa != "wavelength":
        with pytest.raises(DataError, match="not a wavelength"):
            dataset.vacuum_wavelength()
    else:
        wavelength = dataset.vacuum_wavelength()
        stated = ~np.isnan(wavelength)
        assert np.all((wavelength[stated] > 0) | np.isinf(wavelength[stated]))
    assert np.all(np.isfinite(dataset.scaled_values()))
    if dataset.quantity in POLARIZABILITY_QUANTITIES and not dataset.parameters.get(
        "beyond_resonance", False
    ):
        alpha = reduction.polarizability(dataset) / AU_POLARIZABILITY
        # Between hydrogen's 0.67 a.u. atom and SF6's ~45 a.u. molecule, with room for
        # the rise near a resonance: a gross unit or scale error lands outside. A
        # cation is more tightly bound (He+ is exactly 9/32 = 0.28 a.u.).
        lower = 0.2 if materials()[dataset.material].kind == "ion" else 0.5
        assert np.all((alpha > lower) & (alpha < 400)), dataset.key


@pytest.mark.parametrize("dataset", datasets(), ids=lambda d: d.key)
def test_reproduced_values_are_not_primary(dataset):
    if dataset.reproduced_from:
        assert dataset.use in {"reproduced", "check", "excluded"}


def test_huber_first_value_round_trips_to_refractivity():
    (experiment,) = [d for d in datasets("He") if d.key == "huber_1974/experiment"]
    dn = reduction.refractivity_at(experiment, ONE_ATM_0C)
    # The table's state is the evaluation state, so n - 1 comes back unchanged up to
    # the Lorentz-Lorenz round trip (1e-13) and the unknown-density convention (real).
    np.testing.assert_allclose(dn[0], 4.85e-5, rtol=1e-12)
    assert number_density("He", ONE_ATM_0C) > 0


def test_uncertainty_levels_convert_to_one_sigma():
    (experiment,) = [d for d in datasets("He") if d.key == "huber_1974/experiment"]
    sigma = experiment.standard_uncertainty()
    np.testing.assert_allclose(
        sigma / experiment.scaled_values(), 0.03 * 1.482602218505602, rtol=1e-12
    )


def test_quantity_without_polarizability_raises():
    (theory,) = [d for d in datasets("He") if d.key == "huber_1974/victor_theory"]
    broken = type(theory)(**{**theory.__dict__, "quantity": "oscillator_strength"})
    with pytest.raises(DataError, match="no polarizability of its own"):
        reduction.polarizability(broken)


def test_rydberg_photon_energy_is_half_a_hartree():
    rydberg = _vacuum_wavelength(np.array([2.0]), "rydberg", "vacuum")
    hartree = _vacuum_wavelength(np.array([1.0]), "hartree", "vacuum")
    np.testing.assert_allclose(rydberg, hartree, rtol=1e-15)


def test_excluded_rows_are_dropped_and_recorded(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "polarizability"\nunit = "au"\n'
        'abscissa_unit = "nm"\nexcluded_rows = { "2" = "marked extrapolated" }\n'
        'points = [["1", "500", "11.2"], ["2", "300", "11.5"], ["3", "inf", "11.08"]]\n'
    )
    (dataset,) = load_source(record).datasets
    assert dataset.rows == ("1", "3")
    assert dataset.excluded_rows == {"2": "marked extrapolated"}


def test_excluded_rows_must_exist(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "polarizability"\nunit = "au"\n'
        'abscissa_unit = "nm"\nexcluded_rows = { "9" = "no such row" }\n'
        'points = [["1", "500", "11.2"]]\n'
    )
    with pytest.raises(DataError, match="excluded_rows"):
        load_source(record)


_FIGURE_RECORD = (
    'citation = "demo"\nkind = "theory"\n[[datasets]]\nid = "d"\n'
    'material = "Ar"\nquantity = "polarizability"\nunit = "au"\n'
    'abscissa_unit = "nm"\npoints = [["1", "500", "11.2"]]\n'
)


@pytest.mark.parametrize(
    ("entry", "match"),
    [
        # A figure's values name the figure, state their uncertainty, and an image's
        # are never fitted (TC-699).
        ('digitised = "by eye"\n', "unknown digitisation"),
        (
            'digitised = "vector"\ntable = "tbl_1"\n'
            "[datasets.uncertainty]\nrelative = 0.02\n",
            "must name the library figure",
        ),
        ('digitised = "vector"\ntable = "fig_1"\n', "must state their uncertainty"),
        (
            'digitised = "raster"\ntable = "fig_1"\nuse = "primary"\n'
            "[datasets.uncertainty]\nrelative = 0.05\n",
            "never fitted",
        ),
    ],
)
def test_figure_values_follow_their_rules(tmp_path, entry, match):
    record = tmp_path / "demo_2000.toml"
    record.write_text(_FIGURE_RECORD + entry)
    with pytest.raises(DataError, match=match):
        load_source(record)


def test_a_vector_figure_may_be_checked_or_fitted(tmp_path):
    record = tmp_path / "demo_2000.toml"
    for use in ("check", "primary"):
        record.write_text(
            _FIGURE_RECORD + f'digitised = "vector"\ntable = "fig_1"\nuse = "{use}"\n'
            "[datasets.uncertainty]\nrelative = 0.02\n"
        )
        (dataset,) = load_source(record).datasets
        assert dataset.digitised == "vector"
        assert dataset.use == use


def test_rows_without_a_stated_uncertainty_are_nan_not_zero(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "polarizability"\nunit = "au"\n'
        'abscissa_unit = "nm"\n'
        'points = [["1", "500", "11.2"], ["2", "300", "11.5"]]\n'
        '[datasets.uncertainty]\nrelative_rows = { "1" = 0.01 }\n'
    )
    (dataset,) = load_source(record).datasets
    sigma = dataset.standard_uncertainty()
    # Row 1 has its stated 1 %; row 2 has no statement at all: NaN, not 0 (a 0 would
    # make the value look exact in a fit).
    assert sigma[0] == pytest.approx(0.112)
    assert np.isnan(sigma[1])


def test_a_third_order_dataset_must_state_its_process_and_part():
    template = datasets()[0]
    with pytest.raises(DataError, match="must state parameters"):
        dataclasses.replace(
            template,
            quantity="hyperpolarizability",
            unit="au",
            parameters={},
            state=None,
        )
    with pytest.raises(DataError, match="unknown process"):
        dataclasses.replace(
            template,
            quantity="hyperpolarizability",
            unit="au",
            parameters={"process": "spm", "convention": "T", "part": "electronic"},
        )


def test_the_documented_vocabulary_is_the_modules():
    from pathlib import Path  # noqa: PLC0415

    from optmat.vocabulary import QUANTITIES  # noqa: PLC0415

    page = (
        Path(__file__).parents[1] / "docs" / "the-data" / "quantities.md"
    ).read_text()
    table = page.split("<!-- quantities:start -->")[1]
    table = table.split("<!-- quantities:end -->")[0]
    rows = {
        cells[0].strip("`"): cells
        for line in table.strip().splitlines()[2:]
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }
    assert set(rows) == set(QUANTITIES)
    for name, spec in QUANTITIES.items():
        cells = rows[name]
        assert cells[1] == spec.family
        units = {u.strip("`") for u in cells[3].split(", ")}
        assert units == ({"as printed"} if spec.units is None else set(spec.units))
        assert {a.strip("`") for a in cells[4].split(", ")} == set(spec.abscissae)
        assert (cells[5] == "yes") == spec.needs_state
        required = {r.strip("`") for r in cells[6].split(", ") if r}
        assert required == set(spec.required)
        assert {c.strip("`") for c in cells[7].split(", ") if c} == set(spec.columns)


def test_an_interval_quantity_needs_its_upper_edges(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "theory"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "interval_oscillator_strength"\n'
        'abscissa_unit = "eV"\n'
        'points = [["1", "15.76", "0.05"]]\n'
    )
    with pytest.raises(DataError, match="needs the columns"):
        load_source(record)
    record.write_text(
        'citation = "demo"\nkind = "theory"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "interval_oscillator_strength"\n'
        'abscissa_unit = "eV"\ncolumns = { x_upper = "eV" }\n'
        'points = [["1", "15.76", "0.05", "", "15.94"],\n'
        '          ["2", "15.94", "0.5", "", "inf"]]\n'
    )
    (dataset,) = load_source(record).datasets
    np.testing.assert_array_equal(dataset.column_values("x_upper"), [15.94, np.inf])


def test_a_scale_erratum_corrects_a_misprinted_column_factor(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "polarizability"\nunit = "au"\nscale = 10.0\n'
        'abscissa_unit = "nm"\npoints = [["1", "500", "1.12"]]\n'
        '[[errata]]\ndataset = "d"\nrow = ""\ncolumn = "scale"\nprinted = "10.0"\n'
        'corrected = "1.0"\nreason = "r"\nevidence = "e"\n'
    )
    (dataset,) = load_source(record).datasets
    assert dataset.scaled_values()[0] == pytest.approx(1.12)


def test_a_non_wavelength_abscissa_converts_to_si(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "Ar"\nquantity = "ionization_rate"\nunit = "1/s"\n'
        'abscissa = "intensity"\nabscissa_unit = "W/cm2"\n'
        'points = [["1", "1e14", "3e9"]]\n'
        '[datasets.parameters]\nwavelength_nm = 800.0\npolarization = "linear"\n'
    )
    (dataset,) = load_source(record).datasets
    assert dataset.abscissa_values()[0] == pytest.approx(1e18)
    with pytest.raises(DataError, match="not a wavelength"):
        dataset.vacuum_wavelength()


def test_a_field_abscissa_converts_to_volts_per_metre(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "theory"\n[[datasets]]\nid = "d"\n'
        'material = "He"\nquantity = "ionization_rate"\nunit = "au"\n'
        'abscissa = "field"\nabscissa_unit = "au"\n'
        'points = [["1", "0.1", "1e-6"]]\n'
        '[datasets.parameters]\nwavelength_nm = inf\npolarization = "static"\n'
    )
    (dataset,) = load_source(record).datasets
    # One atomic unit of field is E_h / (e a0) = 5.142206747e11 V/m (CODATA 2018).
    assert dataset.abscissa_values()[0] == pytest.approx(5.142206747e10, rel=1e-9)


def test_an_ion_is_its_own_material_and_not_its_gas(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "theory"\n[[datasets]]\nid = "d"\n'
        'material = "He+"\nquantity = "polarizability"\nunit = "au"\n'
        'abscissa_unit = "nm"\npoints = [["1", "inf", "0.28125"]]\n'
    )
    (dataset,) = load_source(record).datasets
    assert dataset.material == "He+"
    assert materials()["He+"].kind == "ion"
    assert not materials()["He+"].is_solid
    assert all(d.material == "He" for d in datasets("He"))


def test_a_label_abscissa_has_no_values_and_none_has_empty_cells(tmp_path):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        'material = "O3"\nquantity = "ionization_energy"\nunit = "eV"\n'
        'abscissa = "none"\n'
        'points = [["adiabatic", "12", "12.52495"]]\n'
        '[datasets.parameters]\nkind = "adiabatic"\n'
    )
    with pytest.raises(DataError, match="empty cells"):
        load_source(record)


def test_a_named_formula_coefficient_is_corrected_on_loading():
    (entry,) = [
        f for f in sources()["shelton_1990"].formulas if f["id"] == "table2_He_si"
    ]
    coefficients = entry["coefficients"]
    assert isinstance(coefficients, dict)
    assert coefficients["C"] == "6.043e-32"  # printed 6.043e-31


def test_a_wavelength_the_source_does_not_state_is_nan():
    (dataset,) = [
        d for d in datasets("N2") if d.key == "lin_1976/dalpha_hirschfelder_N2"
    ]
    assert np.isnan(dataset.vacuum_wavelength()[0])

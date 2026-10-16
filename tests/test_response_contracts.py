# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Response-model limits and record validation independent of the curated fits."""

import dataclasses

import numpy as np
import pytest

import optmat.gas.raman_response as gas
import optmat.solid.raman_response as solid
from optmat.datasets import load_source
from optmat.exceptions import DataError
from optmat.gas.linewidths import AMAGAT, Broadening, LineTable, LinewidthLaw
from optmat.solid.models.formulas import FormulaModel


@pytest.mark.parametrize("value", [None, [], [1], [1, 2, 3], [0, 1], [2, 1], [-1, 2]])
def test_response_windows_are_positive_increasing_pairs(value):
    with pytest.raises(DataError, match="must"):
        gas._pair(value, "synthetic")
    assert gas._pair([1, 2], "synthetic") == (1.0, 2.0)


def test_linewidth_decomposition_checks_the_fourth_density():
    densities = [*gas._PROBE_DENSITIES, gas._CHECK_DENSITY]
    c, d, b = np.array([2, 3]), np.array([5, 7]), np.array([11, 13])
    widths = tuple(np.asarray(c + d / n + b * n, dtype=float) for n in densities)
    actual = gas._width_coefficients(widths)
    # A four-point independent evaluation of c + d/rho + b*rho; solve rounding only.
    np.testing.assert_allclose(actual, [c, d, b], rtol=1e-12, atol=0)
    bad = (*widths[:3], widths[3] * 1.01)
    with pytest.raises(DataError, match="not of the form"):
        gas._width_coefficients(bad)


def test_formula_metadata_can_include_non_numeric_descriptions(tmp_path, monkeypatch):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[formulas]]\nid="shift"\n'
        'coefficients={C0="2",C1="3",C2="4",T0="300",range="descriptive"}\n'
    )
    source = load_source(path)
    for module in (gas, solid):
        monkeypatch.setattr(module, "sources_of", lambda: {source.id: source})
    assert gas._formula_coefficients("synthetic/shift") == {
        "C0": 2,
        "C1": 3,
        "C2": 4,
        "T0": 300,
    }
    assert solid._formula("synthetic/shift")["range"] == "descriptive"
    assert (
        gas._shift({"law": "polynomial", "formula": "synthetic/shift"}, 302.0) == 24e6
    )
    with pytest.raises(DataError, match="unknown shift law"):
        gas._shift({"law": "missing"}, 302.0)
    with pytest.raises(DataError, match="unknown formula"):
        gas._formula_coefficients("synthetic/missing")
    with pytest.raises(DataError, match="unknown formula"):
        solid._formula("synthetic/missing")
    assert solid._strings([]) == {}
    with pytest.raises(DataError, match="not a time"):
        solid._seconds("10 metres")
    assert solid._seconds("10 fs") == pytest.approx(1e-14, rel=2e-15, abs=0)


def test_unknown_solid_response_forms_and_missing_gas_choices_fail(monkeypatch):
    with pytest.raises(DataError, match="unknown Raman form"):
        solid.SilicaRamanResponse("synthetic", "missing", 0.2, False, {})
    monkeypatch.setattr(gas, "_config", lambda: ({}, {}))
    with pytest.raises(DataError, match="has no choices"):
        gas.raman_response("H2", 300.0)


def test_zero_abundance_needs_no_broadening_coefficient():
    law = LinewidthLaw(broadening={"self": Broadening(2e6, 300)})
    np.testing.assert_array_equal(
        law.fwhm(np.array([0, 2]), 300.0, AMAGAT, {"self": 1.0, "missing": 0.0}),
        [2e6, 2e6],
    )
    with pytest.raises(DataError, match="one coefficient and one exponent"):
        LineTable((0, 2), (1.0,), 300.0, (0.0, 0.0))


def test_solid_index_derivative_follows_a_single_sellmeier_pole():
    model = FormulaModel(
        material="SiO2",
        name="synthetic",
        citation="synthetic",
        wavelength_range=(200e-9, 2e-6),
        kind="published",
        form="n2_sellmeier",
        terms={"B": (1.0,), "C": (0.01,)},
    )
    wavelength = np.array([400e-9, 800e-9])
    # n^2 = 1 + l^2/(l^2-C), with l in um; d(l^2)/domega = -2 l^2/omega.
    from optmat.constants import C  # noqa: PLC0415

    l2 = (wavelength * 1e6) ** 2
    omega = 2 * np.pi * C / wavelength
    index = np.sqrt(1 + l2 / (l2 - 0.01))
    expected = 0.01 * l2 / (omega * (l2 - 0.01) ** 2 * index)
    np.testing.assert_allclose(
        model.dindex_domega(wavelength), expected, rtol=2e-15, atol=0
    )


def test_udm_file_loading_and_missing_parameters_are_explicit(tmp_path):
    from optmat.solid import udm  # noqa: PLC0415

    path = tmp_path / "synthetic.par"
    # A zero-strength interband term gives vacuum, regardless of the band edges.
    path.write_text("ex=0\nEg=8\nEh=20\nNvc=0\nsynthetic description\n")
    model = udm.from_file(path)
    np.testing.assert_array_equal(
        model.permittivity(np.array([1.0, 2.0])), [1 + 0j, 1 + 0j]
    )
    with pytest.raises(DataError, match="UDM parameter Eg is missing"):
        udm.from_parameters((), {})


def test_validation_reports_require_the_named_source_formula(tmp_path, monkeypatch):
    from optmat.gas import ionisation_validation  # noqa: PLC0415
    from optmat.solid import raman_response_validation  # noqa: PLC0415

    path = tmp_path / "synthetic.toml"
    path.write_text('citation="synthetic"\nkind="theory"\n')
    source = load_source(path)
    monkeypatch.setattr(ionisation_validation, "sources", lambda: {"synthetic": source})
    with pytest.raises(DataError, match="has no formula"):
        ionisation_validation._coefficients("synthetic", "missing")
    source = dataclasses.replace(
        source,
        formulas=(
            {"id": "unrelated", "coefficients": {"x": "1"}},
            {"id": "target", "coefficients": {"x": "2"}},
        ),
    )
    assert ionisation_validation._coefficients("synthetic", "target") == {"x": 2.0}
    malformed = dataclasses.replace(
        source, formulas=({"id": "low_frequency_gain", "coefficients": []},)
    )
    monkeypatch.setattr(
        raman_response_validation, "sources_of", lambda: {"stolen_1989": malformed}
    )
    with pytest.raises(DataError, match="has no coefficients"):
        raman_response_validation.silica_low_frequency_gain()

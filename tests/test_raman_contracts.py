# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Synthetic spectra protect Raman record assembly and linewidth conventions."""

import dataclasses
import tomllib

import numpy as np
import pytest

from optmat.datasets import load_source
from optmat.exceptions import DataError
from optmat.gas import raman
from optmat.gas.linewidths import AMAGAT, ExponentialGapLaw


@pytest.fixture
def source(tmp_path, monkeypatch):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n'
        '[[datasets]]\nid="width"\nmaterial="H2"\nquantity="raman_linewidth"\n'
        'unit="MHz/amagat"\nabscissa="label"\nabscissa_unit=""\n'
        'points=[["a","0","2"],["b","2","3"]]\n'
        'parameters={transition="S",width="hwhm",temperature_k=300}\n'
        '[[formulas]]\nid="law"\n[formulas.coefficients]\n'
        'constant_MHz="2"\nB_MHz_per_amagat="3"\n'
        'alpha="0.02"\nbeta="1.5"\ndelta="1.2"\na="1.0"\n'
        'n_experiment="1.3"\nT0_experiment="300"\nrange="descriptive text"\n',
        encoding="utf-8",
    )
    source = load_source(path)
    monkeypatch.setattr(raman, "sources_of", lambda: {source.id: source})
    return source


@pytest.mark.parametrize(
    ("unit", "factor"),
    [
        ("MHz/amagat", 2e6),
        ("GHz/amagat", 2e9),
        ("cm-1/amagat", 2 * 100 * 299792458),
        ("cm-1/atm", 2 * 100 * 299792458 * 300 / 273.15),
    ],
)
def test_half_width_units_convert_to_full_width_per_density(source, unit, factor):
    table = dataclasses.replace(source.datasets[0], unit=unit)
    width, temperature = raman._per_amagat(table)
    # Each synthetic HWHM is doubled; pressure coefficients additionally use N=p/kT.
    np.testing.assert_allclose(width, np.array([2, 3]) * factor, rtol=2e-15, atol=0)
    assert temperature == 300


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"unit": "Hz/Pa"}, "no conversion"),
        (
            {"parameters": {"transition": "S", "width": "fwhm"}},
            "width convention and temperature",
        ),
        (
            {
                "parameters": {
                    "transition": "S",
                    "width": "not stated",
                    "temperature_k": 300,
                }
            },
            "width convention and temperature",
        ),
    ],
)
def test_widths_require_supported_units_convention_and_temperature(
    source, changes, message
):
    with pytest.raises(DataError, match=message):
        raman._per_amagat(dataclasses.replace(source.datasets[0], **changes))


def test_constant_plus_collisional_width_retains_its_density_independent_term(source):
    entry = {
        "type": "law",
        "law": "constant_linear",
        "formula": "synthetic/law",
        "substitutes": {"He": "self"},
    }
    law = raman._widths({}, entry, None)
    width = law.fwhm(np.array([0, 2]), 300.0, 2 * AMAGAT, {"He": 1.0})
    np.testing.assert_array_equal(width, [8e6, 8e6])
    assert law.sources == ("synthetic/law",)
    assert "range" not in raman._coefficients("synthetic/law")


def test_exponential_gap_configuration_preserves_constants_and_substitutions(source):
    constants = raman.RotationalConstants(2.0, 1e-5)
    entry = {
        "type": "exponential_gap",
        "formula": "synthetic/law",
        "substitutes": {"He": "self"},
    }
    law = raman._widths({}, entry, constants)
    assert isinstance(law, ExponentialGapLaw)
    assert (law.b, law.d, law.alpha, law.beta) == (2.0, 1e-5, 0.02, 1.5)
    assert law.foreign_as_self == ("He",)
    # The law is linear in density; a declared substitute has exactly the self width.
    j = np.array([0, 2, 4])
    pure = law.fwhm(j, 300.0, AMAGAT)
    np.testing.assert_allclose(
        law.fwhm(j, 300.0, 2 * AMAGAT, {"He": 1.0}), 2 * pure, rtol=2e-15, atol=0
    )
    with pytest.raises(DataError, match="needs the band's rotational constants"):
        raman._widths({}, entry, None)


def test_ground_state_rotational_constants_are_recovered_from_both_row_orientations(
    source,
):
    # Synthetic E(J) = B J(J+1) - D [J(J+1)]^2; lines span several J to separate B,D.
    b, d = 2.0, 0.003
    lows = np.array([0, 1, 2, 3])
    lo, hi = lows * (lows + 1), (lows + 2) * (lows + 3)
    shifts = b * (hi - lo) - d * (hi**2 - lo**2)
    table = dataclasses.replace(
        source.datasets[0],
        quantity="raman_shift",
        unit="cm-1",
        rows=("2-0", "1-3", "4-2", "3-5"),
        abscissae=("",) * 4,
        values=tuple(str(x) for x in shifts),
        sigmas=("",) * 4,
        parameters={},
    )
    recovered = raman._constants({table.key: table}, {"fit_lines": table.key})
    # Small well-conditioned linear least-squares problem: rounding only.
    assert recovered.b == pytest.approx(b, rel=1e-12)
    assert recovered.d == pytest.approx(d, rel=1e-12)
    assert recovered.h == 0


@pytest.mark.parametrize(
    ("function", "arguments", "message"),
    [
        (raman._table, ([], "test"), "must be a table"),
        (raman._list, ({}, "test"), "must be an array"),
        (raman._get, ({}, "missing"), "unknown dataset"),
        (raman._level, ("unlabelled",), "names no rotational level"),
        (raman._widths, ({}, {"type": "missing"}, None), "unknown width type"),
    ],
)
def test_malformed_raman_configuration_is_rejected(function, arguments, message):
    with pytest.raises(DataError, match=message):
        function(*arguments)


def test_raman_references_require_a_unique_row_and_known_formula(source):
    table = source.datasets[0]
    by_key = {table.key: table}
    with pytest.raises(DataError, match="has no row"):
        raman._get(by_key, table.key + "#missing")
    with pytest.raises(DataError, match="name one as"):
        raman._scalar(by_key, table.key)
    assert raman._scalar(by_key, table.key + "#b") == 3
    with pytest.raises(DataError, match="not a wavenumber"):
        raman._wavenumber(by_key, table.key + "#b")
    with pytest.raises(DataError, match="not a mean derivative"):
        raman._transition_from_derivative(by_key, table.key, 1.0)
    with pytest.raises(DataError, match="unknown formula"):
        raman._coefficients("synthetic/missing")
    with pytest.raises(DataError, match="unknown width law"):
        raman._law({}, {"law": "missing", "formula": "synthetic/law"})
    with pytest.raises(DataError, match="Dicke coefficient"):
        raman._widths(
            by_key, {"type": "table", "perturbers": {}, "dicke": table.key}, None
        )


def test_a_rotational_nonlinear_index_requires_linewidths():
    lines = dataclasses.replace(raman.raman_lines("H2", 300.0, AMAGAT), fwhm=None)
    with pytest.raises(DataError, match="needs the lines' widths"):
        raman.rotational_nonlinear_index(lines, AMAGAT)


def test_dicke_plus_linear_law_and_scalar_transition_assembly(source, monkeypatch):
    import optmat.gas.raman_response as response  # noqa: PLC0415
    from optmat.constants import AU_POLARIZABILITY  # noqa: PLC0415

    formula = dict(source.formulas[0])
    formula["coefficients"] = dict(formula["coefficients"]) | {"A_MHz_amagat": "4"}
    source = dataclasses.replace(source, formulas=(formula,))
    monkeypatch.setattr(raman, "sources_of", lambda: {source.id: source})
    law = raman._law({}, {"law": "dicke_linear", "formula": "synthetic/law"})
    assert law.fwhm(np.array([0]), 300.0, 2 * AMAGAT)[0] == 8e6  # 4/rho + 3*rho MHz
    base = source.datasets[0]
    origin = dataclasses.replace(
        base,
        id="origin",
        quantity="raman_shift",
        unit="cm-1",
        rows=("0",),
        abscissae=("",),
        values=("1000",),
        sigmas=("",),
        parameters={},
    )
    transition = dataclasses.replace(
        origin,
        id="transition",
        quantity="raman_transition_polarizability",
        unit="au",
        abscissa="none",
        values=("2",),
        parameters={"transition": "Q", "invariant": "mean"},
    )
    by_key = {d.key: d for d in (origin, transition)}
    band = raman._vibrational(
        by_key,
        {
            "mode": "synthetic",
            "band_origin": origin.key,
            "transition_polarizability": transition.key,
            "widths": {
                "type": "law",
                "law": "dicke_linear",
                "formula": "synthetic/law",
            },
        },
        None,
    )
    assert band.band_origin == 1000
    assert band.transition_polarizability == 2 * AU_POLARIZABILITY
    np.testing.assert_array_equal(
        response._mean_element(band, np.array([0, 2])),
        np.full(2, 2 * AU_POLARIZABILITY),
    )
    assert band.sources["transition_polarizability"] == transition.key
    with pytest.raises(DataError, match="rotational band needs"):
        raman._parameters(raman.Gas.H2, {"rotational": {}}, {})
    with monkeypatch.context() as isolated:
        isolated.setattr(raman, "_all_parameters", lambda: {})
        with pytest.raises(DataError, match="no Raman parameters"):
            raman.raman_parameters("Ar")
    other = dataclasses.replace(base, id="other", rows=("c", "d"))
    with pytest.raises(DataError, match="cover different levels"):
        raman._components(
            {base.key: base, other.key: other},
            {"parallel": base.key, "perpendicular": other.key},
        )


def test_explicit_q_branch_cutoff_preserves_a_constant_transition_frequency():
    constants = raman.RotationalConstants(2.0)
    lines = raman.vibrational_lines(
        1000.0,
        1e-41,
        temperature=300.0,
        lower=constants,
        upper=constants,
        spin=raman.spin_statistics("H2"),
        j_max=100,
    )
    # Identical upper/lower rotational ladders give no rotational splitting.
    np.testing.assert_allclose(
        lines.frequency, 2 * np.pi * 299792458 * 100 * 1000, rtol=2e-15, atol=0
    )
    assert np.all(lines.amplitude > 0)
    assert np.max(lines.j) <= 100


def test_disjoint_temperature_tables_do_not_invent_a_temperature_exponent(source):
    first = dataclasses.replace(source.datasets[0], abscissae=("J=0", "J=2"))
    second = dataclasses.replace(
        first,
        id="hot",
        abscissae=("J=4", "J=6"),
        parameters=first.parameters | {"temperature_k": 600},
    )
    law, used = raman._line_table(
        {first.key: first, second.key: second},
        {"dataset": first.key, "temperatures": [second.key]},
    )
    # No level has two temperatures: exponents remain zero, and hot-only levels use
    # their measured widths at the reference temperature without an invented slope.
    assert law.exponent == (0, 0, 0, 0)
    assert law.j == (0, 2, 4, 6)
    # exp(log(width)) introduces floating-point rounding.
    assert law.coefficient == pytest.approx((4e6, 6e6, 4e6, 6e6), rel=2e-15)
    assert used == [first.key, second.key]


@pytest.mark.parametrize(
    "resolution,use_lower,message",
    [
        ({}, True, r"resolved.*upper.*fit_lines.*b_difference"),
        ({"resolved": True}, True, r"resolved.*upper.*fit_lines.*b_difference"),
        ({"resolved": True}, False, r"resolved.*lower"),
    ],
)
def test_resolved_band_cannot_silently_drop_missing_ladders(
    resolution, use_lower, message
):
    parameters = raman.raman_parameters("N2")
    by_key = {
        d.key: d for source in raman.sources_of().values() for d in source.datasets
    }
    with raman.RAMAN_FILE.open("rb") as handle:
        original = tomllib.load(handle)["N2"]["vibrational"][0]
    assert parameters.rotational is not None
    lower = parameters.rotational.constants
    valid = raman._vibrational(by_key, original, lower)
    assert valid.lower == lower
    assert valid.upper is not None
    # Remove the upper ladder only in a local copy of an otherwise valid record.
    entry = {
        k: v
        for k, v in original.items()
        if k not in {"upper", "fit_lines", "b_difference"}
    }
    entry.update(resolution)
    with pytest.raises(DataError, match=message):
        raman._vibrational(by_key, entry, lower if use_lower else None)
    unresolved = raman._vibrational(by_key, entry | {"resolved": False}, lower)
    assert unresolved.lower is None
    assert unresolved.upper is None
    assert unresolved.band_origin == valid.band_origin
    assert unresolved.transition_polarizability == valid.transition_polarizability

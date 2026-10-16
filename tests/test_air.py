# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Air as the sum of its constituents, its compositions and its validator (OEP 1, §4.4).

Each air model is its composition's weighted sum of the constituents' recommended
models, value and derivative; the reduced compositions are CIPM-2007's, renormalised;
the full composition agrees with Ciddor's standard dry air within his stated
uncertainty; the pseudo-pure air equation lies above the partial-pressure sum by the
amounts OEP 1 states.
"""

import math
import warnings

import numpy as np
import pytest

from optmat import _validity
from optmat.conditions import State
from optmat.constants import AU_ANGULAR_FREQUENCY, C
from optmat.datasets import sources
from optmat.exceptions import DataError, ModelRangeWarning, OptmatError
from optmat.gas import air, raman_response
from optmat.gas.air import AIR_COMPOSITIONS, DEFAULT_AIR, air_model
from optmat.gas.composition import composition
from optmat.gas.eos import air_equation_of_state, number_density
from optmat.gas.infrared import gap_cells
from optmat.gas.models.registry import default_model

WAVELENGTHS = np.array([400e-9, 633e-9, 800e-9, 1.2e-6])


@pytest.mark.parametrize("name", AIR_COMPOSITIONS)
def test_an_air_model_is_the_sum_of_its_constituents(name):
    model = air_model(name)
    assert math.isclose(sum(c.mole_fraction for c in model.constituents), 1.0)
    for method in ("complex_polarizability", "dpolarizability_domega"):
        expected = sum(
            c.mole_fraction * getattr(default_model(c.material), method)(WAVELENGTHS)
            for c in model.constituents
        )
        got = getattr(model, method)(WAVELENGTHS)
        # The same terms in the same order: the sum agrees to its own rounding.
        np.testing.assert_allclose(got, expected, rtol=1e-14)
    assert np.all(
        model.polarizability(WAVELENGTHS)
        == np.real(model.complex_polarizability(WAVELENGTHS))
    )


def test_the_reduced_compositions_are_cipm_2007_renormalised():
    full = {c.material: c.mole_fraction for c in composition("dry air (CIPM-2007)")}
    for name, kept in (
        ("dry air (N2, O2, Ar)", ("N2", "O2", "Ar")),
        ("dry air (N2, O2)", ("N2", "O2")),
    ):
        reduced = {c.material: c.mole_fraction for c in composition(name)}
        total = sum(full[s] for s in kept)
        assert set(reduced) == set(kept)
        for s in kept:
            assert reduced[s] == pytest.approx(full[s] / total, rel=1e-15)


def test_the_reduced_compositions_against_solprops_recipes():
    # solprop's air is 0.7808, 0.2095, 0.0093 (and 0.7808, 0.2095 without argon),
    # renormalised. CIPM-2007's O2, 0.20939, lies 1.1e-4 below solprop's 0.2095, beyond
    # its four-digit rounding (TC-759); the others agree within 7.3e-5.
    recipes = {
        "dry air (N2, O2, Ar)": {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093},
        "dry air (N2, O2)": {"N2": 0.7808, "O2": 0.2095},
    }
    for name, recipe in recipes.items():
        total = sum(recipe.values())
        ours = {c.material: c.mole_fraction for c in composition(name)}
        for s, x in recipe.items():
            assert abs(ours[s] - x / total) < 1.2e-4


def test_the_full_composition_agrees_with_ciddor():
    # Ciddor (1996), Eqs. (1)-(2): standard dry air at 15 C, 101325 Pa, with CIPM-2007's
    # 400 umol/mol of CO2. His stated uncertainty is 2-5e-8 in n (p. 1568); the model's
    # own (the constituents' fits, a few 1e-10) and the unmodelled pair terms (about
    # 2e-9, N2's B_R rho / A_R of 8e-6) are far smaller.
    (formula,) = sources()["ciddor_1996"].formulas
    coefficients = formula["coefficients"]
    assert isinstance(coefficients, dict)
    k = {key: float(v) for key, v in coefficients.items()}
    model = air_model("dry air (CIPM-2007)")
    low, high = model.intervals[0]
    wavelength = np.geomspace(low, high, 20)
    s = 1e-6 / wavelength
    ciddor = 1e-8 * (k["k1"] / (k["k0"] - s**2) + k["k3"] / (k["k2"] - s**2))
    ciddor = ciddor * (1 + k["co2"] * (400 - 450))
    ours = model.refractivity(wavelength, State(101325.0, 288.15))
    assert np.all(np.abs(ours - ciddor) < 5e-8)


def test_the_omission_of_the_reduced_compositions():
    # N2/O2/Ar air leaves out CIPM-2007's CO2 and traces: its polarizability is about
    # 1.9e-4 below the full composition's; N2/O2 air, without argon, 1.5-1.7e-4 above.
    full = air_model("dry air (CIPM-2007)").polarizability(WAVELENGTHS)
    three = air_model("dry air (N2, O2, Ar)").polarizability(WAVELENGTHS)
    two = air_model("dry air (N2, O2)").polarizability(WAVELENGTHS)
    assert np.all((three / full - 1 > -1.96e-4) & (three / full - 1 < -1.87e-4))
    assert np.all((two / full - 1 > 1.49e-4) & (two / full - 1 < 1.87e-4))


@pytest.mark.parametrize(
    ("pressure", "expected"), [(1e5, 1.9e-4), (1e6, 1.8e-3), (5e6, 6.0e-3)]
)
def test_the_air_equation_holds_the_unlike_pairs(pressure, expected):
    # The pseudo-pure dry air (Lemmon et al. 2000) against the constituents' densities
    # at their partial pressures (solprop's mixture rule), N2/O2/Ar at 294 K.
    state = State(pressure, 294.0)
    pure = air_equation_of_state().number_density(state)
    partial = sum(
        number_density(c.material, State(pressure * c.mole_fraction, 294.0))
        for c in composition(DEFAULT_AIR)
    )
    assert 1 - partial / pure == pytest.approx(expected, rel=0.05)


def test_the_intervals_are_the_constituents_intersection_widened():
    full = air_model("dry air (CIPM-2007)")
    sets = [default_model(c.material).intervals for c in full.constituents]
    base = air._intersection(sets)
    lows = [s[0][0] for s in sets]
    highs = [s[-1][1] for s in sets]
    assert base[0][0] == max(lows)
    assert base[-1][1] <= min(highs)
    # The trace CH4 and N2O split the intersection where they absorb as pure gases
    # (OEP 3), while the reduced air has one interval. Every base interval lies inside
    # a widened one (OEP 3, Q13).
    assert len(base) > 1
    assert len(air_model().intervals) == 1
    assert air_model().intervals == air._intersection(
        [default_model(c.material).intervals for c in air_model().constituents]
    )
    for a, b in base:
        assert any(c <= a and b <= d for c, d in full.intervals)
    assert full.intervals != base
    assert set(full.identity_record()["extras"]["gap_tables"]) == {
        "CO",
        "CO2",
        "CH4",
        "N2O",
    }
    with pytest.warns(ModelRangeWarning, match="air model"):
        full.polarizability(max(lows) / 2)
    with pytest.raises(OptmatError, match="share no wavelength"):
        air._intersection([((1.0, 2.0),), ((3.0, 4.0),)])


def test_the_default_and_the_names():
    assert air_model("dry air") is air_model()
    assert air_model().name == DEFAULT_AIR
    with pytest.raises(OptmatError, match="no air model"):
        air_model("moist air")


def test_the_identity_binds_the_constituents():
    model = air_model()
    record = model.identity_record()
    digests = [m["digest"] for m in record["fields"]["models"]]
    assert digests == [m.identity.digest for m in model.models]
    assert model.identity.name == "air/dry air (N2, O2, Ar)"
    assert model.identity.sources == ("picard_2008",)
    assert "4.27e-04" in model.status_reason


def test_the_index_at_a_state():
    model = air_model()
    state = State(1e5, 294.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        real = model.refractivity(800e-9, state)
        ideal = model.refractivity(800e-9, state, ideal=True)
        n = model.refractive_index(800e-9, state)
    assert float(np.real(n)) - 1 == pytest.approx(float(real), rel=1e-12)
    # Air is attractive at 294 K: the real gas is denser than the ideal, by about 4e-4.
    assert 0 < float(real / ideal) - 1 < 1e-3


class TestTheRamanResponse:
    def test_dry_air_is_the_first_reduced_composition(self):
        a = raman_response.raman_response("dry air", 295.0)
        b = raman_response.raman_response("dry air (N2, O2, Ar)", 295.0)
        assert a.identity.digest == b.identity.digest
        two = raman_response.raman_response("dry air (N2, O2)", 295.0)
        assert dict(two.composition) == pytest.approx(
            {c.material: c.mole_fraction for c in composition("dry air (N2, O2)")}
        )

    def test_the_full_composition_is_refused(self):
        with pytest.raises(DataError, match="no Raman response"):
            raman_response.raman_response("dry air (CIPM-2007)", 295.0)

    def test_a_mixture_entry_must_be_its_composition(self):
        mixtures = {"dry air (N2, O2)": {"species": ["N2", "O2", "Ar"], "reason": ""}}
        with pytest.raises(DataError, match="differs from its composition"):
            raman_response._composition("dry air (N2, O2)", mixtures)


def test_no_constituent_pole_lies_inside_the_widened_intervals():
    # The analytic check of OEP 3, Q13, independent of the cells.
    full = air_model("dry air (CIPM-2007)")
    for model in full.models:
        energies = np.asarray(getattr(model, "energies", ()))
        wavelength = 2 * np.pi * C / (energies * AU_ANGULAR_FREQUENCY)
        assert np.all(_validity.outside(wavelength, full.intervals)), model.name


def test_the_widened_pieces_pass_the_mixtures_checks():
    full = air_model("dry air (CIPM-2007)")
    base = air._intersection([m.intervals for m in full.models])
    fractions = np.array([c.mole_fraction for c in full.constituents])
    for low, high in full.intervals:
        if any(a <= low and high <= b for a, b in base):
            continue
        wavelength = np.geomspace(low, high, 41)[1:-1]
        wavenumber = 1e-2 / wavelength
        absorption = np.zeros(wavelength.size)
        for x, model in zip(fractions, full.models, strict=True):
            unqualified = _validity.outside(wavelength, model.intervals)
            if not np.any(unqualified):
                continue
            cells = gap_cells(model.gas)
            assert cells is not None
            index = np.searchsorted(cells["high_cm"], wavenumber)
            assert not np.any(cells["pole"][index][unqualified])
            absorption += np.where(unqualified, x * cells["absorption_m"][index], 0.0)
        assert np.all(absorption < air.ABSORPTION_THRESHOLD_M)
        assert 1e-2 / low - 1e-2 / high >= air.MIN_WIDTH_CM


def test_the_widening_needs_each_unqualified_constituents_table(monkeypatch):
    full = air_model("dry air (CIPM-2007)")
    base = air._intersection([m.intervals for m in full.models])
    # Without CH4's table, nothing beyond its intervals can qualify.
    real = air.gap_cells

    def without_ch4(gas):
        return None if str(gas) == "CH4" else real(gas)

    monkeypatch.setattr(air, "gap_cells", without_ch4)
    narrowed = air._widened(full.constituents, full.models, base)
    ch4 = full.models[[c.material for c in full.constituents].index("CH4")]
    for low, high in narrowed:
        assert not np.any(_validity.outside(np.geomspace(low, high, 50), ch4.intervals))
    # A tolerance nothing meets keeps the intersection alone.
    monkeypatch.setattr(air, "gap_cells", real)
    monkeypatch.setattr(air, "ALPHA_FRACTION", 0.0)
    assert air._widened(full.constituents, full.models, base) == base
    # An absorbing trace leaves nothing to test against the tolerances.
    monkeypatch.setattr(air, "ALPHA_FRACTION", 0.1)

    def opaque(gas):
        table = real(gas)
        if table is None:
            return None
        return table | {"absorption_m": table["absorption_m"] + 1e9}

    monkeypatch.setattr(air, "gap_cells", opaque)
    assert air._widened(full.constituents, full.models, base) == base
    # The merge drops a new piece narrower than the minimum width.
    assert air._merge([(2e-6, 2.001e-6)], ((1e-6, 1.5e-6),)) == ((1e-6, 1.5e-6),)
    assert air._merge([(1.5e-6, 1.6e-6)], ((1e-6, 1.5e-6),)) == ((1e-6, 1.6e-6),)

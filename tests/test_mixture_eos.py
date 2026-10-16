# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The equations of state of gas mixtures (OEP 4, Phase 2; its sections 4.2, 7.4 and 8).

Every method preserves the composition and round-trips pressure to density to pressure;
multi-fluid reduces to the pure equation at x_s = 1; the fallback reduces to the pure
equations and its pair term is the one section 4.2 writes; the pairs are checked against
CoolProp's own models; dry air agrees with CIPM-2007's compressibility.
"""

import math
import warnings

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import K_B, N_A
from optmat.exceptions import DataError, ModelRangeWarning, OptmatError, StateError
from optmat.gas import mixture_eos
from optmat.gas.eos import air_equation_of_state, equation_of_state
from optmat.gas.mixture_eos import (
    METHODS,
    TRUNCATION_BOUND,
    mixture_equation_of_state,
)

AMBIENT = State(pressure=1e6, temperature=294.0)
"""10 bar, 294 K: the cross terms are 1e-3 of Z, far above the solves' tolerance."""


def _quiet(composition, **kwargs):
    """Build an equation of state, silencing the provisional warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        return mixture_equation_of_state(composition, **kwargs)


@pytest.mark.parametrize("gas", ["He", "Ne", "Ar", "N2", "CO2"])
def test_one_species_is_the_pure_equation(gas):
    eos = mixture_equation_of_state({gas: 1.0})
    assert (eos.method, eos.qualified, eos.pairs) == ("multi-fluid", True, ())
    pure = equation_of_state(gas)
    state = State(pressure=2e6, temperature=320.0)
    # The same CoolProp equation, one fluid: the same Z to the flash's rounding.
    assert eos.number_density(state) == pytest.approx(
        pure.number_density(state), rel=1e-13
    )
    for method in ("partial-density",):
        fallback = mixture_equation_of_state({gas: 1.0}, method=method)
        # The fallback at x = 1 is the pure equation's direct (rho, T) evaluation.
        assert fallback.number_density(state) == pytest.approx(
            pure.number_density(state), rel=1e-10
        )


@pytest.mark.parametrize("method", METHODS[::2])
def test_the_composition_is_preserved(method):
    eos = _quiet({"He": 0.3, "Ar": 0.7}, method=method)
    partial = eos.partial_densities(AMBIENT)
    total = eos.number_density(AMBIENT)
    assert partial == {"He": 0.3 * total, "Ar": 0.7 * total}


@pytest.mark.parametrize("method", METHODS[::2])
def test_the_round_trip(method):
    eos = _quiet({"Ne": 0.5, "Ar": 0.5}, method=method)
    # Within the fallback's truncation bound, which Ne/Ar reaches near 45 bar.
    for pressure in (1e5, 1e6, 3e6):
        density = eos.number_density(State(pressure=pressure, temperature=300.0))
        back = float(eos.pressure(density, 300.0))
        # brentq to 4 eps on the forward law: the round trip to a few ulps.
        assert back == pytest.approx(pressure, rel=1e-12)
    out = eos.pressure(np.array([1e25, 2e25]), np.array([[280.0], [300.0]]))
    assert out.shape == (2, 2)


def test_the_fallback_is_section_4_2s():
    eos = _quiet({"He": 0.5, "Ar": 0.5})
    # The default; provisional, since it misses Blancett et al.'s Z (TC-779).
    assert (eos.method, eos.qualified) == ("virial", False)
    density = 2.4e26
    temperature = 294.0
    pure = sum(
        float(equation_of_state(g).pressure(0.5 * density, temperature))
        for g in ("He", "Ar")
    )
    # B_12 from Brewer and Vaughn's E with the reference equations' pure B's, linear
    # in T between their 273.15 and 298.15 K points.
    e = np.interp(temperature, [273.15, 298.15], [21.89e-6, 20.13e-6])
    b_he = mixture_eos._pure_second_virial("Helium", temperature)
    b_ar = mixture_eos._pure_second_virial("Argon", temperature)
    cross = e + (b_he + b_ar) / 2
    rho = density / N_A
    expected = pure + K_B * N_A * temperature * rho**2 * 2 * 0.25 * cross
    # Each pure equation at its partial density (its direct (rho, T) evaluation
    # against the inverse's solve, 1e-10), and the cross term of section 4.2.
    assert float(eos.pressure(density, temperature)) == pytest.approx(
        expected, rel=1e-10
    )


def test_the_partial_density_rule_differs_from_multi_fluid():
    # He-N2 has a fitted multi-fluid model (Kunz and Wagner's) and no curated cross
    # virial, so its partial-density rule has no cross term.
    multi = _quiet({"He": 0.5, "N2": 0.5}, method="multi-fluid")
    partial = _quiet({"He": 0.5, "N2": 0.5}, method="partial-density")
    state = State(pressure=1e5, temperature=294.0)
    ratio = multi.number_density(state) / partial.number_density(state)
    # Section 4.2: rho/rho_pp = 1 - rho sum_{s != t} x_s x_t B_st, with B_st the
    # multi-fluid model's effective cross coefficient.
    b_mix = multi.second_virial(294.0)
    b_like = partial.second_virial(294.0)
    rho = partial.number_density(state) / N_A
    expected = 1 - rho * (b_mix - b_like)
    # At 1 bar the cross term is 1.5e-4; the neglected ones, the third virials' and
    # (rho B)^2, are about 2e-6 (1e-4 at 10 bar).
    assert ratio == pytest.approx(expected, abs=5e-6)
    assert ratio < 1 - 1e-4


def test_the_second_virial_is_the_low_density_limit():
    eos = _quiet({"He": 0.5, "Ar": 0.5}, method="multi-fluid")
    state = State(pressure=100.0, temperature=294.0)
    z = eos.compressibility(state)
    rho = eos.number_density(state) / N_A
    # Z - 1 = B rho + C rho^2: at 100 Pa C rho is about 2e-5 of B; the rounding of
    # Z - 1 (1e-7) about 1e-9.
    assert (z - 1) / rho == pytest.approx(eos.second_virial(294.0), rel=5e-5)


def test_dry_air_against_cipm_2007():
    # CIPM-2007's Z (picard_2008, Eq. A1.4) for dry air (x_v = 0), 600-1100 hPa and
    # 15-27 C, u(Z) = 1.5e-5. The reduced composition leaves out CIPM-2007's CO2 and
    # traces, which move Z by about 1.4e-6 (computed with CO2 added).
    a0, a1, a2, d = 1.58123e-6, -2.9331e-8, 1.1043e-10, 1.83e-11
    eos = _quiet("dry air (N2, O2, Ar)", method="multi-fluid")
    worst = 0.0
    for pressure in (6e4, 8e4, 101325.0, 1.1e5):
        for celsius in (15.0, 20.0, 27.0):
            temperature = celsius + 273.15
            cipm = (
                1
                - pressure / temperature * (a0 + a1 * celsius + a2 * celsius**2)
                + (pressure / temperature) ** 2 * d
            )
            got = eos.compressibility(State(pressure=pressure, temperature=temperature))
            worst = max(worst, abs(got - cipm))
    assert worst < 2 * math.hypot(1.5e-5, 1.4e-6)
    assert worst < 1e-5


def test_dry_air_against_the_pseudo_pure_fluid():
    eos = _quiet("dry air (N2, O2, Ar)", method="multi-fluid")
    state = State(pressure=1e6, temperature=294.0)
    ratio = eos.number_density(state) / air_equation_of_state().number_density(state)
    # OEP 4, section 7.4: they differ by 2.0e-4 at 10 bar.
    assert ratio - 1 == pytest.approx(-2.0e-4, abs=0.05e-4)


class TestTheMolarMass:
    def test_one_species_is_the_pure_equations(self):
        assert (
            mixture_equation_of_state({"Ar": 1.0}).molar_mass
            == equation_of_state("Ar").molar_mass
        )

    def test_it_is_the_mole_fractions_mean(self):
        m = {g: equation_of_state(g).molar_mass for g in ("He", "Ar")}
        # Two products and a sum: rounding only.
        assert _quiet({"He": 0.3, "Ar": 0.7}).molar_mass == pytest.approx(
            0.3 * m["He"] + 0.7 * m["Ar"], rel=1e-15
        )

    def test_dry_air_is_its_compositions_not_the_pseudo_pure_fluids(self):
        eos = _quiet("dry air (N2, O2, Ar)")
        mean = math.fsum(
            x * equation_of_state(s).molar_mass for s, x in eos.composition
        )
        assert eos.molar_mass == mean
        # 28.9597 g/mol against CoolProp's pseudo-pure 28.96546 (TC-832).
        assert eos.molar_mass != air_equation_of_state().molar_mass


def test_the_pairs_are_coolprops():
    eos = _quiet({"He": 0.4, "Ne": 0.3, "Ar": 0.3}, method="multi-fluid")
    assert [(p.species, p.partner) for p in eos.pairs] == [
        ("Ar", "He"),
        ("Ar", "Ne"),
        ("He", "Ne"),
    ]
    assert {p.model for p in eos.pairs} == {"Tkaczuk-JPCRD-2020"}
    assert all(p.kind == "fitted" and p.curated for p in eos.pairs)
    assert eos.identity.sources == ("brewer_1969", "tkaczuk_2020")


def test_a_placeholder_and_a_missing_pair():
    eos = _quiet({"He": 0.5, "O2": 0.5})
    (pair,) = eos.pairs
    assert (pair.kind, pair.curated, pair.qualified) == ("placeholder", False, False)
    assert eos.method == "partial-density"
    assert eos.status == "provisional"
    with pytest.raises(OptmatError, match="no fitted multi-fluid model for He-O2"):
        mixture_equation_of_state({"He": 0.5, "O2": 0.5}, method="multi-fluid")
    with pytest.raises(OptmatError, match=r"Ar-Xe \(none\)"):
        mixture_equation_of_state({"Ar": 0.5, "Xe": 0.5}, method="multi-fluid")
    with pytest.raises(OptmatError, match="no curated cross virial for Ar-Xe"):
        mixture_equation_of_state({"Ar": 0.5, "Xe": 0.5}, method="virial")
    fitted = _quiet({"He": 0.5, "N2": 0.5}, method="multi-fluid")
    assert (
        fitted.pairs[0].reason == "CoolProp's binary model is not curated or compared"
    )
    assert fitted.species == ("He", "N2")


def test_a_provisional_equation_warns_when_built():
    with pytest.warns(ModelRangeWarning, match="provisional"):
        mixture_equation_of_state({"Ar": 0.5, "Xe": 0.5})


@pytest.mark.parametrize(
    ("composition", "match"),
    [
        ({}, "distinct species"),
        ({"Ar": 0.5, "He": 0.6}, "sum to"),
        ({"Ar": 1.0, "He": 0.0}, r"in \(0, 1\]"),
        ("wet air", "no mixture composition named"),
    ],
)
def test_the_composition_is_checked(composition, match):
    with pytest.raises(OptmatError, match=match):
        mixture_equation_of_state(composition)


def test_an_unknown_method():
    with pytest.raises(OptmatError, match="unknown density method"):
        mixture_equation_of_state({"Ar": 1.0}, method="ideal")


def test_the_phase_is_checked():
    co2 = _quiet({"CO2": 0.9, "Ar": 0.1}, method="multi-fluid")
    with pytest.raises(StateError, match="no equation-of-state solution"):
        co2.number_density(State(pressure=8e6, temperature=260.0))
    fallback = _quiet({"CO2": 0.9, "Ar": 0.1}, method="partial-density")
    with pytest.raises(StateError, match="beyond its gas phase"):
        fallback.pressure(1.2e28, 260.0)
    with pytest.raises(StateError, match="outside its"):
        co2.number_density(State(pressure=1e5, temperature=5000.0))
    with pytest.raises(StateError, match="density > 0"):
        co2.pressure(0.0, 300.0)


def test_the_truncation_bound_warns():
    # A cross term set to make rho sum x x B exceed the bound at 1e27 m^-3.
    eos = _quiet({"He": 0.5, "Ar": 0.5}, method="partial-density")
    big = TRUNCATION_BOUND / (1e27 / N_A) / 0.5 * 4
    object.__setattr__(eos, "_cross", lambda temperature: {("Ar", "He"): big})
    with pytest.warns(ModelRangeWarning, match="truncation bound"):
        eos.pressure(1e27, 294.0)


def test_the_pair_file_is_checked(monkeypatch, tmp_path):
    path = tmp_path / "pairs.toml"
    path.write_text('[pairs.He-Ar]\nmodel = "x"\nsource = "x"\nqualified = false\n')
    monkeypatch.setattr(mixture_eos, "PAIRS_FILE", path)
    mixture_eos._pair_records.cache_clear()
    try:
        with pytest.raises(DataError, match="sorted order"):
            mixture_eos._pair_records()
        path.write_text(
            '[pairs.Ar-He]\nmodel = "x"\nsource = "x"\nqualified = false\n'
            'reason = "r"\ncross_virial = "achtermann_1991/A_R_this_work_Ar"\n'
        )
        mixture_eos._pair_records.cache_clear()
        with pytest.raises(DataError, match="not a cross virial of it"):
            mixture_eos._pair_records()
        path.write_text(
            '[pairs.Ar-He]\nmodel = "x"\nsource = "x"\nqualified = false\n'
            'reason = "r"\n'
        )
        mixture_eos._pair_records.cache_clear()
        with pytest.raises(DataError, match="but CoolProp's is 'Tkaczuk-JPCRD-2020'"):
            mixture_equation_of_state({"He": 0.5, "Ar": 0.5})
    finally:
        mixture_eos._pair_records.cache_clear()


def test_the_ideal_limit():
    eos = _quiet({"He": 0.5, "Ar": 0.5}, method="multi-fluid")
    state = State(pressure=10.0, temperature=294.0)
    ideal = state.pressure / (K_B * state.temperature)
    # At 10 Pa, B rho is about 1e-8.
    assert eos.number_density(state) == pytest.approx(ideal, rel=1e-7)


@pytest.mark.parametrize(
    ("composition", "printed", "tolerance"),
    [
        # Tkaczuk et al. 2020, Table 8: p at 200 K and 10 mol/dm^3, equimolar. The He-Ne
        # value agrees to 3.4e-7; both argon pairs lie 2.7e-6 below, the same offset,
        # whose cause is not established (TC-779).
        ({"He": 0.5, "Ne": 0.5}, 18430775.292601, 5e-7),
        ({"He": 0.5, "Ar": 0.5}, 17128034.388363, 3e-6),
        ({"Ne": 0.5, "Ar": 0.5}, 15905875.375781, 3e-6),
    ],
)
def test_coolprop_reproduces_tkaczuks_check_values(composition, printed, tolerance):
    eos = _quiet(composition, method="multi-fluid")
    backend = eos._backend()
    backend.update(mixture_eos._coolprop().DmolarT_INPUTS, 10000.0, 200.0)
    assert float(backend.p()) == pytest.approx(printed, rel=tolerance)


def test_the_fallback_brackets_far_from_the_ideal_gas():
    # Supercritical CO2 at 310 K and 100 bar has Z near 0.3, helium at 3000 bar near
    # 2: the solve widens its bracket both ways and still round-trips.
    for composition, state in (
        ({"CO2": 0.99, "Ar": 0.01}, State(pressure=1e7, temperature=310.0)),
        ({"He": 0.99, "Ne": 0.01}, State(pressure=3e8, temperature=294.0)),
    ):
        eos = _quiet(composition, method="partial-density")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ModelRangeWarning)
            density = eos.number_density(state)
            back = float(eos.pressure(density, state.temperature))
        assert back == pytest.approx(state.pressure, rel=1e-10)
        z = state.pressure / (density * K_B * state.temperature)
        assert not 0.5 < z < 2


def test_a_failed_flash_is_a_state_error(monkeypatch):
    eos = _quiet({"He": 0.5, "Ar": 0.5}, method="multi-fluid")

    class Refusing:
        def update(self, *args):
            raise ValueError("refused")

    monkeypatch.setattr(eos.__class__, "_backend", lambda self: Refusing())
    with pytest.raises(StateError, match="refused"):
        eos.number_density(AMBIENT)


def test_a_measured_cross_coefficient_as_the_fallbacks(monkeypatch):
    # Vatter et al.'s B_12 of Ar-O2 (296-463 K) in place of Martin et al.'s E.
    records = dict(mixture_eos._pair_records())
    records["Ar", "O2"] = dict(records["Ar", "O2"]) | {
        "cross_virial": "vatter_1996/table1_B12_Ar_O2"
    }
    monkeypatch.setattr(mixture_eos, "_pair_records", lambda: records)
    eos = _quiet({"Ar": 0.5, "O2": 0.5}, method="virial")
    assert eos.temperature_range == (296.15, 463.15)
    assert eos._cross_sum(296.15) == pytest.approx(2 * 0.25 * -16.4e-6, rel=1e-12)

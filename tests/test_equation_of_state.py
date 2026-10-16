# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The real-gas equation of state, its inverse and its phase domain (OEP 1, §4.2, §6.5).

The inverse must round-trip on every gas's phase domain, keep the CODATA convention and
recover the ideal gas at vanishing pressure; condensed states and states beyond the
equation's pressure are refused forward and inverse alike.
"""

import dataclasses
import math

import numpy as np
import pytest

from optmat.conditions import ONE_ATM_0C, State
from optmat.constants import K_B, N_A
from optmat.exceptions import StateError
from optmat.gas import eos
from optmat.gas.eos import (
    ROUND_TRIP_TOLERANCE,
    air_equation_of_state,
    equation_of_state,
)
from optmat.gas.gases import Gas, info
from optmat.gas.transport import TRANSPORT_GASES, transport_model

GASES = [g for g in Gas if info(g).fluid is not None]


def _temperatures(e):
    low, high = e.temperature_range
    candidates = [low * 1.2, 0.9 * e.critical_temperature, 294.0, 400.0]
    return [t for t in candidates if low < t < high]


@pytest.mark.parametrize("gas", GASES, ids=lambda g: g.value)
def test_the_inverse_round_trips_on_the_phase_domain(gas):
    e = equation_of_state(gas)
    for temperature in _temperatures(e):
        _, top = e.phase_domain(temperature)
        # From 0.01 bar to the domain's top below the vapour pressure, or to 100 bar.
        high = top * (1 - 1e-9) if temperature < e.critical_temperature else 1e7
        for p in np.geomspace(1e3, min(high, 1e7), 7):
            n = e.number_density(State(float(p), temperature))
            back = float(e.pressure(n, temperature))
            again = e.number_density(State(back, temperature))
            # The solve's stated bound; the pressure to the same accuracy over the
            # compressibility's slope, which is at most a few near saturation.
            assert abs(again / n - 1) <= ROUND_TRIP_TOLERANCE
            assert back == pytest.approx(float(p), rel=1e-9)


EQUATIONS = [equation_of_state(g) for g in GASES] + [eos.air_equation_of_state()]
BOUNDS = r"(vapour|melting) pressure|equation's pressure"


def _melting_line(e):
    """Return CoolProp's melting line of an equation's fluid, or None."""
    backend = eos._abstract_state(e.fluid)
    return backend if backend.has_melting_line() else None


def _scan_temperatures(e):
    low, high = e.temperature_range
    top = min(400.0, high)
    temperatures = [low * (1 + 1e-9), e.critical_temperature]
    temperatures += [e.critical_temperature * (1 + 1e-9), top]
    temperatures += list(np.arange(math.ceil(low), top, 4.0))
    line = _melting_line(e)
    if line is not None:
        # Either side of where the melting line crosses the equation's pressure.
        cp = eos._coolprop()
        crossing = line.melting_line(cp.iT, cp.iP, e.maximum_pressure)
        temperatures += [crossing - 0.1, crossing + 0.1]
    return sorted(t for t in temperatures if low < t <= top)


@pytest.mark.parametrize("e", EQUATIONS, ids=lambda e: e.material)
def test_the_inverse_round_trips_across_the_domain(e):
    # Every temperature from just above the triple point to 400 K, and every density
    # from 1e20 m^-3 to just below the domain's top: the solve's stated bound holds,
    # and the density at the top, and above it, is refused by name (TC-811).
    for temperature in _scan_temperatures(e):
        top = e.phase_domain(temperature)[1]
        excluded = temperature < e.critical_temperature and top < (
            e._melting_pressure(temperature) or math.inf
        )
        dense = e._density_at(top * (1 - 1e-12) if excluded else top, temperature)
        for n in np.geomspace(1e20, dense * (1 - 1e-9), 6):
            back = float(e.pressure(n, temperature))
            again = e.number_density(State(back, temperature))
            assert abs(again / n - 1) <= ROUND_TRIP_TOLERANCE
        with pytest.raises(StateError, match=BOUNDS):
            e.pressure(dense * 1.001, temperature)


@pytest.mark.parametrize("e", EQUATIONS, ids=lambda e: e.material)
def test_the_vapour_lies_below_the_melting_line(e):
    # Below the critical temperature the vapour pressure binds, not the melting line:
    # asserted, not assumed. CO2's equation puts its triple point at 517 964 Pa and its
    # melting line at 517 950 Pa, 2.7e-5 apart, so the comparison allows 3e-5.
    for temperature in _scan_temperatures(e):
        saturation = e.saturation_pressure(temperature)
        melting = e._melting_pressure(temperature)
        if saturation is not None and melting is not None:
            assert saturation < melting * (1 + 3e-5)


@pytest.mark.parametrize(
    "e", [e for e in EQUATIONS if _melting_line(e)], ids=lambda e: e.material
)
def test_the_melting_lines_reach_the_equations_pressure(e):
    # Above a line's own temperatures its pressure exceeds the line's top, so the
    # equation's pressure binds there only if the line reaches it.
    cp = eos._coolprop()
    line = _melting_line(e)
    assert line is not None
    assert line.melting_line(cp.iP_max, -1, 0.0) >= e.maximum_pressure


class TestTheMeltingLine:
    @pytest.mark.parametrize(
        ("gas", "temperature"), [("Ar", 250.0), ("N2", 280.0), ("CH4", 255.0)]
    )
    def test_the_reported_states_round_trip(self, gas, temperature):
        # solprop SEP 3, Q2: every density raised at these temperatures.
        e = equation_of_state(gas)
        n = e.number_density(State(1e5, temperature))
        back = float(e.pressure(n, temperature))
        assert abs(e.number_density(State(back, temperature)) / n - 1) <= (
            ROUND_TRIP_TOLERANCE
        )
        assert back == pytest.approx(1e5, rel=1e-9)

    def test_the_melting_pressure_bounds_the_supercritical_fluid(self):
        nitrogen = equation_of_state("N2")
        melting = nitrogen._melting_pressure(280.0)
        assert melting is not None
        top = nitrogen.phase_domain(280.0)[1]
        assert top == pytest.approx(melting * (1 - 1e-5), rel=1e-15)
        assert top < nitrogen.maximum_pressure
        solid = State((melting + nitrogen.maximum_pressure) / 2, 280.0)
        with pytest.raises(StateError, match="solid"):
            nitrogen.number_density(solid)
        dense = nitrogen.number_density(State(top, 280.0))
        assert float(nitrogen.pressure(dense, 280.0)) == pytest.approx(top, rel=1e-9)
        with pytest.raises(StateError, match="its melting pressure"):
            nitrogen.pressure(dense * 1.001, 280.0)

    def test_without_a_line_the_equations_pressure_binds(self):
        xenon = equation_of_state("Xe")
        assert xenon._melting_pressure(300.0) is None
        assert xenon.phase_domain(300.0)[1] == xenon.maximum_pressure
        # Krypton's line ends at 170 K, at 205 MPa, above its equation's 200 MPa.
        krypton = equation_of_state("Kr")
        assert krypton._melting_pressure(250.0) is None
        assert krypton.phase_domain(250.0)[1] == krypton.maximum_pressure

    def test_the_triple_point_is_refused(self):
        argon = equation_of_state("Ar")
        triple = argon.temperature_range[0]
        with pytest.raises(StateError, match="triple point"):
            argon.phase_domain(triple)
        with pytest.raises(StateError, match="triple point"):
            argon.number_density(State(1e3, triple))
        with pytest.raises(StateError, match="triple point"):
            argon.pressure(1e23, triple)


@pytest.mark.parametrize("gas", GASES, ids=lambda g: g.value)
def test_the_ideal_gas_is_the_limit_at_vanishing_pressure(gas):
    e = equation_of_state(gas)
    temperature = _temperatures(e)[-1]
    state = State(1.0, temperature)
    # Z - 1 is proportional to the pressure: B p / (R T) with |B| below 1e-3 m^3/mol for
    # every gas here, so at 1 Pa below 1e-6.
    assert e.compressibility(state) == pytest.approx(1.0, abs=1e-6)


def test_the_inverse_keeps_the_codata_convention():
    argon = equation_of_state("Ar")
    state = State(1e5, 294.0)
    n = argon.number_density(state)
    assert n == pytest.approx(
        state.pressure / (argon.compressibility(state) * K_B * 294.0)
    )
    # The equation's own gas constant (Tegeler et al.: 8.31451 J/(mol K)) puts its molar
    # density 5.7e-6 from the CODATA one; the inverse follows N = p / (Z k_B T).
    assert abs(argon.gas_constant / (N_A * K_B) - 1) == pytest.approx(5.7e-6, rel=0.05)
    assert float(argon.pressure(n, 294.0)) == pytest.approx(1e5, rel=1e-12)


class TestTheDomain:
    def test_a_liquid_is_refused_forward(self):
        # SF6 at 294 K condenses near 21.8 bar, CO2 at 280 K near 41.6 bar (TC-731).
        with pytest.raises(StateError, match="condensed"):
            equation_of_state("SF6").number_density(State(25e5, 294.0))
        with pytest.raises(StateError, match="condensed"):
            eos.number_density("CO2", State(50e5, 280.0))

    def test_a_supercritical_liquid_is_refused(self):
        # Below T_c (304.1 K) and above p_c (73.8 bar): a compressed liquid.
        with pytest.raises(StateError, match="condensed"):
            eos.number_density("CO2", State(100e5, 290.0))

    def test_the_saturated_vapour_and_denser_states_are_refused_inverse(self):
        sf6 = equation_of_state("SF6")
        saturation = sf6.saturation_pressure(294.0)
        assert saturation is not None
        _, top = sf6.phase_domain(294.0)
        assert top == pytest.approx(saturation * (1 - eos.SATURATION_MARGIN))
        vapour = sf6.number_density(State(top * (1 - 1e-9), 294.0))
        with pytest.raises(StateError, match="top of its phase domain"):
            sf6.pressure(vapour * 1.001, 294.0)
        assert float(sf6.pressure(vapour, 294.0)) < top
        with pytest.raises(StateError, match="condensed"):
            sf6.number_density(State(saturation * (1 - 1e-7), 294.0))

    def test_the_equations_pressure_bounds_the_supercritical_fluid(self):
        oxygen = equation_of_state("O2")  # p_max 80 MPa
        with pytest.raises(StateError, match="above its equation's pressure"):
            oxygen.number_density(State(oxygen.maximum_pressure * 1.1, 300.0))
        dense = oxygen.number_density(State(oxygen.maximum_pressure, 300.0))
        with pytest.raises(StateError, match="top of its phase domain"):
            oxygen.pressure(dense * 1.01, 300.0)

    def test_temperatures_outside_the_equation_are_refused(self):
        with pytest.raises(StateError, match="outside its"):
            eos.number_density("He", State(1e5, 1.0))
        with pytest.raises(StateError, match="outside its"):
            equation_of_state("He").pressure(1e25, 1.0)

    def test_the_vacuum_is_refused_by_the_inverse(self):
        with pytest.raises(StateError, match="density > 0"):
            eos.pressure("Ar", 0.0, 294.0)

    def test_the_domain_below_and_above_the_critical_point(self):
        co2 = equation_of_state("CO2")
        assert co2.phase_domain(280.0)[1] == pytest.approx(41.6e5, rel=2e-3)
        # Above 327.7 K CO2's melting line is above its equation's 800 MPa; between
        # T_c and there, the line binds (TC-811).
        assert co2.phase_domain(340.0) == (0.0, co2.maximum_pressure)
        assert co2.phase_domain(320.0)[1] < co2.maximum_pressure
        assert co2.saturation_pressure(320.0) is None


def test_the_inverse_takes_arrays():
    argon = equation_of_state("Ar")
    n = np.array([[1e24, 1e25], [2e25, 3e25]])
    p = argon.pressure(n, 294.0)
    assert p.shape == (2, 2)
    assert np.all(np.diff(p.ravel()) > 0)
    ideal = eos.pressure("Ar", n, 294.0, ideal=True)
    np.testing.assert_allclose(ideal, n * K_B * 294.0)


def test_the_solve_widens_a_poor_start(monkeypatch):
    argon = equation_of_state("Ar")
    n = argon.number_density(State(5e5, 294.0))
    original = argon._start
    monkeypatch.setattr(
        type(argon), "_start", lambda self, d, t: 3 * original(d, t), raising=True
    )
    assert float(argon.pressure(n, 294.0)) == pytest.approx(5e5, rel=1e-12)


def test_a_solve_that_misses_the_round_trip_raises(monkeypatch):
    monkeypatch.setattr(eos, "ROUND_TRIP_TOLERANCE", -1.0)
    with pytest.raises(StateError, match="missed the round trip"):
        eos.pressure("Ar", 2.5e25, 294.0)


def test_the_speed_of_sound_tends_to_the_ideal_gas():
    argon = equation_of_state("Ar")
    temperature = 294.0
    ideal = math.sqrt(5 / 3 * N_A * K_B * temperature / argon.molar_mass)
    # At 1 kPa the real gas's second virial moves it by below 1e-5.
    assert argon.speed_of_sound(State(1e3, temperature)) == pytest.approx(
        ideal, rel=1e-4
    )


class TestTheMolarMass:
    def test_it_is_the_equations_own(self):
        # Tegeler et al.'s 39.948 g/mol, as CoolProp holds it.
        assert equation_of_state("Ar").molar_mass == 0.039948

    @pytest.mark.parametrize("gas", TRANSPORT_GASES)
    def test_it_is_the_transport_models(self, gas):
        assert equation_of_state(gas).molar_mass == transport_model(gas).molar_mass

    def test_air_is_the_pseudo_pure_fluids_not_lemmons_printed_one(self):
        # CoolProp's pseudo-pure air holds 28.96546 g/mol; Lemmon et al. (2000) print
        # 28.9586 g/mol (fixed points, p. 334), their composition's mean (TC-832).
        air = air_equation_of_state().molar_mass
        assert air == 0.02896546
        assert air / 0.0289586 - 1 == pytest.approx(2.37e-4, rel=1e-2)

    def test_the_speed_of_sound_uses_it(self):
        # w^2 = (c_p / c_v) (dp/drho)_T with rho the mass density, rho_molar * M: the
        # backend's speed of sound is consistent with this molar mass, to rounding.
        air = air_equation_of_state()
        backend = eos._abstract_state(air.fluid)
        coolprop = eos._coolprop()
        backend.update(coolprop.PT_INPUTS, 101325.0, 294.0)
        slope = backend.first_partial_deriv(coolprop.iP, coolprop.iDmolar, coolprop.iT)
        ratio = backend.cpmolar() / backend.cvmolar()
        assert air.speed_of_sound(State(101325.0, 294.0)) == pytest.approx(
            math.sqrt(ratio * slope / air.molar_mass), rel=1e-14
        )

    def test_it_is_not_in_the_identity(self):
        # Resolved from the backend, which the fluid and the backend's version fix.
        assert "molar_mass" not in air_equation_of_state().identity_record()["fields"]


def test_the_identity_names_the_backend_and_the_equation():
    argon = equation_of_state("Ar")
    identity = argon.identity
    assert identity.kind == "equation_of_state"
    assert identity.status == "validated"
    assert argon.reference == "Tegeler-JPCRD-1999"
    assert argon.backend.startswith("CoolProp ")
    other = dataclasses.replace(argon, backend="CoolProp 0.0")
    assert other.identity.digest != identity.digest


def test_a_state_the_backend_cannot_solve_raises():
    with pytest.raises(StateError, match="no equation-of-state solution"):
        eos._compressibility("Argon", -1.0, 294.0)


def test_ozone_has_no_equation_of_state():
    with pytest.raises(StateError, match="no equation of state for O3"):
        equation_of_state("O3")


def test_the_forward_law_is_the_one_every_reduction_used():
    # The reductions of the curated data read the same Z (the golden values of every
    # model are unchanged); spot-check the familiar value.
    assert eos.compressibility("He", ONE_ATM_0C) == pytest.approx(1.00053, abs=1e-5)


def test_without_coolprop_the_error_says_what_to_install(monkeypatch):
    import sys  # noqa: PLC0415

    eos._coolprop.cache_clear()
    monkeypatch.setitem(sys.modules, "CoolProp", None)
    monkeypatch.setitem(sys.modules, "CoolProp.CoolProp", None)
    try:
        with pytest.raises(ImportError, match=r"install optmat\[gas\]"):
            eos._coolprop()
    finally:
        monkeypatch.undo()
        eos._coolprop.cache_clear()
    assert eos.number_density("Ar", State(1e5, 300.0), ideal=True) > 0


def test_the_cached_states_are_released(monkeypatch):
    eos.number_density("Ar", State(1e5, 300.0))
    eos._release()
    assert eos._abstract_state.cache_info().currsize == 0

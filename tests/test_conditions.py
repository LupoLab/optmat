# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from optmat.conditions import (
    ONE_ATM_0C,
    ONE_ATM_15C,
    ONE_ATM_20C,
    ONE_BAR_0C,
    State,
    standard_air_refractivity,
    vacuum_wavelength_from_standard_air,
)
from optmat.constants import K_B
from optmat.exceptions import StateError
from optmat.gas.eos import compressibility, ideal_number_density, number_density
from optmat.gas.gases import Gas


def test_named_states():
    assert State(101325.0, 273.15) == ONE_ATM_0C
    assert State(1e5, 273.15) == ONE_BAR_0C
    assert ONE_ATM_15C.temperature == pytest.approx(288.15, abs=1e-12)
    assert ONE_ATM_20C.temperature == pytest.approx(293.15, abs=1e-12)


@pytest.mark.parametrize(
    ("pressure", "temperature"),
    [(0.0, 273.15), (-1.0, 273.15), (1e5, 0.0), (np.nan, 273.15), (1e5, np.inf)],
)
def test_state_rejects_unphysical_values(pressure, temperature):
    with pytest.raises(StateError, match="finite and positive"):
        State(pressure, temperature)


def test_ideal_density_is_p_over_kt():
    state = State(5e4, 300.0)
    assert ideal_number_density(state) == pytest.approx(
        5e4 / (K_B * 300.0), rel=1e-15, abs=0.0
    )
    assert number_density("Ar", state, ideal=True) == ideal_number_density(state)


@pytest.mark.parametrize("gas", [g for g in Gas if g is not Gas.O3])
def test_real_density_approaches_ideal_at_low_pressure(gas):
    # The second virial coefficient sets 1 - Z ~ B p / (R T): at 10 Pa it is below
    # 2e-6 for every gas here (|B| < 4e-4 m^3/mol even for SF6 at 273 K).
    state = State(10.0, 273.15)
    ratio = number_density(gas, state) / ideal_number_density(state)
    assert ratio == pytest.approx(1.0, abs=2e-6)


def test_compressibility_signs_at_stp():
    # Helium is repulsive (B > 0, Z > 1) at 0 degrees C; xenon strongly attractive.
    assert compressibility("He", ONE_ATM_0C) > 1.0
    assert compressibility("Xe", ONE_ATM_0C) < 0.995


def test_ozone_has_no_real_gas_density():
    with pytest.raises(StateError, match="no equation of state for O3"):
        number_density("O3", ONE_ATM_0C)
    assert number_density("O3", ONE_ATM_0C, ideal=True) > 0


def test_state_outside_the_equation_of_state_domain_raises():
    with pytest.raises(StateError, match="no equation-of-state solution"):
        number_density("He", State(1e14, 1.0))


def test_standard_air_refractivity_is_normal_and_of_the_right_size():
    # The coefficients themselves are queued for checking against Edlen 1966 (TC-2);
    # here only the physics: n - 1 of air near 2.8e-4, falling with wavelength.
    wavelength = np.linspace(250e-9, 2e-6, 50)
    dn = standard_air_refractivity(1 / wavelength)
    assert np.all(np.diff(dn) < 0)
    assert 2.7e-4 < dn[-1] < dn[0] < 3.1e-4


@given(st.floats(min_value=200e-9, max_value=2e-6))
def test_air_to_vacuum_conversion_is_self_consistent(air):
    vacuum = float(vacuum_wavelength_from_standard_air(air))
    residual = vacuum - air * (1 + float(standard_air_refractivity(1 / vacuum)))
    # Four fixed-point iterations contract by ~dn/dlnlambda ~ 1e-5 each: the residual
    # is at double-precision round-off of a ~1e-6 m number.
    assert abs(residual) < 1e-20
    assert vacuum > air


def test_air_to_vacuum_keeps_shape():
    air = np.array([[400e-9, 500e-9], [600e-9, 700e-9]])
    assert vacuum_wavelength_from_standard_air(air).shape == (2, 2)


def test_air_refractivity_is_standard_air_at_15c_and_760_torr():
    from optmat.conditions import ONE_ATM_15C, air_refractivity  # noqa: PLC0415

    wavenumber = 1 / np.array([300e-9, 633e-9, 1500e-9])
    # 720.775 is Edlen's rounded density factor of standard air; his Eq. (12) returns
    # the standard refractivity to its fifth digit (4e-7 relative, 1e-10 in n - 1).
    ratio = air_refractivity(wavenumber, ONE_ATM_15C) / standard_air_refractivity(
        wavenumber
    )
    assert ratio == pytest.approx(1.0, abs=1e-6)


def test_air_refractivity_scales_with_density():
    from optmat.conditions import State, air_refractivity  # noqa: PLC0415

    wavenumber = 1 / 633e-9
    cold = air_refractivity(wavenumber, State(101325.0, 278.15))
    warm = air_refractivity(wavenumber, State(101325.0, 303.15))
    # The ideal-gas ratio 303.15/278.15 = 1.0899; Edlen's non-ideal term,
    # p (0.817 - 0.0133 t) 1e-6, adds 760 torr x 0.0133 x 25 K x 1e-6 = 2.5e-4.
    expected = 303.15 / 278.15 * (1 + 760 * 0.0133 * 25 * 1e-6)
    assert float(cold / warm) == pytest.approx(expected, rel=1e-5, abs=0.0)
    half = air_refractivity(wavenumber, State(50662.5, 293.15))
    full = air_refractivity(wavenumber, State(101325.0, 293.15))
    assert float(full / half) == pytest.approx(2.0, rel=1e-3, abs=0.0)

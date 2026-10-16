# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import pytest

from optmat import constants


def test_atomic_unit_of_polarizability():
    # 4 pi eps0 a0^3: CODATA 2022's atomic unit of electric polarizability,
    # 1.64877727212(51)e-41 C^2 m^2 J^-1; tolerance its 12 printed digits.
    assert (
        pytest.approx(1.64877727212e-41, rel=1e-10, abs=0.0)
        == constants.AU_POLARIZABILITY
    )


def test_atomic_unit_of_angular_frequency():
    assert (
        pytest.approx(4.134137333518e16, rel=1e-11, abs=0.0)
        == constants.AU_ANGULAR_FREQUENCY
    )


def test_pressure_units():
    assert constants.STD_ATMOSPHERE == 101325.0
    assert constants.BAR == 1e5
    assert pytest.approx(133.322368, rel=1e-8, abs=0.0) == constants.TORR
    assert math.isclose(constants.ZERO_CELSIUS, 273.15)

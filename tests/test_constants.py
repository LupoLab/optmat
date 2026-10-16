# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import pytest

from optmat import constants


def test_atomic_unit_of_polarizability():
    # 4 pi eps0 a0^3 = 1.64877727e-41 C m^2/V (CODATA 2018); tolerance the 9 digits.
    assert pytest.approx(1.64877727436e-41, rel=1e-9) == constants.AU_POLARIZABILITY


def test_atomic_unit_of_angular_frequency():
    assert pytest.approx(4.134137333518e16, rel=1e-11) == constants.AU_ANGULAR_FREQUENCY


def test_pressure_units():
    assert constants.STD_ATMOSPHERE == 101325.0
    assert constants.BAR == 1e5
    assert pytest.approx(133.322368, rel=1e-8) == constants.TORR
    assert math.isclose(constants.ZERO_CELSIUS, 273.15)

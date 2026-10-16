# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Physical constants and unit conversions, SI throughout.

Thin, named wrappers over :mod:`scipy.constants` (CODATA) so that no module spells a
constant by hand, plus the atomic-unit conversions the polarizability data need. The
historical pressure units (atmosphere, torr) and the named reference temperatures are
here too, because the literature reduces its measurements to states written in them.
"""

import math

from scipy import constants as _constants

__all__ = [
    "ATOMIC_MASS",
    "AU_ANGULAR_FREQUENCY",
    "AU_ENERGY",
    "AU_FIELD",
    "AU_POLARIZABILITY",
    "BAR",
    "BOHR_RADIUS",
    "CLASSICAL_ELECTRON_RADIUS",
    "EPS0",
    "E_CHARGE",
    "HBAR",
    "K_B",
    "M_E",
    "N_A",
    "STD_ATMOSPHERE",
    "TORR",
    "ZERO_CELSIUS",
    "C",
]

C: float = _constants.c
"""Speed of light in vacuum, m/s."""

EPS0: float = _constants.epsilon_0
"""Vacuum permittivity, F/m."""

HBAR: float = _constants.hbar
"""Reduced Planck constant, J s."""

E_CHARGE: float = _constants.e
"""Elementary charge, C."""

M_E: float = _constants.m_e
"""Electron mass, kg."""

K_B: float = _constants.k
"""Boltzmann constant, J/K."""

N_A: float = _constants.N_A
"""Avogadro constant, 1/mol."""

ATOMIC_MASS: float = _constants.physical_constants["atomic mass constant"][0]
"""Atomic mass constant (unified atomic mass unit), kg."""

BOHR_RADIUS: float = _constants.physical_constants["Bohr radius"][0]
"""Bohr radius a0, m."""

CLASSICAL_ELECTRON_RADIUS: float = _constants.physical_constants[
    "classical electron radius"
][0]
"""Classical electron radius r_e = e^2 / (4 pi eps0 m_e c^2), m."""

AU_ENERGY: float = _constants.physical_constants["Hartree energy"][0]
"""Hartree energy, J: the atomic unit of energy."""

AU_ANGULAR_FREQUENCY: float = AU_ENERGY / HBAR
"""Atomic unit of angular frequency, rad/s (one hartree divided by hbar)."""

AU_POLARIZABILITY: float = 4 * math.pi * EPS0 * BOHR_RADIUS**3
"""Atomic unit of dipole polarizability, C m^2/V.

An atomic-unit polarizability ``alpha_au`` (in units of a0^3 in the Gaussian sense) is
``alpha_au * AU_POLARIZABILITY`` in SI.
"""

AU_FIELD: float = AU_ENERGY / (E_CHARGE * BOHR_RADIUS)
"""Atomic unit of electric field, V/m (E_h / (e a0), 5.14e11 V/m)."""

STD_ATMOSPHERE: float = _constants.atm
"""Standard atmosphere, Pa (101 325 exactly)."""

BAR: float = _constants.bar
"""Bar, Pa (100 000 exactly)."""

TORR: float = STD_ATMOSPHERE / 760
"""Torr (mm of mercury in the older literature, to the precision that matters), Pa."""

ZERO_CELSIUS: float = _constants.zero_Celsius
"""The ice point, 273.15 K."""

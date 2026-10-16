# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""States (pressure and temperature) and air-to-vacuum wavelengths.

Measurements are reported at, or reduced to, a stated pressure and temperature ("0
degrees C and 760 Torr", "15 degrees C and 1 atm", "1 bar and 0 degrees C"); this module
holds the state and the named states. What a state implies for a medium (the number
density of a gas, :mod:`optmat.gas.eos`) belongs to that medium's subpackage.

Pressures here are in pascals, like everything else in the package.

Much of the older literature (Leonard 1974, Cuthbertson) gives wavelengths in *standard
air* (15 degrees C, 101 325 Pa); :func:`vacuum_wavelength_from_standard_air` converts
them with Edlen's 1966 dispersion formula for standard air. Solids measured by minimum
deviation in room air give an index relative to that air; :func:`air_refractivity`
gives its refractivity at the room's temperature and pressure.
"""

from dataclasses import dataclass
from typing import Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat.constants import BAR, STD_ATMOSPHERE, TORR, ZERO_CELSIUS
from optmat.exceptions import StateError

__all__ = [
    "EVALUATION_VERSION",
    "ONE_ATM_0C",
    "ONE_ATM_15C",
    "ONE_ATM_20C",
    "ONE_BAR_0C",
    "State",
    "air_refractivity",
    "standard_air_refractivity",
    "vacuum_wavelength_from_standard_air",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""


@dataclass(frozen=True)
class State:
    """A thermodynamic state: a pressure and a temperature.

    Attributes
    ----------
    pressure : float
        Absolute pressure in Pa, > 0.
    temperature : float
        Temperature in K, > 0.

    Raises
    ------
    StateError
        If either value is not finite and positive.

    Examples
    --------
    >>> State(101325.0, 273.15).pressure
    101325.0
    """

    pressure: float
    temperature: float

    def __post_init__(self) -> None:
        """Reject a non-physical state."""
        for name in ("pressure", "temperature"):
            value = getattr(self, name)
            if not (np.isfinite(value) and value > 0):
                msg = f"{name} must be finite and positive, got {value!r}"
                raise StateError(msg)


ONE_ATM_0C: Final = State(STD_ATMOSPHERE, ZERO_CELSIUS)
"""0 degrees C and 1 atm (760 Torr): the "STP" or "NTP" of most refractivity papers."""

ONE_BAR_0C: Final = State(BAR, ZERO_CELSIUS)
"""0 degrees C and 1 bar: Borzsonyi et al. (2008) and Luna's normalisation."""

ONE_ATM_15C: Final = State(STD_ATMOSPHERE, ZERO_CELSIUS + 15.0)
"""15 degrees C and 1 atm: the atmospheric-science convention (He et al. 2021)."""

ONE_ATM_20C: Final = State(STD_ATMOSPHERE, ZERO_CELSIUS + 20.0)
"""20 degrees C and 1 atm: the frequency-comb measurements (Zhang et al. 2008)."""


def standard_air_refractivity(
    vacuum_wavenumber: ArrayLike | float,
) -> Float[np.ndarray, "*shape"]:
    """Return ``n - 1`` of standard air from Edlen's 1966 formula.

    Standard air is dry air with 0.03 % CO2 at 15 degrees C and 101 325 Pa. The formula
    is ``1e8 (n - 1) = 8342.13 + 2406030 / (130 - s^2) + 15997 / (38.9 - s^2)`` with
    ``s`` the vacuum wavenumber in um^-1 (B. Edlen, Metrologia 2, 71 (1966)).

    Parameters
    ----------
    vacuum_wavenumber : array_like or float
        Vacuum wavenumber ``1 / lambda_vac`` in m^-1; the formula is fitted over
        roughly 200 nm to 2 um.

    Returns
    -------
    numpy.ndarray
        ``n - 1`` of standard air, the shape of the input.

    Examples
    --------
    >>> f"{float(standard_air_refractivity(1 / 546.2e-9)):.4e}"
    '2.7790e-04'
    """
    s2 = (np.asarray(vacuum_wavenumber, dtype=float) * 1e-6) ** 2
    return np.asarray(
        1e-8 * (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2))
    )


def air_refractivity(
    vacuum_wavenumber: ArrayLike | float, state: State
) -> Float[np.ndarray, "*shape"]:
    """Return ``n - 1`` of dry air at a temperature and pressure, from Edlen 1966.

    Edlen's Eq. (12) (Metrologia 2, 71 (1966), p. 78) scales the standard-air
    refractivity by the density factor of dry air,
    ``(n - 1)_tp = p (n - 1)_s / 720.775 * [1 + p (0.817 - 0.0133 t) 1e-6] /
    (1 + 0.0036610 t)``, with ``t`` in degrees C and ``p`` in torr, for 5 to 30
    degrees C. It returns :func:`standard_air_refractivity` at 15 degrees C and 760
    torr. Water vapour and a CO2 content other than 0.03 % are not included (Edlen's
    corrections for them are some 1e-7 at room conditions).

    Parameters
    ----------
    vacuum_wavenumber : array_like or float
        Vacuum wavenumber ``1 / lambda_vac`` in m^-1 (roughly 200 nm to 2 um, the range
        of the standard-air formula).
    state : State
        The air's pressure (Pa) and temperature (K); the formula is stated for 5 to 30
        degrees C and laboratory pressures.

    Returns
    -------
    numpy.ndarray
        ``n - 1`` of the air, the shape of the input.

    Examples
    --------
    >>> f"{float(air_refractivity(1 / 800e-9, ONE_ATM_20C)):.4e}"
    '2.7033e-04'
    """
    t = state.temperature - ZERO_CELSIUS
    p = state.pressure / TORR
    factor = p / 720.775 * (1 + p * (0.817 - 0.0133 * t) * 1e-6) / (1 + 0.0036610 * t)
    return np.asarray(standard_air_refractivity(vacuum_wavenumber) * factor)


def vacuum_wavelength_from_standard_air(
    air_wavelength: ArrayLike | float,
) -> Float[np.ndarray, "*shape"]:
    """Convert wavelengths measured in standard air to vacuum wavelengths.

    Solves ``lambda_vac = n(lambda_vac) lambda_air`` with
    :func:`standard_air_refractivity` by fixed-point iteration; the refractivity changes
    by ~1e-10 per iteration after the first, so four iterations reach double precision.

    Parameters
    ----------
    air_wavelength : array_like or float
        Wavelength in standard air, m.

    Returns
    -------
    numpy.ndarray
        Vacuum wavelength, m, the shape of the input.

    Examples
    --------
    >>> f"{float(vacuum_wavelength_from_standard_air(546.0740e-9)) * 1e9:.4f}"
    '546.2258'
    """
    air = np.asarray(air_wavelength, dtype=float)
    vacuum = air.copy()
    for _ in range(4):
        vacuum = air * (1.0 + standard_air_refractivity(1.0 / vacuum))
    return np.asarray(vacuum)

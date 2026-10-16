# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The widths of Raman lines against density, temperature and collision partner.

A Raman line's full width at half maximum (FWHM) in a gas of number density :math:`n`
(in amagat) is written

.. math::

   \Delta\nu = c + \frac{A(T)}{n} + n \sum_p x_p\, \gamma_p(J, T),

with :math:`c` a density-independent width, :math:`A` the Dicke (motional-narrowing)
coefficient of the light molecules, :math:`x_p` the mole fraction of the collision
partner :math:`p` ("self" for the gas itself) and :math:`\gamma_p` its broadening
coefficient. Three forms of :math:`\gamma_p` are used, as the curated data allow:

- :class:`LinewidthLaw`: one coefficient for every line of a band, with a temperature
  law :math:`\gamma(T) = [\gamma_0 + c_1 (T - T_0) + c_2 (T - T_0)^2]
  (T_0 / T)^N` (:class:`Broadening`), and :math:`A(T) = A_0 (T / T_A)^{n_A}`;
- :class:`LineTable`: a coefficient per rotational level, tabulated at one temperature
  with a per-level exponent, interpolated in :math:`J`;
- :class:`ExponentialGapLaw`: the modified exponential-gap law of Koszykowski et al.
  (1987), in which a Q-branch line's width is the sum of the rates of rotational
  energy transfer out of its level.

Widths are in Hz, densities in amagat (:data:`AMAGAT`), temperatures in K. The laws are
per line; the overlap of lines at high density (line mixing, the collapse of a Q
branch) is not modelled.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
from jaxtyping import Float, Int

from optmat.conditions import ONE_ATM_0C
from optmat.constants import HBAR, K_B, STD_ATMOSPHERE, C
from optmat.exceptions import DataError

__all__ = [
    "AMAGAT",
    "WAVENUMBER_TO_HZ",
    "Broadening",
    "ExponentialGapLaw",
    "LineTable",
    "LineTableLaw",
    "LinewidthLaw",
]

AMAGAT: float = ONE_ATM_0C.pressure / (K_B * ONE_ATM_0C.temperature)
"""One amagat, m^-3: the ideal-gas number density at 1 atm and 0 degrees C."""

WAVENUMBER_TO_HZ: float = C * 100
"""Factor from a wavenumber in cm^-1 to a frequency in Hz."""


def _amagat(density: float) -> float:
    """Return a density in amagat, refusing a non-positive one.

    Parameters
    ----------
    density : float
        Number density, m^-3.

    Returns
    -------
    float
        The density in amagat.

    Raises
    ------
    DataError
        If the density is not positive.
    """
    if not density > 0:
        msg = f"a linewidth needs a positive density, got {density}"
        raise DataError(msg)
    return density / AMAGAT


@dataclass(frozen=True)
class Broadening:
    """A broadening coefficient and its temperature law.

    ``gamma(T) = [coefficient + linear (T - T0) + quadratic (T - T0)^2]
    (T0 / T)^exponent``
    in Hz/amagat (FWHM).

    Attributes
    ----------
    coefficient : float
        ``gamma(T0)``, Hz/amagat.
    reference_temperature : float
        ``T0``, K.
    exponent : float
        The power-law exponent ``N``.
    linear, quadratic : float
        Polynomial terms, Hz/amagat/K and Hz/amagat/K^2 (Bischel & Dyer's law).
    """

    coefficient: float
    reference_temperature: float
    exponent: float = 0.0
    linear: float = 0.0
    quadratic: float = 0.0

    def at(self, temperature: float) -> float:
        """Return the coefficient at a temperature.

        Parameters
        ----------
        temperature : float
            Temperature, K, > 0.

        Returns
        -------
        float
            ``gamma(T)``, Hz/amagat.
        """
        dt = temperature - self.reference_temperature
        polynomial = self.coefficient + self.linear * dt + self.quadratic * dt**2
        return polynomial * (self.reference_temperature / temperature) ** self.exponent


def _perturbers(
    fractions: Mapping[str, float] | None,
    known: Mapping[str, object],
    substitutes: Mapping[str, str],
) -> list[tuple[str, float]]:
    """Resolve the mole fractions to the perturbers that have coefficients.

    Parameters
    ----------
    fractions : mapping of str to float or None
        Mole fractions by perturber; ``{"self": 1}`` by default.
    known : mapping
        The perturbers with coefficients.
    substitutes : mapping of str to str
        Perturbers without coefficients, and the one whose coefficient stands in.

    Returns
    -------
    list of tuple of (str, float)
        The perturbers to use and their fractions.

    Raises
    ------
    DataError
        If a perturber has neither a coefficient nor a substitute.
    """
    resolved = []
    for perturber, fraction in (fractions or {"self": 1.0}).items():
        if fraction == 0:
            continue
        name = perturber if perturber in known else substitutes.get(perturber)
        if name is None or name not in known:
            msg = f"no broadening coefficient for perturber {perturber!r}"
            raise DataError(msg)
        resolved.append((name, fraction))
    return resolved


@dataclass(frozen=True)
class LinewidthLaw:
    """One width law for every line of a band.

    Attributes
    ----------
    constant : float
        Density-independent width, Hz.
    dicke : float
        Dicke coefficient ``A`` at ``dicke_temperature``, Hz amagat.
    dicke_temperature : float
        The temperature of ``dicke``, K.
    dicke_exponent : float
        ``n_A`` in ``A(T) = A (T / T_A)^n_A``.
    broadening : dict of str to Broadening
        By perturber (``self``, a gas id, or ``air``).
    substitutes : dict of str to str
        Perturbers without coefficients, and the perturber whose coefficient stands in
        for each (a declared approximation).
    sources : tuple of str
        The curated datasets and formulas the law comes from.
    """

    constant: float = 0.0
    dicke: float = 0.0
    dicke_temperature: float = 298.0
    dicke_exponent: float = 0.0
    broadening: Mapping[str, Broadening] = field(default_factory=dict)
    substitutes: Mapping[str, str] = field(default_factory=dict)
    sources: tuple[str, ...] = ()

    def fwhm(
        self,
        j: Int[np.ndarray, " n"],
        temperature: float,
        density: float,
        fractions: Mapping[str, float] | None = None,
    ) -> Float[np.ndarray, " n"]:
        """Return the linewidths at a state.

        Parameters
        ----------
        j : numpy.ndarray
            The lines' rotational levels (the width is the same for all).
        temperature : float
            Temperature, K, > 0.
        density : float
            Total number density, m^-3, > 0.
        fractions : mapping of str to float, optional
            Mole fractions by perturber; the pure gas by default.

        Returns
        -------
        numpy.ndarray
            FWHM of each line, Hz.
        """
        amagat = _amagat(density)
        dicke = (
            self.dicke * (temperature / self.dicke_temperature) ** self.dicke_exponent
        )
        width = self.constant + dicke / amagat
        if self.broadening:  # a law with no collisional term takes any partner
            for name, fraction in _perturbers(
                fractions, self.broadening, self.substitutes
            ):
                width += fraction * amagat * self.broadening[name].at(temperature)
        return np.full(np.shape(j), width, dtype=float)


@dataclass(frozen=True)
class LineTable:
    """Broadening coefficients per rotational level, from one tabulation.

    Between tabulated levels the coefficient and exponent are interpolated linearly in
    ``J``; beyond the ends they are held at the end values (the linewidths of the
    highest levels change slowly, and their populations are small).

    Attributes
    ----------
    j : tuple of int
        The tabulated levels.
    coefficient : tuple of float
        ``gamma(J, T0)``, Hz/amagat.
    reference_temperature : float
        ``T0``, K.
    exponent : tuple of float
        The power-law exponent of each level (``(T0 / T)^N``).
    """

    j: tuple[int, ...]
    coefficient: tuple[float, ...]
    reference_temperature: float
    exponent: tuple[float, ...]

    def __post_init__(self) -> None:
        """Check that the columns agree and the levels increase.

        Raises
        ------
        DataError
            If they do not.
        """
        if not len(self.j) == len(self.coefficient) == len(self.exponent) > 0:
            msg = "a line table needs one coefficient and one exponent per level"
            raise DataError(msg)
        if np.any(np.diff(self.j) <= 0):
            msg = "a line table's levels must increase"
            raise DataError(msg)

    def at(
        self, j: Int[np.ndarray, " n"], temperature: float
    ) -> Float[np.ndarray, " n"]:
        """Return the coefficients of levels ``j`` at a temperature.

        Parameters
        ----------
        j : numpy.ndarray
            Rotational levels.
        temperature : float
            Temperature, K, > 0.

        Returns
        -------
        numpy.ndarray
            ``gamma(J, T)``, Hz/amagat.
        """
        levels = np.asarray(j, dtype=float)
        gamma = np.interp(levels, self.j, self.coefficient)
        exponent = np.interp(levels, self.j, self.exponent)
        return gamma * (self.reference_temperature / temperature) ** exponent


@dataclass(frozen=True)
class LineTableLaw:
    """Per-line widths from per-level tables, one per perturber.

    Attributes
    ----------
    tables : dict of str to LineTable
        By perturber.
    substitutes : dict of str to str
        Perturbers without a table, and the perturber whose table stands in.
    dicke : float
        A Dicke coefficient ``A`` for every line, Hz amagat (H2's rotational lines).
    sources : tuple of str
        The curated datasets the tables come from.
    """

    tables: Mapping[str, LineTable]
    substitutes: Mapping[str, str] = field(default_factory=dict)
    dicke: float = 0.0
    sources: tuple[str, ...] = ()

    def fwhm(
        self,
        j: Int[np.ndarray, " n"],
        temperature: float,
        density: float,
        fractions: Mapping[str, float] | None = None,
    ) -> Float[np.ndarray, " n"]:
        """Return the linewidths at a state.

        Parameters
        ----------
        j : numpy.ndarray
            The lines' lower rotational levels.
        temperature : float
            Temperature, K, > 0.
        density : float
            Total number density, m^-3, > 0.
        fractions : mapping of str to float, optional
            Mole fractions by perturber; the pure gas by default.

        Returns
        -------
        numpy.ndarray
            FWHM of each line, Hz.
        """
        amagat = _amagat(density)
        width = np.full(np.shape(j), self.dicke / amagat)
        for name, fraction in _perturbers(fractions, self.tables, self.substitutes):
            width = width + fraction * amagat * self.tables[name].at(j, temperature)
        return width


@dataclass(frozen=True)
class ExponentialGapLaw:
    """Q-branch widths from the modified exponential-gap law (Koszykowski et al. 1987).

    The rate of transfer up from level ``i`` to ``j`` (``E_j > E_i``, the same nuclear
    spin symmetry) is ``gamma_ji = alpha p [(1 + a E_i / (k T delta)) / (1 + a E_i /
    (k T))]^2 exp(-beta (E_j - E_i) / (k T))``, the rate down follows by detailed
    balance, ``gamma_ij = (2 J_i + 1) / (2 J_j + 1) gamma_ji exp((E_j - E_i) / (k T))``,
    and the FWHM of line ``J`` is ``2 sum_{i != J} gamma_iJ`` (their eqs. 4.1-4.4).
    ``alpha`` scales with temperature as ``(T0 / T)^n``; ``p`` is the pressure in atm.

    Attributes
    ----------
    alpha : float
        ``alpha`` at ``T0``, cm^-1/atm.
    beta, delta, a : float
        The law's constants.
    exponent : float
        ``n``.
    reference_temperature : float
        ``T0``, K.
    b, d : float
        The rotational constants of the levels, cm^-1.
    j_max : int
        The highest level summed over.
    foreign_as_self : tuple of str
        Collision partners whose broadening is taken as the self-broadening (a declared
        approximation, where no foreign coefficients are curated).
    sources : tuple of str
        The formula record and the constants' datasets.
    """

    alpha: float
    beta: float
    delta: float
    a: float
    exponent: float
    reference_temperature: float
    b: float
    d: float = 0.0
    j_max: int = 60
    foreign_as_self: tuple[str, ...] = ()
    sources: tuple[str, ...] = ()

    def fwhm(
        self,
        j: Int[np.ndarray, " n"],
        temperature: float,
        density: float,
        fractions: Mapping[str, float] | None = None,
    ) -> Float[np.ndarray, " n"]:
        """Return the linewidths at a state.

        Parameters
        ----------
        j : numpy.ndarray
            The lines' rotational levels.
        temperature : float
            Temperature, K, > 0.
        density : float
            Total number density, m^-3, > 0 (converted to a pressure with the ideal
            gas law).
        fractions : mapping of str to float, optional
            The pure gas, or partners listed in :attr:`foreign_as_self`: the law is for
            self-broadening.

        Returns
        -------
        numpy.ndarray
            FWHM of each line, Hz.

        Raises
        ------
        DataError
            If a foreign perturber is asked for.
        """
        allowed = {"self", *self.foreign_as_self}
        if fractions and any(k not in allowed and v for k, v in fractions.items()):
            msg = "the exponential-gap law is for self-broadening only"
            raise DataError(msg)
        pressure = _amagat(density) * AMAGAT * K_B * temperature / STD_ATMOSPHERE
        levels = np.arange(self.j_max + 1)
        x = levels * (levels + 1.0)
        energy = self.b * x - self.d * x**2  # cm^-1
        kt = K_B * temperature / (2 * math.pi * HBAR * C * 100)  # cm^-1
        alpha = self.alpha * (self.reference_temperature / temperature) ** self.exponent
        widths = []
        for level in np.asarray(j, dtype=int):
            total = 0.0
            for other in range(level % 2, self.j_max + 1, 2):
                if other == level:
                    continue
                low, high = min(level, other), max(level, other)
                gap = energy[high] - energy[low]
                factor = (1 + self.a * energy[low] / (kt * self.delta)) / (
                    1 + self.a * energy[low] / kt
                )
                up = alpha * pressure * factor**2 * math.exp(-self.beta * gap / kt)
                if other > level:
                    total += up
                else:
                    total += (2 * low + 1) / (2 * high + 1) * up * math.exp(gap / kt)
            widths.append(2 * total * WAVENUMBER_TO_HZ)
        return np.array(widths)

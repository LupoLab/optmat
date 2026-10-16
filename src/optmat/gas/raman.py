# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The Raman lines of the molecular gases, ready for a propagation code (TC-439).

A molecule's delayed (Raman) third-order response is a sum of damped oscillators, one
per Raman line, each driven by the intensity. This module turns curated spectroscopic
data into those lines: the frequency of every line, its coupling amplitude and its
linewidth at a stated temperature and density. The tensor response is assembled by
:mod:`optmat.gas.raman_response`; a propagation
code chooses how to sample or integrate it.

**The convention of the amplitudes.** For a field linearly polarised along
:math:`x`, the Raman polarisation of a gas of number density :math:`N` is

.. math::

   P_x(t) = \varepsilon_0 N E_x(t) \int_0^\infty h(\tau) E_x^2(t - \tau)\, d\tau,
   \qquad h(\tau) = \sum_k A_k \sin(\Omega_k \tau)\, e^{-\Gamma_k \tau},

with :math:`A_k` in m^5 V^-2 s^-1, :math:`\Omega_k` the line's angular frequency and
:math:`\Gamma_k = \pi \Delta\nu_k` its dephasing rate (:math:`\Delta\nu_k` the full
width at half maximum in Hz). A propagation adapter maps the
lines one to one.

- **Rotational S branch** (:math:`J \to J + 2` of a linear molecule):

  .. math::

     A_J = \frac{2}{15} \frac{\Delta\alpha^2}{\hbar \varepsilon_0}
           \frac{(J + 1)(J + 2)}{2J + 3} (n_J - n_{J+2}),

  with :math:`\Delta\alpha` the polarizability anisotropy (C m^2/V) and
  :math:`n_J = g_J e^{-E_J / k_B T} / Z` the population of *one* of the
  :math:`2J + 1` sublevels (Placzek-Teller; equivalently,
  :math:`p_J/(2J+1) = n_J`).
- **Vibrational Q branch** (:math:`v = 0 \to 1`, one line per rotational level):

  .. math::

     A_J = \frac{|\langle 0 | \bar\alpha | 1 \rangle|^2}{\hbar \varepsilon_0} P_J,

  with :math:`\langle 0 | \bar\alpha | 1 \rangle` the transition polarizability of the
  mean (isotropic) polarizability (C m^2/V) and :math:`P_J` the fraction of molecules
  in level :math:`J` (all sublevels). For a harmonic mode
  :math:`|\langle 0|\bar\alpha|1\rangle|^2 = (\partial\bar\alpha/\partial q)^2 \hbar /
  (2 \mu \Omega)`, so :math:`A = (\partial\bar\alpha/\partial q)^2 / (2\mu\Omega
  \varepsilon_0)`, with a consistent vibrational coordinate and reduced mass.

Energies come from the curated constants in cm^-1: the rotational levels
:math:`E_J = B J(J+1) - D [J(J+1)]^2 + H [J(J+1)]^3`, and the Q-branch line
:math:`\nu_J = \nu_0 + (B_1 - B_0) J(J+1) - (D_1 - D_0) [J(J+1)]^2`.
"""

import dataclasses
import functools
import math
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from jaxtyping import Float, Int

from optmat import _catalogue
from optmat._catalogue import records as sources_of
from optmat.conditions import ONE_ATM_0C
from optmat.constants import ATOMIC_MASS, EPS0, HBAR, K_B, STD_ATMOSPHERE, C
from optmat.datasets import Dataset
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.gas.linewidths import (
    WAVENUMBER_TO_HZ,
    Broadening,
    ExponentialGapLaw,
    LineTable,
    LineTableLaw,
    LinewidthLaw,
)
from optmat.gas.reduction import POLARIZABILITY_UNITS

__all__ = [
    "RAMAN_FILE",
    "AnisotropyTable",
    "BandLines",
    "RamanLines",
    "RamanParameters",
    "RotationalBand",
    "RotationalConstants",
    "SpinStatistics",
    "VibrationalBand",
    "band_lines",
    "default_j_max",
    "raman_lines",
    "raman_parameters",
    "rotational_lines",
    "rotational_nonlinear_index",
    "rotational_populations",
    "spin_statistics",
    "stokes_gain",
    "vibrational_lines",
]

_WAVENUMBER_TO_ANGULAR_FREQUENCY: Final[float] = 2 * math.pi * C * 100
"""Factor from a wavenumber in cm^-1 to an angular frequency in rad/s."""


@dataclass(frozen=True)
class SpinStatistics:
    """The nuclear-spin weights of the even and odd rotational levels.

    Attributes
    ----------
    even : int
        Statistical weight of the even-J (or even-N) levels.
    odd : int
        Statistical weight of the odd levels.
    reason : str
        Where the weights come from.
    """

    even: int
    odd: int
    reason: str

    def weight(self, j: Int[np.ndarray, " n"]) -> Int[np.ndarray, " n"]:
        """Return the weight of each rotational level.

        Parameters
        ----------
        j : numpy.ndarray
            Rotational quantum numbers.

        Returns
        -------
        numpy.ndarray
            ``even`` or ``odd`` for each.
        """
        return np.where(np.asarray(j) % 2 == 0, self.even, self.odd)


def _homonuclear(spin_twice: int, electronic_plus: bool, reason: str) -> SpinStatistics:
    """Return the weights of a homonuclear or centrosymmetric linear molecule.

    For two identical nuclei of spin ``I`` there are ``(I + 1)(2I + 1)`` symmetric and
    ``I (2I + 1)`` antisymmetric spin states. The total wave function must be symmetric
    under exchange of bosons and antisymmetric for fermions; in a Sigma_g^+ electronic
    state the rotational factor has the parity of J, and a Sigma_g^- state reverses it.

    Parameters
    ----------
    spin_twice : int
        ``2 I`` of the exchanged nuclei.
    electronic_plus : bool
        True for a Sigma_g^+ ground state, False for Sigma_g^- (O2).
    reason : str
        The source of the nuclear spin and the electronic state.

    Returns
    -------
    SpinStatistics
        The weights of the even and odd levels.
    """
    symmetric = (spin_twice + 2) * (spin_twice + 1) // 2
    antisymmetric = spin_twice * (spin_twice + 1) // 2
    boson = spin_twice % 2 == 0
    even, odd = (symmetric, antisymmetric) if boson else (antisymmetric, symmetric)
    if not electronic_plus:
        even, odd = odd, even
    return SpinStatistics(even, odd, reason)


# The nuclear spins are I(1H) = 1/2, I(2H) = 1, I(14N) = 1, I(16O) = 0 (quoted from
# memory, standard nuclear data; D2's is Murphy & Johnston's, murphy_1934, TC-173), and
# the ground states are Sigma_g^+ except O2's Sigma_g^- (TC-442).
_SPIN: Final[dict[Gas, SpinStatistics]] = {
    Gas.H2: _homonuclear(1, electronic_plus=True, reason="I(1H) = 1/2, X 1Sigma_g+"),
    Gas.D2: _homonuclear(2, electronic_plus=True, reason="I(2H) = 1 (murphy_1934)"),
    Gas.N2: _homonuclear(2, electronic_plus=True, reason="I(14N) = 1, X 1Sigma_g+"),
    Gas.O2: _homonuclear(0, electronic_plus=False, reason="I(16O) = 0, X 3Sigma_g-"),
    Gas.CO2: _homonuclear(0, electronic_plus=True, reason="I(16O) = 0, X 1Sigma_g+"),
    Gas.N2O: SpinStatistics(1, 1, "N-N-O has no exchange symmetry"),
    Gas.CO: SpinStatistics(1, 1, "heteronuclear"),
}


def spin_statistics(gas: Gas | str) -> SpinStatistics:
    """Return the nuclear-spin weights of a linear molecule's rotational levels.

    Parameters
    ----------
    gas : Gas or str
        A linear molecule.

    Returns
    -------
    SpinStatistics
        The weights.

    Raises
    ------
    DataError
        If the gas is not a linear molecule (an atom, or a spherical top, which has no
        rotational Raman spectrum).

    Examples
    --------
    >>> weights = spin_statistics("N2")
    >>> weights.even, weights.odd
    (6, 3)
    >>> weights = spin_statistics("O2")
    >>> weights.even, weights.odd
    (0, 1)
    """
    species = as_gas(gas)
    if species not in _SPIN:
        msg = f"{species.value} is not a linear molecule: no rotational Raman spectrum"
        raise DataError(msg)
    return _SPIN[species]


@dataclass(frozen=True)
class RotationalConstants:
    """The constants of one vibrational level's rotational ladder, cm^-1.

    Attributes
    ----------
    b : float
        Rotational constant B.
    d : float
        Centrifugal distortion D.
    h : float
        Sextic distortion H.
    """

    b: float
    d: float = 0.0
    h: float = 0.0

    def energy(self, j: Int[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
        """Return the rotational energies of levels J.

        Parameters
        ----------
        j : numpy.ndarray
            Rotational quantum numbers.

        Returns
        -------
        numpy.ndarray
            ``B x - D x^2 + H x^3`` with ``x = J (J + 1)``, cm^-1.
        """
        x = np.asarray(j, dtype=float) * (np.asarray(j, dtype=float) + 1)
        return self.b * x - self.d * x**2 + self.h * x**3


def default_j_max(constants: RotationalConstants, temperature: float) -> int:
    """Return the rotational level at about 50 k_B T, where the population is ~1e-20.

    Parameters
    ----------
    constants : RotationalConstants
        The rotational constants, cm^-1.
    temperature : float
        Temperature, K.

    Returns
    -------
    int
        ``sqrt(50 k_B T / B) + 5``, at least 10.

    Examples
    --------
    >>> default_j_max(RotationalConstants(2.0), 295.0)
    76
    """
    kt = K_B * temperature / (HBAR * _WAVENUMBER_TO_ANGULAR_FREQUENCY)
    return max(10, int(math.sqrt(50 * kt / constants.b)) + 5)


def rotational_populations(
    constants: RotationalConstants,
    spin: SpinStatistics,
    temperature: float,
    j_max: int,
) -> tuple[Int[np.ndarray, " n"], Float[np.ndarray, " n"]]:
    """Return the fraction of molecules in each rotational level.

    A truncated distortion series (``B x - D x^2``) turns over at high J, where it no
    longer describes the molecule; the ladder is cut at its last increasing level,
    which the callers then require to be negligibly populated.

    Parameters
    ----------
    constants : RotationalConstants
        The ground vibrational level's constants.
    spin : SpinStatistics
        The nuclear-spin weights.
    temperature : float
        Temperature, K, > 0.
    j_max : int
        The highest level wanted.

    Returns
    -------
    tuple of numpy.ndarray
        ``J = 0 ... min(j_max, the turning level)`` and the fraction ``P_J`` of
        molecules in each (all ``2J + 1`` sublevels), summing to one.
    """
    j = np.arange(j_max + 1)
    energy = constants.energy(j)
    falling = np.flatnonzero(np.diff(energy) <= 0)
    if falling.size:
        j, energy = j[: falling[0] + 1], energy[: falling[0] + 1]
    energy_j = energy * _WAVENUMBER_TO_ANGULAR_FREQUENCY * HBAR
    weights = spin.weight(j) * (2 * j + 1) * np.exp(-energy_j / (K_B * temperature))
    return j, weights / np.sum(weights)


def _check_negligible(
    j: Int[np.ndarray, " n"], population: Float[np.ndarray, " n"], floor: float
) -> None:
    """Refuse a ladder whose highest kept levels are populated.

    Parameters
    ----------
    j : numpy.ndarray
        The levels kept.
    population : numpy.ndarray
        Their populations.
    floor : float
        The population regarded as negligible.

    Raises
    ------
    DataError
        If either of the two highest levels (one of each parity) holds more than
        ``floor``.
    """
    if np.max(population[-2:]) > floor:
        msg = (
            f"the population at J = {int(j[-1])} is {np.max(population[-2:]):.1e}: "
            "raise j_max, or the distortion series ends before the population does"
        )
        raise DataError(msg)


@dataclass(frozen=True)
class RamanLines:
    """Raman lines in the convention of the module docstring.

    Attributes
    ----------
    branch : tuple of str
        ``S`` (rotational) or ``Q`` (vibrational), per line.
    j : numpy.ndarray
        The lower level's rotational quantum number (-1 where the line is not resolved
        by J).
    frequency : numpy.ndarray
        Angular frequency ``Omega``, rad/s.
    amplitude : numpy.ndarray
        ``A``, m^5 V^-2 s^-1.
    fwhm : numpy.ndarray or None
        The full width at half maximum of each line at the state asked for, Hz
        (``Gamma = pi FWHM``); None where no state was given.
    """

    branch: tuple[str, ...]
    j: Int[np.ndarray, " n"]
    frequency: Float[np.ndarray, " n"]
    amplitude: Float[np.ndarray, " n"]
    fwhm: Float[np.ndarray, " n"] | None = None

    def __len__(self) -> int:
        """Return the number of lines.

        Returns
        -------
        int
            The number of lines.
        """
        return len(self.branch)


@dataclass(frozen=True)
class AnisotropyTable:
    """A polarizability element that depends on the rotational level.

    For H2 the anisotropy of the J -> J + 2 transition rises by 10 % from J = 0 to
    J = 6 (centrifugal stretching), and the Q(J) transition polarizability by 2 % from
    J = 0 to 5. Between tabulated levels it is interpolated linearly; beyond them it is
    held at the end values.

    Attributes
    ----------
    j : tuple of int
        The lower levels, increasing.
    value : tuple of float
        The anisotropy of each, C m^2/V.
    """

    j: tuple[int, ...]
    value: tuple[float, ...]

    def __call__(self, j: Int[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
        """Return the anisotropy of the transitions from levels ``j``.

        Parameters
        ----------
        j : numpy.ndarray
            Lower rotational levels.

        Returns
        -------
        numpy.ndarray
            The anisotropy, C m^2/V.
        """
        return np.interp(np.asarray(j, dtype=float), self.j, self.value)


def rotational_lines(
    constants: RotationalConstants,
    spin: SpinStatistics,
    anisotropy: float | AnisotropyTable,
    temperature: float,
    *,
    j_max: int | None = None,
    population_floor: float = 1e-12,
) -> RamanLines:
    """Return the rotational S-branch lines of a linear molecule.

    Parameters
    ----------
    constants : RotationalConstants
        The ground vibrational level's constants, cm^-1.
    spin : SpinStatistics
        The nuclear-spin weights.
    anisotropy : float or AnisotropyTable
        Polarizability anisotropy ``Delta alpha``, C m^2/V, or its value per lower
        level.
    temperature : float
        Temperature, K, > 0.
    j_max : int, optional
        The highest lower level considered; by default the level near 50 k_B T
        (:func:`default_j_max`).
    population_floor : float, optional
        Lines whose lower and upper populations are both below this fraction are
        dropped.

    Returns
    -------
    RamanLines
        One line per ``J -> J + 2`` with a nonzero population difference.

    Raises
    ------
    DataError
        If the population at the highest level kept is not negligible.

    Examples
    --------
    With illustrative constants close to nitrogen's (B = 2 cm^-1 and a
    polarizability anisotropy of 4.6 au), the line of largest amplitude at 295 K
    starts from J = 8 (an even level, which nitrogen's spin statistics favour 2:1):

    >>> from optmat.constants import AU_POLARIZABILITY
    >>> anisotropy = 4.6 * AU_POLARIZABILITY
    >>> lines = rotational_lines(
    ...     RotationalConstants(2.0), spin_statistics("N2"), anisotropy, 295.0
    ... )
    >>> int(lines.j[np.argmax(lines.amplitude)])
    8
    """
    if j_max is None:
        j_max = default_j_max(constants, temperature)
    j, population = rotational_populations(constants, spin, temperature, j_max + 2)
    _check_negligible(j, population, population_floor)
    per_sublevel = population / (2 * j + 1)
    lower = j[:-2]
    difference = per_sublevel[:-2] - per_sublevel[2:]
    frequency = (
        constants.energy(lower + 2) - constants.energy(lower)
    ) * _WAVENUMBER_TO_ANGULAR_FREQUENCY
    delta = anisotropy(lower) if isinstance(anisotropy, AnisotropyTable) else anisotropy
    amplitude = (
        2
        / 15
        * delta**2
        / (HBAR * EPS0)
        * (lower + 1)
        * (lower + 2)
        / (2 * lower + 3)
        * difference
    )
    keep = (population[:-2] > population_floor) | (population[2:] > population_floor)
    keep &= difference != 0
    return RamanLines(
        branch=("S",) * int(np.sum(keep)),
        j=lower[keep],
        frequency=frequency[keep],
        amplitude=amplitude[keep],
    )


def vibrational_lines(  # noqa: PLR0913 - a band, its two ladders and the state
    band_origin: float,
    transition_polarizability: float | AnisotropyTable,
    *,
    lower: RotationalConstants | None = None,
    upper: RotationalConstants | None = None,
    spin: SpinStatistics | None = None,
    temperature: float = 296.0,
    j_max: int | None = None,
    population_floor: float = 1e-12,
) -> RamanLines:
    """Return the vibrational Q-branch lines of a mode.

    Parameters
    ----------
    band_origin : float
        The band origin ``nu_0``, cm^-1.
    transition_polarizability : float or AnisotropyTable
        ``<0| alpha_mean |1>``, C m^2/V, or its value per rotational level (the
        table's ``J`` is the level of the Q(J) line).
    lower, upper : RotationalConstants, optional
        The rotational constants of ``v = 0`` and ``v = 1``. Without them (a spherical
        top, or an unresolved band) the branch is one line at the band origin.
    spin : SpinStatistics, optional
        The nuclear-spin weights; needed with the constants.
    temperature : float, optional
        Temperature, K.
    j_max : int, optional
        The highest rotational level considered; by default the level near 50 k_B T
        (:func:`default_j_max`).
    population_floor : float, optional
        Levels below this fraction are dropped.

    Returns
    -------
    RamanLines
        The Q(J) lines, or one line.

    Raises
    ------
    DataError
        If constants are given without spin weights, or the population at the highest
        level kept is not negligible.

    Examples
    --------
    With illustrative constants close to hydrogen's (B = 60 cm^-1 in v = 0 and 57
    cm^-1 in v = 1), Q(1) lies 2 (B_1 - B_0) = -6 cm^-1 from Q(0):

    >>> lines = vibrational_lines(
    ...     4160.0, 1e-41, lower=RotationalConstants(60.0),
    ...     upper=RotationalConstants(57.0), spin=spin_statistics("H2"))
    >>> f"{(lines.frequency[1] - lines.frequency[0]) / (2 * np.pi * C * 100):.1f}"
    '-6.0'
    """
    if lower is None or upper is None:
        element = (
            float(np.mean(transition_polarizability.value))
            if isinstance(transition_polarizability, AnisotropyTable)
            else transition_polarizability
        )
        return RamanLines(
            branch=("Q",),
            j=np.array([-1]),
            frequency=np.array([band_origin * _WAVENUMBER_TO_ANGULAR_FREQUENCY]),
            amplitude=np.array([element**2 / (HBAR * EPS0)]),
        )
    if spin is None:
        msg = "a rotationally resolved Q branch needs the spin weights"
        raise DataError(msg)
    if j_max is None:
        j_max = default_j_max(lower, temperature)
    j, population = rotational_populations(lower, spin, temperature, j_max)
    _check_negligible(j, population, population_floor)
    keep = population > population_floor
    shift = upper.energy(j) - lower.energy(j)
    element = (
        transition_polarizability(j)
        if isinstance(transition_polarizability, AnisotropyTable)
        else np.full(len(j), transition_polarizability)
    )
    strength = element**2 / (HBAR * EPS0)
    return RamanLines(
        branch=("Q",) * int(np.sum(keep)),
        j=j[keep],
        frequency=(band_origin + shift[keep]) * _WAVENUMBER_TO_ANGULAR_FREQUENCY,
        amplitude=strength[keep] * population[keep],
    )


# The curated Raman parameters of each gas (TC-492): data/gas/raman.toml names, per gas
# and band, the datasets and formulas each parameter comes from, with the reason.

RAMAN_FILE = GAS_DATA_DIRECTORY / "raman.toml"
"""The per-gas choices of the Raman parameters, each with its reason."""

_LinewidthModel = LinewidthLaw | LineTableLaw | ExponentialGapLaw
"""A width law of a band: one law, per-line tables, or the exponential-gap law."""

_PRESSURES_ATM: Final[dict[str, float]] = {"atm": 1.0, "bar": 1e5 / STD_ATMOSPHERE}


@dataclass(frozen=True)
class RotationalBand:
    """The rotational Raman band of a linear molecule.

    Attributes
    ----------
    constants : RotationalConstants
        The ground vibrational level's constants, cm^-1.
    spin : SpinStatistics
        The nuclear-spin weights.
    anisotropy : float or AnisotropyTable
        The polarizability anisotropy, C m^2/V.
    widths : LinewidthLaw or LineTableLaw or ExponentialGapLaw
        The width law of the lines.
    sources : dict of str to str
        The dataset or formula behind each parameter.
    """

    constants: RotationalConstants
    spin: SpinStatistics
    anisotropy: float | AnisotropyTable
    widths: _LinewidthModel
    sources: Mapping[str, str]


@dataclass(frozen=True)
class VibrationalBand:
    """One vibrational Raman band (its Q branch).

    Attributes
    ----------
    mode : str
        The mode, as the curated data describe it.
    band_origin : float
        ``nu_0``, cm^-1.
    transition_polarizability : float
        ``<0| alpha_mean |1>``, C m^2/V.
    lower, upper : RotationalConstants or None
        The rotational constants of the two vibrational levels (None for an
        unresolved band).
    widths : LinewidthLaw or LineTableLaw or ExponentialGapLaw
        The width law of the lines.
    sources : dict of str to str
        The dataset or formula behind each parameter.
    resonance : float or None
        The effective resonance ``nu_i`` of a single-resonance dispersion law,
        cm^-1: the element scales as ``1 / (nu_i^2 - nu_p^2)`` with the pump
        wavenumber ``nu_p`` (Bischel & Black 1983). None where no law is curated.
    reference_wavenumber : float
        The pump wavenumber at which ``transition_polarizability`` holds, cm^-1 (0
        for a static value).
    transition_anisotropy : AnisotropyTable or None
        ``<0J| alpha_par - alpha_perp |1J>`` per rotational level, C m^2/V, where the
        band's elements are curated as components (H2, D2); None otherwise.
    depolarization : float or None
        The depolarization ratio of the Q branch for linearly polarised light at 90
        degrees, ``rho = 3 g^2 / (45 a^2 + 4 g^2)``, where one is curated (the value
        that reduced the cross sections, or the band's ``depolarization`` key); None
        where none is.
    key : str
        A short, stable name of the band (``v1``, ``nu3``, ``fermi_1388``): the band of
        the Raman response's lines and of their identifiers; the mode where none is
        curated.
    """

    mode: str
    band_origin: float
    transition_polarizability: float | AnisotropyTable
    lower: RotationalConstants | None
    upper: RotationalConstants | None
    widths: _LinewidthModel
    sources: Mapping[str, str]
    resonance: float | None = None
    reference_wavenumber: float = 0.0
    transition_anisotropy: AnisotropyTable | None = None
    depolarization: float | None = None
    key: str = ""

    def dispersion_factor(self, wavelength: float | None) -> float:
        """Return the factor on the transition polarizability at a pump wavelength.

        Parameters
        ----------
        wavelength : float or None
            Vacuum wavelength of the pump, m; None for the element as curated.

        Returns
        -------
        float
            ``(nu_i^2 - nu_ref^2) / (nu_i^2 - nu_p^2)``, or 1 without a law or a
            wavelength.

        Raises
        ------
        DataError
            If the pump is at or beyond the effective resonance.
        """
        if wavelength is None or self.resonance is None:
            return 1.0
        pump = 1 / (wavelength * 100)
        if pump >= self.resonance:
            msg = f"the pump at {pump:.0f} cm-1 is beyond the resonance of {self.mode}"
            raise DataError(msg)
        return (self.resonance**2 - self.reference_wavenumber**2) / (
            self.resonance**2 - pump**2
        )


@dataclass(frozen=True)
class RamanParameters:
    """The curated Raman parameters of a gas.

    Attributes
    ----------
    gas : Gas
        The gas.
    rotational : RotationalBand or None
        The rotational band; None for a spherical top.
    vibrational : tuple of VibrationalBand
        The Raman-active vibrational bands curated.
    reasons : dict of str to str
        Why each choice was made (from ``raman.toml``).
    missing : str
        What is known to be missing, or an empty string.
    """

    gas: Gas
    rotational: RotationalBand | None
    vibrational: tuple[VibrationalBand, ...]
    reasons: Mapping[str, str] = field(default_factory=dict)
    missing: str = ""


def _table(value: object, name: str) -> dict[str, object]:
    """Return a TOML sub-table with string keys.

    Parameters
    ----------
    value : object
        The value read from ``raman.toml``.
    name : str
        Its key, for the error message.

    Returns
    -------
    dict of str to object
        The table.

    Raises
    ------
    DataError
        If the value is not a table.
    """
    if not isinstance(value, Mapping):
        msg = f"raman.toml: {name} must be a table"
        raise DataError(msg)
    return {str(k): v for k, v in value.items()}  # pyright: ignore[reportUnknownVariableType]


def _strings(value: object, name: str) -> dict[str, str]:
    """Return a TOML sub-table of strings.

    Parameters
    ----------
    value : object
        The value read from ``raman.toml``.
    name : str
        Its key, for the error message.

    Returns
    -------
    dict of str to str
        The table's keys and values as strings.
    """
    return {k: str(v) for k, v in _table(value, name).items()}


def _list(value: object, name: str) -> list[object]:
    """Return a TOML array.

    Parameters
    ----------
    value : object
        The value read from ``raman.toml``.
    name : str
        Its key, for the error message.

    Returns
    -------
    list of object
        The array.

    Raises
    ------
    DataError
        If the value is not an array.
    """
    if not isinstance(value, list):
        msg = f"raman.toml: {name} must be an array"
        raise DataError(msg)
    return list(value)  # pyright: ignore[reportUnknownArgumentType]


def _get(by_key: Mapping[str, Dataset], reference: str) -> tuple[Dataset, int | None]:
    """Return a dataset and a row index from ``key`` or ``key#row``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        ``source/dataset`` or ``source/dataset#row``.

    Returns
    -------
    tuple of (Dataset, int or None)
        The dataset and the row's index (None when no row is named).

    Raises
    ------
    DataError
        If the dataset or row does not exist.
    """
    key, _, row = reference.partition("#")
    if key not in by_key:
        msg = f"raman.toml names an unknown dataset {key!r}"
        raise DataError(msg)
    dataset = by_key[key]
    if not row:
        return dataset, None
    if row not in dataset.rows:
        msg = f"{key} has no row {row!r}"
        raise DataError(msg)
    return dataset, dataset.rows.index(row)


def _scalar(by_key: Mapping[str, Dataset], reference: str) -> float:
    """Return one value of a dataset, times its scale.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        ``source/dataset`` (of one row) or ``source/dataset#row``.

    Returns
    -------
    float
        The value, in the dataset's unit.

    Raises
    ------
    DataError
        If a dataset of several rows is named without a row.
    """
    dataset, index = _get(by_key, reference)
    if index is None:
        if len(dataset) != 1:
            msg = f"{dataset.key} has {len(dataset)} rows: name one as #<row>"
            raise DataError(msg)
        index = 0
    return float(dataset.scaled_values()[index])


def _wavenumber(by_key: Mapping[str, Dataset], reference: str) -> float:
    """Return a constant or a line position in cm^-1.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        The dataset (and row).

    Returns
    -------
    float
        The value, cm^-1.

    Raises
    ------
    DataError
        If its unit is not a wavenumber or a frequency.
    """
    dataset, _ = _get(by_key, reference)
    factors = {
        "cm-1": 1.0,
        "MHz": 1e6 / WAVENUMBER_TO_HZ,
        "GHz": 1e9 / WAVENUMBER_TO_HZ,
    }
    if dataset.unit not in factors:
        msg = f"{dataset.key}: {dataset.unit!r} is not a wavenumber"
        raise DataError(msg)
    return _scalar(by_key, reference) * factors[dataset.unit]


def _polarizability(by_key: Mapping[str, Dataset], reference: str) -> float:
    """Return a polarizability (an anisotropy, a transition element) in C m^2/V.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        The dataset (and row).

    Returns
    -------
    float
        The value, C m^2/V.
    """
    dataset, _ = _get(by_key, reference)
    return _scalar(by_key, reference) * POLARIZABILITY_UNITS[dataset.unit]


def _transition_from_derivative(
    by_key: Mapping[str, Dataset], reference: str, band_origin: float
) -> float:
    """Return a transition polarizability from a normal-coordinate derivative.

    For a harmonic mode with the mass-weighted normal coordinate ``Q``,
    ``|<0|Q|1>|^2 = hbar / (2 omega)``, so ``<0|alpha|1> = (d alpha / d Q)
    (hbar / (2 omega))^(1/2)``: no reduced mass or bond convention enters.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        A ``polarizability_derivative`` of the mean polarizability with respect to the
        mass-weighted normal coordinate, in angstrom2/amu^(1/2) (a polarizability
        volume per unit coordinate).
    band_origin : float
        The mode's wavenumber, cm^-1.

    Returns
    -------
    float
        ``<0| alpha_mean |1>``, C m^2/V.

    Raises
    ------
    DataError
        If the derivative is not of that kind.
    """
    dataset, _ = _get(by_key, reference)
    if (
        dataset.unit != "angstrom2/amu^(1/2)"
        or dataset.parameters.get("component") != "mean"
    ):
        msg = f"{dataset.key} is not a mean derivative in angstrom2/amu^(1/2)"
        raise DataError(msg)
    derivative = (
        _scalar(by_key, reference)
        * POLARIZABILITY_UNITS["angstrom3"]
        * 1e10  # angstrom2 = angstrom3 / angstrom
        / math.sqrt(ATOMIC_MASS)
    )
    omega = band_origin * _WAVENUMBER_TO_ANGULAR_FREQUENCY
    return derivative * math.sqrt(HBAR / (2 * omega))


def _mean_from_components(
    by_key: Mapping[str, Dataset], components: Mapping[str, str]
) -> AnisotropyTable:
    """Return the mean transition polarizability per level from its components.

    ``<alpha_mean> = (<alpha_par> + 2 <alpha_perp>) / 3``, level by level.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    components : mapping of str to str
        ``parallel`` and ``perpendicular``: ``raman_transition_polarizability``
        datasets with a row per lower level ``J``.

    Returns
    -------
    AnisotropyTable
        The mean element per level, C m^2/V.

    Raises
    ------
    DataError
        If the two datasets do not cover the same levels.
    """
    levels, parallel, perpendicular, unit = _components(by_key, components)
    return _per_level(levels, (parallel + 2 * perpendicular) / 3 * unit)


def _anisotropy_from_components(
    by_key: Mapping[str, Dataset], components: Mapping[str, str]
) -> AnisotropyTable:
    """Return the transition anisotropy per level from its components.

    ``<gamma> = <alpha_par> - <alpha_perp>``, level by level: the element of the
    anisotropic (rank-2) part of the Q(J) line (optmat.gas.raman_response).

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    components : mapping of str to str
        ``parallel`` and ``perpendicular``, as for :func:`_mean_from_components`.

    Returns
    -------
    AnisotropyTable
        The anisotropy per level, C m^2/V.
    """
    levels, parallel, perpendicular, unit = _components(by_key, components)
    return _per_level(levels, (parallel - perpendicular) * unit)


def _components(
    by_key: Mapping[str, Dataset], components: Mapping[str, str]
) -> tuple[list[int], Float[np.ndarray, " n"], Float[np.ndarray, " n"], float]:
    """Return the parallel and perpendicular transition elements per level.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    components : mapping of str to str
        ``parallel`` and ``perpendicular``: ``raman_transition_polarizability``
        datasets with a row per lower level ``J``.

    Returns
    -------
    tuple of (list of int, numpy.ndarray, numpy.ndarray, float)
        The levels (in the datasets' order), the two elements in the datasets' unit
        and that unit in C m^2/V.

    Raises
    ------
    DataError
        If the two datasets do not cover the same levels.
    """
    parallel, _ = _get(by_key, components["parallel"])
    perpendicular, _ = _get(by_key, components["perpendicular"])
    if parallel.rows != perpendicular.rows:
        msg = f"{parallel.key} and {perpendicular.key} cover different levels"
        raise DataError(msg)
    return (
        [_level(label) for label in parallel.rows],
        parallel.scaled_values(),
        perpendicular.scaled_values(),
        POLARIZABILITY_UNITS[parallel.unit],
    )


def _per_level(levels: list[int], value: Float[np.ndarray, " n"]) -> AnisotropyTable:
    """Return values per level as a table ordered by level.

    Parameters
    ----------
    levels : list of int
        The rotational level of each value.
    value : numpy.ndarray
        The values, C m^2/V.

    Returns
    -------
    AnisotropyTable
        The table.
    """
    order = np.argsort(levels)
    return AnisotropyTable(
        tuple(int(levels[i]) for i in order), tuple(float(value[i]) for i in order)
    )


def _transition_from_ratios(
    by_key: Mapping[str, Dataset],
    ratios: list[str],
    depolarization: str,
    n2_reference: str,
) -> float:
    """Return a transition polarizability from normalised Raman cross sections.

    Schrötter & Klöckner's normalised cross section Sigma_j is the band's
    ``dsigma/dOmega (1 - exp(-hc nu_j/kT)) / (nu_0 - nu_j)^4`` over N2's, whose
    absolute value is ``C = 5.05e-48 cm^6/sr`` (their Table 4.1). For a band with
    transition elements ``a`` (mean) and ``g`` (anisotropy),
    ``dsigma/dOmega = (2 pi nu_s)^4 (a^2 + 7 g^2 / 45)`` summed over the scattered
    polarizations, and the depolarization ratio for linear polarization is
    ``rho = 3 g^2 / (45 a^2 + 4 g^2)``. So
    ``a^2 = Sigma C / (16 pi^4) / (1 + 7 rho / (3 - 4 rho))``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    ratios : list of str
        ``raman_cross_section_ratio`` datasets of Sigma_j; every value is averaged.
    depolarization : str
        A ``depolarization_ratio`` dataset (and row) for linear polarization, or a
        number (as a string) where no ratio is curated.
    n2_reference : str
        The formula of N2's absolute normalised cross section (coefficient ``C``).

    Returns
    -------
    float
        ``<0| alpha_mean |1>``, C m^2/V.
    """
    values = np.concatenate([_get(by_key, key)[0].scaled_values() for key in ratios])
    sigma = float(np.mean(values))
    rho = _depolarization(by_key, depolarization)
    c = _coefficients(n2_reference)["C"]
    mean_squared = sigma * c / (16 * math.pi**4) / (1 + 7 * rho / (3 - 4 * rho))
    return math.sqrt(mean_squared) * POLARIZABILITY_UNITS["cm3"]


def _level(label: str) -> int:
    """Return the rotational level of a row label: ``S(8)``, ``Q(J)``, ``N=7``, ``J=2``.

    Parameters
    ----------
    label : str
        The row label.

    Returns
    -------
    int
        The level.

    Raises
    ------
    DataError
        If the label names no level.
    """
    match = re.search(r"(?:\(|=)\s*(\d+)", label)
    if match is None:
        msg = f"row label {label!r} names no rotational level"
        raise DataError(msg)
    return int(match.group(1))


def _per_amagat(
    dataset: Dataset, temperature: float | None = None
) -> tuple[Float[np.ndarray, " n"], float]:
    """Return a table of broadening coefficients as FWHM per amagat.

    Parameters
    ----------
    dataset : Dataset
        A ``raman_linewidth`` coefficient per pressure (atm or bar, at the dataset's
        ``temperature_k``) or per amagat (MHz, GHz or cm^-1), as a FWHM or a HWHM.
    temperature : float, optional
        The table's temperature, K, where the dataset states none as a number (a range,
        or "room temperature"; ``raman.toml`` says which value is taken and why).

    Returns
    -------
    tuple of (numpy.ndarray, float)
        The coefficients, Hz/amagat (FWHM), and their temperature, K.

    Raises
    ------
    DataError
        If the unit or the width convention is not one of these, or the temperature
        is not known.
    """
    if temperature is None:
        temperature = float(str(dataset.parameters.get("temperature_k", "nan")))
    width = {"fwhm": 1.0, "hwhm": 2.0}.get(str(dataset.parameters.get("width")))
    if width is None or not math.isfinite(temperature):
        msg = f"{dataset.key}: a width table needs its width convention and temperature"
        raise DataError(msg)
    values = dataset.scaled_values() * width
    unit = dataset.unit
    if unit in {"MHz/amagat", "GHz/amagat"}:
        return values * (1e6 if unit == "MHz/amagat" else 1e9), temperature
    wavenumber, _, per = unit.partition("/")
    if wavenumber == "cm-1" and per == "amagat":
        return values * WAVENUMBER_TO_HZ, temperature
    if wavenumber == "cm-1" and per in _PRESSURES_ATM:
        # Per unit pressure at T: one atm at T is ONE_ATM_0C.temperature / T amagat.
        per_atm = values * WAVENUMBER_TO_HZ / _PRESSURES_ATM[per]
        return per_atm * temperature / ONE_ATM_0C.temperature, temperature
    msg = f"{dataset.key}: no conversion of {unit!r} to Hz/amagat"
    raise DataError(msg)


def _line_table(
    by_key: Mapping[str, Dataset], entry: Mapping[str, object]
) -> tuple[LineTable, list[str]]:
    """Build one perturber's line table from its entry in ``raman.toml``.

    The exponents come from a formula's coefficients ``<prefix><J>`` (as a per-amagat
    power law ``(T / T0)^gamma``, so ``N = -gamma``), or from a least-squares fit of
    the log coefficients across the tables at other temperatures, or are zero.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``dataset`` (the table), and optionally ``exponent_formula`` with
        ``exponent_prefix``, or ``temperatures`` (the other tables).

    Returns
    -------
    tuple of (LineTable, list of str)
        The table and the sources used.
    """
    dataset, _ = _get(by_key, str(entry["dataset"]))
    levels = [_level(label) for label in dataset.abscissae]
    stated = entry.get("temperature_k")
    coefficient, temperature = _per_amagat(
        dataset, None if stated is None else float(str(stated))
    )
    order = np.argsort(levels)
    j = np.asarray(levels)[order]
    gamma = coefficient[order]
    sources = [dataset.key]
    exponent = np.zeros(len(j))
    if "exponent_formula" in entry:
        source_id, _, formula_id = str(entry["exponent_formula"]).partition("/")
        (formula,) = [
            f for f in sources_of()[source_id].formulas if f["id"] == formula_id
        ]
        coefficients = _strings(formula["coefficients"], "coefficients")
        prefix = str(entry["exponent_prefix"])
        pairs = sorted(
            (int(name[len(prefix) :]), -float(str(value)))
            for name, value in coefficients.items()
            if name.startswith(prefix) and name[len(prefix) :].isdigit()
        )
        exponent = np.interp(j, [p[0] for p in pairs], [p[1] for p in pairs])
        sources.append(str(entry["exponent_formula"]))
    elif "temperatures" in entry:
        others = [
            _get(by_key, str(key))[0]
            for key in _list(entry["temperatures"], "temperatures")
        ]
        tables = [(dict(zip(levels, coefficient, strict=True)), temperature)]
        for other in others:
            values, t = _per_amagat(other)
            tables.append(
                (dict(zip(map(_level, other.abscissae), values, strict=True)), t)
            )
            sources.append(other.key)
        # Each level of the reference table gets a power law fitted through its
        # coefficients at every temperature (the mean exponent where it has only one).
        # The hot tables reach higher levels: those keep the exponent of the nearest
        # fitted level and take their coefficient at the reference temperature from the
        # geometric mean of their measurements, so a slope is never fitted to the hot
        # points alone (Hölzer et al.'s scatter by 20-40 % between temperatures).
        for i, level in enumerate(j):
            points = [
                (math.log(temperature / t), math.log(v[level]))
                for v, t in tables
                if level in v
            ]
            if len(points) >= 2:
                x, y = np.array(points).T
                exponent[i] = float(np.polyfit(x, y, 1)[0])
        fitted = [
            e
            for e, level in zip(exponent, j, strict=True)
            if sum(level in v for v, _ in tables) >= 2
        ]
        if fitted:
            exponent = np.where(
                [sum(level in v for v, _ in tables) >= 2 for level in j],
                exponent,
                np.mean(fitted),
            )
        extra = sorted({level for v, _ in tables[1:] for level in v} - set(levels))
        if extra:
            exponents = np.interp(extra, j, exponent)
            coefficients = [
                math.exp(
                    np.mean(
                        [
                            math.log(v[level]) - n * math.log(temperature / t)
                            for v, t in tables[1:]
                            if level in v
                        ]
                    )
                )
                for level, n in zip(extra, exponents, strict=True)
            ]
            order = np.argsort(np.concatenate([j, extra]))
            j = np.concatenate([j, extra])[order]
            gamma = np.concatenate([gamma, coefficients])[order]
            exponent = np.concatenate([exponent, exponents])[order]
    return (
        LineTable(
            tuple(int(x) for x in j),
            tuple(map(float, gamma)),
            temperature,
            tuple(map(float, exponent)),
        ),
        sources,
    )


def _coefficients(reference: str) -> dict[str, float]:
    """Return a formula record's coefficients as numbers.

    Parameters
    ----------
    reference : str
        ``source/formula``.

    Returns
    -------
    dict of str to float
        The numeric coefficients (non-numeric ones, such as ranges, are left out).

    Raises
    ------
    DataError
        If the formula does not exist.
    """
    source_id, _, formula_id = reference.partition("/")
    matches = [f for f in sources_of()[source_id].formulas if f["id"] == formula_id]
    if len(matches) != 1:
        msg = f"raman.toml names an unknown formula {reference!r}"
        raise DataError(msg)
    numbers = {}
    for name, value in _strings(matches[0]["coefficients"], "coefficients").items():
        try:
            numbers[str(name)] = float(str(value))
        except ValueError:
            continue
    return numbers


def _law(by_key: Mapping[str, Dataset], entry: Mapping[str, object]) -> LinewidthLaw:
    """Build a single width law from its entry in ``raman.toml``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``law``: ``bischel_dyer`` (their Eq. (2): the Dicke term ``A (T/T0)^n`` and a
        quadratic broadening in T, MHz), ``dicke_linear`` (``A/rho + B rho``, MHz),
        ``constant_linear`` (``c + B rho``, MHz), or ``width`` (a measured FWHM taken
        as density-independent); ``formula`` or ``dataset``.

    Returns
    -------
    LinewidthLaw
        The law.

    Raises
    ------
    DataError
        If the law is unknown.
    """
    kind = str(entry["law"])
    substitutes = _strings(entry.get("substitutes", {}), "substitutes")
    if kind == "dephasing":
        # A dephasing time T2 measured at a stated pressure and temperature: the FWHM
        # 1 / (pi T2), taken as collisional (proportional to the density).
        reference = str(entry["dataset"])
        dataset, _ = _get(by_key, reference)
        seconds = {"ps": 1e-12, "ns": 1e-9, "fs": 1e-15, "s": 1.0}[dataset.unit]
        fwhm = 1 / (math.pi * _scalar(by_key, reference) * seconds)
        temperature = float(str(entry["temperature_k"]))
        pressure = float(str(entry["pressure_atm"]))
        amagat = pressure * ONE_ATM_0C.temperature / temperature
        return LinewidthLaw(
            broadening={"self": Broadening(fwhm / amagat, temperature)},
            substitutes=substitutes,
            sources=(dataset.key,),
        )
    if kind == "width":
        reference = str(entry["dataset"])
        dataset, _ = _get(by_key, reference)
        factor = {"cm-1": WAVENUMBER_TO_HZ, "MHz": 1e6}[dataset.unit]
        return LinewidthLaw(
            constant=_scalar(by_key, reference) * factor, sources=(dataset.key,)
        )
    reference = str(entry["formula"])
    c = _coefficients(reference)
    if kind == "bischel_dyer":
        return LinewidthLaw(
            dicke=c["A298"] * 1e6,
            dicke_temperature=c["T0"],
            dicke_exponent=c["A_exponent"],
            broadening={
                "self": Broadening(
                    c["B0"] * 1e6,
                    c["T0"],
                    linear=c["B1"] * 1e6,
                    quadratic=c["B2"] * 1e6,
                )
            },
            substitutes=substitutes,
            sources=(reference,),
        )
    if kind == "dicke_linear":
        return LinewidthLaw(
            dicke=c["A_MHz_amagat"] * 1e6,
            broadening={"self": Broadening(c["B_MHz_per_amagat"] * 1e6, 298.0)},
            substitutes=substitutes,
            sources=(reference,),
        )
    if kind == "constant_linear":
        return LinewidthLaw(
            constant=c["constant_MHz"] * 1e6,
            broadening={"self": Broadening(c["B_MHz_per_amagat"] * 1e6, 298.0)},
            substitutes=substitutes,
            sources=(reference,),
        )
    msg = f"raman.toml: unknown width law {kind!r}"
    raise DataError(msg)


def _widths(
    by_key: Mapping[str, Dataset],
    entry: Mapping[str, object],
    constants: RotationalConstants | None,
) -> _LinewidthModel:
    """Build a band's width model from its entry in ``raman.toml``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``type``: ``law``, ``table`` (one sub-table per perturber) or
        ``exponential_gap`` (with ``formula``).
    constants : RotationalConstants or None
        The band's lower-level constants (for the exponential-gap law).

    Returns
    -------
    LinewidthLaw or LineTableLaw or ExponentialGapLaw
        The model.

    Raises
    ------
    DataError
        If the type is unknown, or the exponential-gap law has no constants.
    """
    kind = str(entry["type"])
    substitutes = _strings(entry.get("substitutes", {}), "substitutes")
    if kind == "law":
        return _law(by_key, entry)
    if kind == "table":
        tables = {}
        used: list[str] = []
        for perturber, sub in _table(entry["perturbers"], "perturbers").items():
            table, sources = _line_table(by_key, _table(sub, perturber))
            tables[str(perturber)] = table
            used += sources
        dicke = 0.0
        if "dicke" in entry:
            dataset, _ = _get(by_key, str(entry["dicke"]))
            if dataset.unit != "MHz amagat":
                msg = f"{dataset.key}: a Dicke coefficient in MHz amagat is expected"
                raise DataError(msg)
            dicke = _scalar(by_key, str(entry["dicke"])) * 1e6
            used.append(dataset.key)
        elif "dicke_formula" in entry:
            # A formula's Dicke coefficient, in MHz amagat (FWHM).
            reference = str(entry["dicke_formula"])
            dicke = _coefficients(reference)[str(entry["dicke_coefficient"])] * 1e6
            used.append(reference)
        return LineTableLaw(tables, substitutes, dicke, tuple(used))
    if kind == "exponential_gap":
        if constants is None:
            msg = "the exponential-gap law needs the band's rotational constants"
            raise DataError(msg)
        reference = str(entry["formula"])
        c = _coefficients(reference)
        return ExponentialGapLaw(
            alpha=c["alpha"],
            beta=c["beta"],
            delta=c["delta"],
            a=c["a"],
            exponent=c[str(entry.get("exponent", "n_experiment"))],
            reference_temperature=c[
                str(entry.get("reference_temperature", "T0_experiment"))
            ],
            b=constants.b,
            d=constants.d,
            foreign_as_self=tuple(k for k, v in substitutes.items() if v == "self"),
            sources=(reference,),
        )
    msg = f"raman.toml: unknown width type {kind!r}"
    raise DataError(msg)


def _transitions(dataset: Dataset) -> list[tuple[int, int, float]]:
    """Return the lines of a table whose rows are ``J'-J''`` pairs.

    Parameters
    ----------
    dataset : Dataset
        A ``raman_shift`` table with rows such as ``3-1`` or ``1-1``.

    Returns
    -------
    list of tuple of (int, int, float)
        The lower and upper J of each line and its position, cm^-1.
    """
    lines = []
    for row, value in zip(dataset.rows, dataset.scaled_values(), strict=True):
        first, second = (int(x) for x in row.split("-"))
        lines.append((min(first, second), max(first, second), float(value)))
    return lines


def _fit_ground(dataset: Dataset) -> RotationalConstants:
    """Fit B and D to a pure rotational S branch, by least squares.

    Parameters
    ----------
    dataset : Dataset
        The lines, rows ``J'-J''`` (Stokes and anti-Stokes lines both count).

    Returns
    -------
    RotationalConstants
        B and D, cm^-1.
    """
    rows, values = [], []
    for low, high, value in _transitions(dataset):
        x_low, x_high = low * (low + 1.0), high * (high + 1.0)
        rows.append([x_high - x_low, -(x_high**2 - x_low**2)])
        values.append(value)
    (b, d), *_ = np.linalg.lstsq(np.array(rows), np.array(values), rcond=None)
    return RotationalConstants(float(b), float(d))


def _fit_upper(
    dataset: Dataset, lower: RotationalConstants, origin: float
) -> RotationalConstants:
    """Fit the upper level's B and D to a Q branch, given the band origin.

    ``Q(J) = nu_0 + (B_1 - B_0) x - (D_1 - D_0) x^2`` with ``x = J (J + 1)``; only the
    table's Q lines (rows ``J-J``) are used.

    Parameters
    ----------
    dataset : Dataset
        The lines.
    lower : RotationalConstants
        The ground level's constants.
    origin : float
        ``nu_0``, cm^-1.

    Returns
    -------
    RotationalConstants
        B_1 and D_1, cm^-1.
    """
    rows, values = [], []
    for low, high, value in _transitions(dataset):
        if low != high:
            continue
        x = low * (low + 1.0)
        rows.append([x, -(x**2)])
        values.append(value - origin)
    (db, dd), *_ = np.linalg.lstsq(np.array(rows), np.array(values), rcond=None)
    return RotationalConstants(lower.b + float(db), lower.d + float(dd))


def _constants(
    by_key: Mapping[str, Dataset], entry: Mapping[str, object]
) -> RotationalConstants:
    """Build rotational constants from references (cm^-1), or fit them to lines.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``b`` and optionally ``d`` and ``h``; or ``fit_lines``, an S-branch table.

    Returns
    -------
    RotationalConstants
        The constants.
    """
    if "fit_lines" in entry:
        return _fit_ground(_get(by_key, str(entry["fit_lines"]))[0])
    return RotationalConstants(
        b=_wavenumber(by_key, str(entry["b"])),
        d=_wavenumber(by_key, str(entry["d"])) if "d" in entry else 0.0,
        h=_wavenumber(by_key, str(entry["h"])) if "h" in entry else 0.0,
    )


def _anisotropy(
    by_key: Mapping[str, Dataset], reference: str
) -> float | AnisotropyTable:
    """Return an anisotropy, or a table of it per lower level when rows name levels.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        A ``polarizability_anisotropy`` dataset (and row).

    Returns
    -------
    float or AnisotropyTable
        The anisotropy, C m^2/V.
    """
    dataset, index = _get(by_key, reference)
    if index is not None or len(dataset) == 1:
        return _polarizability(by_key, reference)
    levels = [_level(label) for label in dataset.rows]
    values = dataset.scaled_values() * POLARIZABILITY_UNITS[dataset.unit]
    order = np.argsort(levels)
    return AnisotropyTable(
        tuple(int(levels[i]) for i in order), tuple(float(values[i]) for i in order)
    )


def _band_origin(by_key: Mapping[str, Dataset], entry: Mapping[str, object]) -> float:
    """Return a band origin: a ``raman_shift``, or a difference of term values.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``band_origin`` (a reference), or ``upper_term`` and ``lower_term``.

    Returns
    -------
    float
        ``nu_0``, cm^-1.
    """
    if "band_origin" in entry:
        return _wavenumber(by_key, str(entry["band_origin"]))
    return _wavenumber(by_key, str(entry["upper_term"])) - _wavenumber(
        by_key, str(entry["lower_term"])
    )


def _vibrational(
    by_key: Mapping[str, Dataset],
    entry: Mapping[str, object],
    lower: RotationalConstants | None,
) -> VibrationalBand:
    """Build a vibrational band from its entry in ``raman.toml``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    entry : mapping
        ``mode``; ``band_origin`` or ``upper_term`` and ``lower_term``;
        ``transition_polarizability`` or ``derivative``; the upper level's constants
        as ``upper`` (a table of ``b``, ``d``, ``h``), ``fit_lines`` (a Q-branch table
        fitted with the band origin) or ``b_difference`` (B0 - B1, with D the lower
        level's); ``resolved`` (default true when lower constants exist); ``widths``;
        optionally ``depolarization`` (the Q branch's, a dataset and row, or a number
        as a string), which the cross-section route requires.
    lower : RotationalConstants or None
        The ground level's constants.

    Returns
    -------
    VibrationalBand
        The band.

    Raises
    ------
    DataError
        If a resolved band lacks either rotational ladder's prescription.
    """
    origin = _band_origin(by_key, entry)
    sources = {
        "band_origin": str(entry.get("band_origin", entry.get("upper_term", "")))
    }
    transition: float | AnisotropyTable
    anisotropy = None
    if "transition_polarizability" in entry:
        reference = str(entry["transition_polarizability"])
        transition = _polarizability(by_key, reference)
    elif "transition_components" in entry:
        components = _strings(entry["transition_components"], "transition_components")
        transition = _mean_from_components(by_key, components)
        anisotropy = _anisotropy_from_components(by_key, components)
        reference = f"{components['parallel']}, {components['perpendicular']}"
    elif "cross_section_ratios" in entry:
        ratios = [str(k) for k in _list(entry["cross_section_ratios"], "ratios")]
        transition = _transition_from_ratios(
            by_key,
            ratios,
            str(entry["depolarization"]),
            str(entry["n2_reference"]),
        )
        reference = ", ".join(ratios)
    else:
        reference = str(entry["derivative"])
        transition = _transition_from_derivative(by_key, reference, origin)
    sources["transition_polarizability"] = reference
    depolarization = None
    if "depolarization" in entry:
        depolarization = _depolarization(by_key, str(entry["depolarization"]))
        sources["depolarization"] = str(entry["depolarization"])
    resonance = None
    if "dispersion" in entry:
        resonance = _coefficients(str(entry["dispersion"]))["nu_i"]
        sources["dispersion"] = str(entry["dispersion"])
    reference_nm = float(str(entry.get("reference_nm", "inf")))
    upper = None
    if entry.get("resolved", lower is not None):
        # A requested resolved band must not silently become one line (TC-735).
        if lower is None:
            msg = f"{entry['mode']}: a resolved band needs lower rotational constants"
            raise DataError(msg)
        if "upper" in entry:
            upper = _constants(by_key, _table(entry["upper"], "upper"))
        elif "fit_lines" in entry:
            upper = _fit_upper(_get(by_key, str(entry["fit_lines"]))[0], lower, origin)
        elif "b_difference" in entry:
            upper = RotationalConstants(
                lower.b - _wavenumber(by_key, str(entry["b_difference"])),
                lower.d,
                lower.h,
            )
        else:
            msg = (
                f"{entry['mode']}: a resolved band needs "
                "upper, fit_lines or b_difference"
            )
            raise DataError(msg)
    return VibrationalBand(
        mode=str(entry["mode"]),
        band_origin=origin,
        transition_polarizability=transition,
        lower=lower if upper is not None else None,
        upper=upper,
        widths=_widths(by_key, _table(entry["widths"], "widths"), lower),
        sources=sources,
        resonance=resonance,
        reference_wavenumber=0.0 if math.isinf(reference_nm) else 1e7 / reference_nm,
        transition_anisotropy=anisotropy,
        depolarization=depolarization,
        key=str(entry.get("key", entry["mode"])),
    )


def _depolarization(by_key: Mapping[str, Dataset], reference: str) -> float:
    """Return a depolarization ratio: a number written as a string, or a dataset's.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        A ``depolarization_ratio`` dataset (and row), or a number.

    Returns
    -------
    float
        The ratio.
    """
    try:
        return float(reference)
    except ValueError:
        return _scalar(by_key, reference)


def _parameters(
    gas: Gas, entry: Mapping[str, object], by_key: Mapping[str, Dataset]
) -> RamanParameters:
    """Build a gas's Raman parameters from its entry in ``raman.toml``.

    Parameters
    ----------
    gas : Gas
        The gas.
    entry : mapping
        Its table.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    RamanParameters
        The parameters.
    """
    rotational = None
    lower = None
    if "constants" in entry:
        lower = _constants(by_key, _table(entry["constants"], "constants"))
    if "rotational" in entry:
        sub = _table(entry["rotational"], "rotational")
        if lower is None:
            msg = f"{gas.value}: a rotational band needs [<gas>.constants]"
            raise DataError(msg)
        rotational = RotationalBand(
            constants=lower,
            spin=spin_statistics(gas),
            anisotropy=_anisotropy(by_key, str(sub["anisotropy"])),
            widths=_widths(by_key, _table(sub["widths"], "widths"), lower),
            sources={"anisotropy": str(sub["anisotropy"])},
        )
    vibrational = tuple(
        _vibrational(by_key, _table(band, "vibrational"), lower)
        for band in _list(entry.get("vibrational", []), "vibrational")
    )
    reasons = _strings(entry.get("reasons", {}), "reasons")
    return RamanParameters(
        gas, rotational, vibrational, reasons, str(entry.get("missing", ""))
    )


@functools.cache
def _all_parameters() -> dict[Gas, RamanParameters]:
    """Build every gas's Raman parameters from ``raman.toml``.

    Returns
    -------
    dict of Gas to RamanParameters
        By gas.
    """
    with RAMAN_FILE.open("rb") as handle:
        table = tomllib.load(handle)
    by_key = _catalogue.by_key()
    return {
        as_gas(name): _parameters(as_gas(name), entry, by_key)
        for name, entry in table.items()
    }


def raman_parameters(gas: Gas | str) -> RamanParameters:
    """Return the curated Raman parameters of a molecular gas.

    Parameters
    ----------
    gas : Gas or str
        The gas.

    Returns
    -------
    RamanParameters
        Its bands, with the source of every parameter.

    Raises
    ------
    DataError
        If the gas has none (an atom, or a molecule not yet curated).
    """
    species = as_gas(gas)
    parameters = _all_parameters()
    if species not in parameters:
        msg = f"no Raman parameters for {species.value}"
        raise DataError(msg)
    return parameters[species]


@dataclass(frozen=True)
class BandLines:
    """The lines of one Raman band of a gas at a state.

    Attributes
    ----------
    band : RotationalBand or VibrationalBand
        The band, with its curated parameters.
    lines : RamanLines
        Its lines, with their widths, in the convention of this module.
    dispersion_factor : float
        The factor applied to the band's transition elements for the pump wavelength
        asked for (1 for the rotational band, which is not dispersed).
    """

    band: RotationalBand | VibrationalBand
    lines: RamanLines
    dispersion_factor: float = 1.0


def band_lines(
    gas: Gas | str,
    temperature: float,
    density: float,
    fractions: Mapping[str, float] | None = None,
    *,
    wavelength: float | None = None,
    population_floor: float = 1e-12,
) -> tuple[BandLines, ...]:
    """Return a gas's Raman lines band by band, with their widths, at a state.

    The parameters are those of :func:`raman_lines`, which joins these lines.

    Parameters
    ----------
    gas : Gas or str
        The gas.
    temperature : float
        Temperature, K, > 0.
    density : float
        Total number density, m^-3, > 0.
    fractions : mapping of str to float, optional
        Mole fractions of the collision partners; the pure gas by default.
    wavelength : float, optional
        Vacuum wavelength of the optical field, m, for the vibrational strengths'
        dispersion.
    population_floor : float, optional
        Lines from levels less populated than this fraction are left out.

    Returns
    -------
    tuple of BandLines
        The rotational band (if any), then each vibrational band.
    """
    p = raman_parameters(gas)
    parts: list[BandLines] = []
    if p.rotational is not None:
        r = p.rotational
        lines = rotational_lines(
            r.constants,
            r.spin,
            r.anisotropy,
            temperature,
            population_floor=population_floor,
        )
        parts.append(
            BandLines(
                r,
                dataclasses.replace(
                    lines, fwhm=r.widths.fwhm(lines.j, temperature, density, fractions)
                ),
            )
        )
    for band in p.vibrational:
        spin = spin_statistics(p.gas) if band.lower is not None else None
        factor = band.dispersion_factor(wavelength)
        element = band.transition_polarizability
        element = (
            AnisotropyTable(element.j, tuple(v * factor for v in element.value))
            if isinstance(element, AnisotropyTable)
            else element * factor
        )
        lines = vibrational_lines(
            band.band_origin,
            element,
            lower=band.lower,
            upper=band.upper,
            spin=spin,
            temperature=temperature,
            population_floor=population_floor,
        )
        parts.append(
            BandLines(
                band,
                dataclasses.replace(
                    lines,
                    fwhm=band.widths.fwhm(lines.j, temperature, density, fractions),
                ),
                factor,
            )
        )
    return tuple(parts)


def raman_lines(
    gas: Gas | str,
    temperature: float,
    density: float,
    fractions: Mapping[str, float] | None = None,
    *,
    wavelength: float | None = None,
    population_floor: float = 1e-12,
) -> RamanLines:
    """Return a gas's Raman lines, with their widths, at a state.

    Parameters
    ----------
    gas : Gas or str
        The gas.
    temperature : float
        Temperature, K, > 0.
    density : float
        Total number density, m^-3, > 0 (the widths depend on it; the amplitudes are
        per molecule).
    fractions : mapping of str to float, optional
        Mole fractions of the collision partners (``self`` for the gas itself); the
        pure gas by default.
    wavelength : float, optional
        Vacuum wavelength of the optical field, m: the vibrational strengths are
        scaled by their dispersion law where one is curated (TC-542). By default the
        strengths are as curated (static, or at the stated reference wavelength).
    population_floor : float, optional
        Lines from levels less populated than this fraction are left out.

    Returns
    -------
    RamanLines
        The rotational S-branch lines, then each vibrational band's Q-branch lines, in
        the convention of this module.
    """
    parts = [
        part.lines
        for part in band_lines(
            gas,
            temperature,
            density,
            fractions,
            wavelength=wavelength,
            population_floor=population_floor,
        )
    ]
    return RamanLines(
        branch=tuple(b for part in parts for b in part.branch),
        j=np.concatenate([part.j for part in parts]),
        frequency=np.concatenate([part.frequency for part in parts]),
        amplitude=np.concatenate([part.amplitude for part in parts]),
        fwhm=np.concatenate([np.asarray(part.fwhm) for part in parts]),
    )


def stokes_gain(
    lines: RamanLines,
    pump_wavelength: float,
    density: float,
) -> Float[np.ndarray, " n"]:
    """Return the steady-state Stokes gain of each line, for a pump at its centre.

    With a real pump field ``E_p`` and Stokes field ``E_s`` polarised alike, the
    response ``H(Delta) = int h(tau) exp(i Delta tau) d tau`` of one line at the Raman
    resonance ``Delta = Omega`` is ``i A / (2 Gamma)``. The Stokes polarisation
    ``(1/2) eps0 N |E_p|^2 E_s H*`` then gives the intensity gain
    ``g = omega_s N A / (2 c^2 eps0 Gamma)``, with ``Gamma = pi FWHM`` and
    ``omega_s = omega_p - Omega``. That is the plane-wave, narrow-band (steady-state)
    gain, ``dI_s / dz = g I_p I_s``.

    Parameters
    ----------
    lines : RamanLines
        Lines with their widths (from :func:`raman_lines`).
    pump_wavelength : float
        Vacuum wavelength of the pump, m.
    density : float
        Number density of the Raman-active gas, m^-3.

    Returns
    -------
    numpy.ndarray
        The gain of each line, m/W.

    Raises
    ------
    DataError
        If the lines carry no widths.
    """
    if lines.fwhm is None:
        msg = "the gain needs the lines' widths: use raman_lines"
        raise DataError(msg)
    omega_s = 2 * math.pi * C / pump_wavelength - lines.frequency
    rate = math.pi * lines.fwhm
    return omega_s * density * lines.amplitude / (2 * C**2 * EPS0 * rate)


def rotational_nonlinear_index(lines: RamanLines, density: float) -> float:
    """Return the rotational Raman nonlinear index for a pulse long against the lines.

    A pulse long against every rotational period sees the time integral of the
    response, ``R0 = sum_J A_J Omega_J / (Omega_J^2 + Gamma_J^2)``. In the convention
    of this module the index change is then ``N R0 I / (2 eps0 c)`` (``n0 = 1``). The
    response is slow, so a probe sees the same shift as the pulse itself (only the
    slowly varying part of ``E^2`` drives it), unlike the instantaneous electronic
    response, whose probe coefficient is twice the self one.

    The integral is that of the coherent (Delta J = +-2) response. The isothermal,
    classical orientational response ``Delta alpha^2 / (22.5 eps0 k_B T)`` is 4/3
    larger at high temperature: it also contains the redistribution of the M
    sublevels within a level, which collisions establish but a coherent drive does not
    (TC-493).

    Parameters
    ----------
    lines : RamanLines
        Lines with their widths (from :func:`raman_lines`); only the S branch is
        used.
    density : float
        Number density of the gas, m^-3.

    Returns
    -------
    float
        ``n2``, m^2/W.

    Raises
    ------
    DataError
        If the lines carry no widths.
    """
    if lines.fwhm is None:
        msg = "the index needs the lines' widths: use raman_lines"
        raise DataError(msg)
    s = np.array(lines.branch) == "S"
    omega = lines.frequency[s]
    rate = math.pi * lines.fwhm[s]
    integral = float(np.sum(lines.amplitude[s] * omega / (omega**2 + rate**2)))
    return density * integral / (2 * EPS0 * C)

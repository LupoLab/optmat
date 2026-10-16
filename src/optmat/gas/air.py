# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Air as a propagation model: the sum of its constituents (OEP 1, section 4.4).

For a named composition of mole fractions ``x_s`` the air polarizability is the sum of
its constituents' recommended models,

    alpha_air(omega) = sum_s x_s alpha_s(omega),

so that it is consistent with the populations a propagation code carries, has a
frequency derivative, and needs no second source. A refractometric formula of air
(Ciddor's, Edlen's, Borzsonyi et al.'s) is an independent source, and serves as its
validator.

The named compositions come from CIPM-2007 (:mod:`optmat.gas.composition`):
``"dry air (N2, O2, Ar)"`` (the default, and the air of the Raman response),
``"dry air (N2, O2)"`` (the argon-free reduction) and ``"dry air (CIPM-2007)"``, all
twelve constituents, renormalised by their printed sum, 1 - 2.9e-6.

The number density of air comes from CoolProp's pseudo-pure dry air (Lemmon et al.
2000), which holds the unlike-pair interactions: the constituents' partial-pressure
densities are 1.9e-4 lower at 1 bar and 294 K. The pair terms of the index itself are
OEP 4's; here Lorentz-Lorenz takes none.
"""

import functools
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat import _validity
from optmat.conditions import ONE_ATM_0C, State
from optmat.constants import AU_ANGULAR_FREQUENCY, C
from optmat.exceptions import ModelRangeWarning, OptmatError
from optmat.gas.composition import REDUCED, Constituent, composition
from optmat.gas.eos import (
    EquationOfState,
    air_equation_of_state,
    ideal_number_density,
)
from optmat.gas.infrared import gap_cells, gap_digest
from optmat.gas.lorentz_lorenz import refractive_index, refractivity
from optmat.gas.models.base import Model
from optmat.gas.models.fitted import FittedModel
from optmat.gas.models.registry import default_model
from optmat.identity import Identified

__all__ = [
    "ABSORPTION_THRESHOLD_M",
    "AIR_COMPOSITIONS",
    "ALPHA_FRACTION",
    "DEFAULT_AIR",
    "GVD_FRACTION",
    "MIN_WIDTH_CM",
    "AirModel",
    "air_model",
]

DEFAULT_AIR: Final = "dry air (N2, O2, Ar)"
"""The default air: CIPM-2007's N2, O2 and Ar, renormalised."""

AIR_COMPOSITIONS: Final = (*REDUCED, "dry air (CIPM-2007)")
"""The named compositions an air model is built for; ``"dry air"`` is the default."""


def _intersection(sets: list[tuple[tuple[float, float], ...]]) -> tuple:
    """Return the intersection of several sets of disjoint intervals.

    Parameters
    ----------
    sets : list of tuple of tuple of float
        Each a model's qualified intervals, m.

    Returns
    -------
    tuple of tuple of float
        The intervals common to every set.

    Raises
    ------
    OptmatError
        If they have no wavelength in common.
    """
    common = sets[0]
    for other in sets[1:]:
        common = tuple(
            (max(a, c), min(b, d))
            for a, b in common
            for c, d in other
            if max(a, c) < min(b, d)
        )
    if not common:
        msg = "the constituents' models share no wavelength"
        raise OptmatError(msg)
    return common


ABSORPTION_THRESHOLD_M: Final = 0.1
"""Transparency: the mixture's smoothed absorption below this, m^-1 (OEP 3, sec. 4)."""

ALPHA_FRACTION: Final = 0.1
"""The mixture's tolerance in alpha: this fraction of its statistical one-sigma."""

GVD_FRACTION: Final = 0.01
"""Its tolerance in ``D_2``: this fraction of the larger of ``|D_2|`` and its electronic
terms' ``|D_2|``."""

MIN_WIDTH_CM: Final = 50.0
"""A widened interval that adds a piece narrower than this, cm^-1, is not kept."""


def _tolerances(
    fractions: Float[np.ndarray, " s"], models: tuple[Model, ...], omega: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return a mixture's tolerances in alpha and D_2 (OEP 3, Q13).

    Parameters
    ----------
    fractions : numpy.ndarray
        The constituents' mole fractions.
    models : tuple of Model
        Their models (fits, with a covariance).
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    tuple of numpy.ndarray
        ``ALPHA_FRACTION`` of the mixture's statistical one-sigma of alpha, the
        constituents' independent sigmas combined; ``GVD_FRACTION`` of the larger of
        the mixture's |D_2| and its electronic terms' |D_2| (atomic units).
    """
    variance = np.zeros_like(omega)
    d2 = np.zeros_like(omega)
    electronic = np.zeros_like(omega)
    for x, model in zip(fractions, models, strict=True):
        assert isinstance(model, FittedModel)
        variance += (x * model.polarizability_sigma_au(omega)) ** 2
        d2 += x * model.dispersion_au(omega, 2)
        n = len(model.strengths) - model.n_fixed
        f = np.asarray(model.strengths[:n])[:, None]
        e = np.asarray(model.energies[:n])[:, None]
        # D_2 = (f / 2) 2! [(E - w)^-3 - (E + w)^-3] per pole.
        electronic += x * np.sum(f * ((e - omega) ** -3 - (e + omega) ** -3), axis=0)
    tol2 = GVD_FRACTION * np.maximum(np.abs(d2), np.abs(electronic))
    return ALPHA_FRACTION * np.sqrt(variance), tol2


def _trace_errors(
    fractions: Float[np.ndarray, " s"],
    models: tuple[Model, ...],
    tables: list[dict[str, np.ndarray] | None],
    middle: Float[np.ndarray, " n"],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the unqualified constituents' sums at the middles of the segments.

    Parameters
    ----------
    fractions : numpy.ndarray
        The constituents' mole fractions.
    models : tuple of Model
        Their models.
    tables : list of dict or None
        Their gap tables.
    middle : numpy.ndarray
        The segments' middles, cm^-1; each lies inside one cell of every table.

    Returns
    -------
    tuple of numpy.ndarray
        Where the segment can be qualified at all (every unqualified constituent has
        a cell there, without a pole); the mole-fraction sums of their absorption
        (m^-1) and of their errors in alpha and D_2 (atomic units).
    """
    wavelength = 1e-2 / middle
    absorption = np.zeros(middle.size)
    error0 = np.zeros(middle.size)
    error2 = np.zeros(middle.size)
    possible = np.ones(middle.size, dtype=bool)
    for x, model, table in zip(fractions, models, tables, strict=True):
        unqualified = _validity.outside(wavelength, model.intervals)
        if table is None:
            possible &= ~unqualified
            continue
        index = np.searchsorted(table["high_cm"], middle).clip(
            0, table["high_cm"].size - 1
        )
        covered = (middle >= table["low_cm"][index]) & (
            middle <= table["high_cm"][index]
        )
        possible &= ~unqualified | (covered & ~table["pole"][index])
        here = unqualified & covered
        absorption += np.where(here, x * table["absorption_m"][index], 0.0)
        error0 += np.where(here, x * table["b0"][index], 0.0)
        error2 += np.where(here, x * table["b2"][index], 0.0)
    return possible, absorption, error0, error2


def _widened(
    constituents: tuple[Constituent, ...],
    models: tuple[Model, ...],
    base: tuple[tuple[float, float], ...],
) -> tuple[tuple[float, float], ...]:
    """Return a mixture's intervals, widened where a trace constituent allows (Q13).

    The intersection of the constituents' intervals widens across a constituent's gap
    only where three checks pass for the mixture at its mole fractions (OEP 3, Q13):
    its smoothed absorption, the mole-fraction sum of each unqualified constituent's,
    is below ``ABSORPTION_THRESHOLD_M``; no unqualified constituent has a pole there;
    and the mole-fraction sums of their representation errors are within the
    mixture's tolerances (:func:`_tolerances`). The errors and absorption come from
    each molecule's gap table (:func:`optmat.gas.infrared.gap_cells`); where an
    unqualified constituent has none, or beyond its table, the mixture is not
    widened. Widening changes where the model warns, not its values.

    Parameters
    ----------
    constituents : tuple of Constituent
        The mixture.
    models : tuple of Model
        Their models.
    base : tuple of tuple of float
        The intersection of their intervals, m.

    Returns
    -------
    tuple of tuple of float
        The widened intervals, m, ascending.
    """
    fractions = np.array([c.mole_fraction for c in constituents])
    tables = [gap_cells(model.gas) for model in models]
    edges: set[float] = set()
    for table in tables:
        if table is not None:
            edges.update(table["low_cm"].tolist())
            edges.update(table["high_cm"].tolist())
    if not edges:
        return base
    span = (min(edges), max(edges))
    for model in models:
        for low, high in model.intervals:
            for wavelength in (low, high):
                wavenumber = 1e-2 / wavelength
                if span[0] < wavenumber < span[1]:
                    edges.add(wavenumber)
    grid = np.array(sorted(edges))
    middle = 0.5 * (grid[:-1] + grid[1:])
    wavelength = 1e-2 / middle
    possible, absorption, error0, error2 = _trace_errors(
        fractions, models, tables, middle
    )
    possible &= absorption < ABSORPTION_THRESHOLD_M
    test = possible & ~np.all(
        [~_validity.outside(wavelength, m.intervals) for m in models], axis=0
    )
    if np.any(test):
        hartree = 2 * np.pi * C / AU_ANGULAR_FREQUENCY
        ends = np.concatenate([grid[:-1][test], grid[1:][test]])
        tol0, tol2 = _tolerances(fractions, models, 1e2 * ends * hartree)
        half = int(np.count_nonzero(test))
        tol0 = np.minimum(tol0[:half], tol0[half:])
        tol2 = np.minimum(tol2[:half], tol2[half:])
        passed = (error0[test] <= tol0) & (error2[test] <= tol2)
        possible[np.flatnonzero(test)[~passed]] = False
    pieces = [
        (1e-2 / high, 1e-2 / low)
        for low, high in zip(grid[:-1][possible], grid[1:][possible], strict=True)
    ]
    return _merge(pieces, base)


def _merge(
    pieces: list[tuple[float, float]], base: tuple[tuple[float, float], ...]
) -> tuple[tuple[float, float], ...]:
    """Return qualified pieces merged with the base intervals.

    Parameters
    ----------
    pieces : list of tuple of float
        Qualified wavelength pieces, m.
    base : tuple of tuple of float
        The intersection of the constituents' intervals, m.

    Returns
    -------
    tuple of tuple of float
        The union, its touching pieces joined; a new interval narrower than
        ``MIN_WIDTH_CM`` that holds no base interval is dropped.
    """
    merged: list[list[float]] = []
    for low, high in sorted([*pieces, *base]):
        if merged and low <= merged[-1][1] * (1 + 1e-12):
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    kept = [
        (low, high)
        for low, high in merged
        if any(a >= low and b <= high for a, b in base)
        or 1e-2 / low - 1e-2 / high >= MIN_WIDTH_CM
    ]
    return tuple((float(a), float(b)) for a, b in kept)


@dataclass(frozen=True)
class AirModel(Identified):
    """A named composition of air as a dispersion model.

    It has the interface of a gas dispersion model
    (:class:`~optmat.gas.models.base.Model`): the complex polarizability per molecule
    of air, its frequency derivative, the index at a state, the qualified intervals
    and an identity. It is bound to its constituents' models, whose digests enter its
    own.

    Attributes
    ----------
    name : str
        The composition's name.
    constituents : tuple of Constituent
        The constituents and their mole fractions, summing to one.
    models : tuple of Model
        Each constituent's recommended dispersion model, in the same order.
    equation_of_state : EquationOfState or None
        The pseudo-pure dry air's (Lemmon et al. 2000); None where CoolProp (the
        ``gas`` extra) is not installed, when only the ideal-gas density is
        available (TC-1092).
    wavelength_range : tuple of float
        The envelope of :attr:`intervals`, m.
    intervals : tuple of tuple of float
        The intersection of the constituents' qualified intervals, widened across a
        trace constituent's gaps where the mixture passes OEP 3's checks (Q13), m.
    status : str
        ``validated``.
    """

    IDENTITY_KIND: ClassVar[str] = "gas_refractive"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.lorentz_lorenz",)

    name: str
    constituents: tuple[Constituent, ...]
    models: tuple[Model, ...]
    equation_of_state: EquationOfState | None
    wavelength_range: tuple[float, float]
    intervals: tuple[tuple[float, float], ...]

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, composition and source records.

        Returns
        -------
        tuple
            ``("air/<composition>", composition, sources)``.
        """
        sources = sorted({c.source.split("/", 1)[0] for c in self.constituents})
        return f"air/{self.name}", self.name, tuple(sources)

    def _check_range(self, wavelength: Float[np.ndarray, "*shape"]) -> None:
        """Warn when any wavelength lies outside every one of :attr:`intervals`.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.
        """
        if np.any(_validity.outside(wavelength, self.intervals)):
            msg = (
                f"air model {self.name} evaluated outside "
                f"{_validity.describe(self.intervals)}, where its constituents' models "
                "are qualified"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)

    def _sum(
        self, method: str, wavelength: ArrayLike | float, warn: bool
    ) -> np.ndarray:
        """Return the mole-fraction sum of a method of the constituents' models.

        Parameters
        ----------
        method : str
            ``complex_polarizability`` or ``dpolarizability_domega``.
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool
            Warn outside the intervals.

        Returns
        -------
        numpy.ndarray
            The sum, complex.
        """
        wavelength = _validity.checked_wavelength(
            wavelength, f"air/{self.name}", static=True
        )
        if warn:
            self._check_range(wavelength)
        return np.asarray(
            sum(
                c.mole_fraction * getattr(m, method)(wavelength, warn=False)
                for c, m in zip(self.constituents, self.models, strict=True)
            ),
            dtype=complex,
        )

    def complex_polarizability(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex polarizability per molecule of air.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; ``inf`` is the static limit.
        warn : bool, optional
            Warn (:class:`~optmat.exceptions.ModelRangeWarning`) outside every one of
            :attr:`intervals`; True by default.

        Returns
        -------
        numpy.ndarray
            ``sum_s x_s alpha_s``, C m^2/V (complex).

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.

        Examples
        --------
        >>> from optmat.constants import AU_POLARIZABILITY
        >>> air = air_model()
        >>> f"{float(air.complex_polarizability(800e-9).real / AU_POLARIZABILITY):.3f}"
        '11.591'
        """
        return self._sum("complex_polarizability", wavelength, warn)

    def polarizability(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Float[np.ndarray, "*shape"]:
        """Return the real polarizability per molecule of air.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside the intervals; True by default.

        Returns
        -------
        numpy.ndarray
            C m^2/V.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        """
        return np.asarray(np.real(self.complex_polarizability(wavelength, warn=warn)))

    def dpolarizability_domega(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the frequency derivative of the polarizability per molecule.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside the intervals; True by default.

        Returns
        -------
        numpy.ndarray
            ``sum_s x_s d alpha_s / d omega``, C m^2 V^-1 s (complex).

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        """
        return self._sum("dpolarizability_domega", wavelength, warn)

    def _density(self, state: State, ideal: bool) -> float:
        """Return the number density of air at a state.

        Parameters
        ----------
        state : State
            Pressure and temperature.
        ideal : bool
            Use the ideal gas instead of the pseudo-pure dry air's equation.

        Returns
        -------
        float
            The number density, m^-3.

        Raises
        ------
        ImportError
            For the real-gas density where CoolProp is not installed.
        """
        if ideal:
            return ideal_number_density(state)
        if self.equation_of_state is None:
            msg = "real-gas densities need CoolProp: install optmat[gas]"
            raise ImportError(msg)
        return self.equation_of_state.number_density(state)

    def refractivity(
        self,
        wavelength: ArrayLike | float,
        state: State = ONE_ATM_0C,
        *,
        ideal: bool = False,
        warn: bool = True,
    ) -> Float[np.ndarray, "*shape"]:
        """Return ``n - 1`` of air at a state, through Lorentz-Lorenz.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State, optional
            Pressure and temperature; 0 degrees C and 1 atm by default.
        ideal : bool, optional
            Use the ideal-gas density instead of the pseudo-pure dry air's.
        warn : bool, optional
            Warn outside the intervals; True by default.

        Returns
        -------
        numpy.ndarray
            The refractivity ``n - 1``, without pair terms (OEP 4's).

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        ImportError
            If ``ideal`` is False and CoolProp is not installed.
        """
        alpha = self.polarizability(wavelength, warn=warn)
        return refractivity(alpha, self._density(state, ideal))

    def refractive_index(
        self,
        wavelength: ArrayLike | float,
        state: State = ONE_ATM_0C,
        *,
        ideal: bool = False,
        warn: bool = True,
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex refractive index of air at a state.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State, optional
            Pressure and temperature; 0 degrees C and 1 atm by default.
        ideal : bool, optional
            Use the ideal-gas density instead of the pseudo-pure dry air's.
        warn : bool, optional
            Warn outside the intervals; True by default.

        Returns
        -------
        numpy.ndarray
            ``n``, complex.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        ImportError
            If ``ideal`` is False and CoolProp is not installed.
        """
        alpha = self.complex_polarizability(wavelength, warn=warn)
        return refractive_index(alpha, self._density(state, ideal))

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the digests of the gap tables the intervals were widened with.

        The models are bound by their own digests, through the fields.

        Returns
        -------
        Mapping
            ``gap_tables``: the SHA-256 of each constituent's gap table (OEP 3, Q13),
            by gas; empty for a composition without one.
        """
        tables = {
            str(m.gas): gap_digest(m.gas)
            for m in self.models
            if gap_cells(m.gas) is not None
        }
        return {"gap_tables": tables} if tables else {}


def air_model(name: str = DEFAULT_AIR) -> AirModel:
    """Return the air model of a named composition.

    Parameters
    ----------
    name : str, optional
        One of :data:`AIR_COMPOSITIONS`; ``"dry air"`` is :data:`DEFAULT_AIR`, as for
        the Raman response.

    Returns
    -------
    AirModel
        The composition's model, one per name.

    Raises
    ------
    OptmatError
        For another name.

    Examples
    --------
    >>> air_model().identity.name
    'air/dry air (N2, O2, Ar)'
    >>> len(air_model("dry air (CIPM-2007)").constituents)
    12
    """
    name = DEFAULT_AIR if name == "dry air" else name
    if name not in AIR_COMPOSITIONS:
        msg = f"no air model named {name!r}; the air models are {AIR_COMPOSITIONS}"
        raise OptmatError(msg)
    return _air_model(name)


def _equation_of_state() -> EquationOfState | None:
    """Return the pseudo-pure dry air's equation of state, or None without CoolProp.

    The equation is needed only for the real-gas density; without the ``gas`` extra
    the air model still gives the polarizability and the ideal-gas index (TC-1092).

    Returns
    -------
    EquationOfState or None
        :func:`~optmat.gas.eos.air_equation_of_state`, or None where CoolProp is
        not installed.
    """
    try:
        return air_equation_of_state()
    except ImportError:
        return None


@functools.cache
def _air_model(name: str) -> AirModel:
    """Build the air model of a named composition, once per name.

    Parameters
    ----------
    name : str
        One of :data:`AIR_COMPOSITIONS`.

    Returns
    -------
    AirModel
        The model.
    """
    constituents = composition(name)
    total = sum(c.mole_fraction for c in constituents)
    # The full composition is renormalised by its printed sum (1 - 2.9e-6); the reduced
    # ones sum to one already.
    constituents = tuple(
        Constituent(
            material=c.material,
            mole_fraction=c.mole_fraction / total,
            sigma=c.sigma / total,
            source=c.source,
            primary=c.primary,
        )
        for c in constituents
    )
    models = tuple(default_model(c.material) for c in constituents)
    intervals = _widened(
        constituents, models, _intersection([m.intervals for m in models])
    )
    full = composition("dry air (CIPM-2007)")
    kept = sum(
        c.mole_fraction
        for c in full
        if c.material in {m.material for m in constituents}
    )
    omitted = (
        ""
        if name == "dry air (CIPM-2007)"
        else f"; it leaves out {sum(c.mole_fraction for c in full) - kept:.2e} of "
        "CIPM-2007's molecules"
    )
    return AirModel(
        name=name,
        constituents=constituents,
        models=models,
        equation_of_state=_equation_of_state(),
        wavelength_range=(intervals[0][0], intervals[-1][1]),
        intervals=intervals,
        status="validated",
        status_reason="the sum of its constituents' recommended models" + omitted,
    )

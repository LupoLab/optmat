# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The model interface: polarizability against vacuum wavelength, with a valid range."""

import abc
import functools
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, ClassVar

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat import _validity
from optmat.conditions import ONE_ATM_0C, State
from optmat.exceptions import ModelAbsorptionWarning, ModelRangeWarning
from optmat.gas.eos import ideal_number_density, number_density
from optmat.gas.gases import Gas, info
from optmat.gas.lorentz_lorenz import refractive_index, refractivity
from optmat.gas.virial import virial_ratio
from optmat.identity import Identified

__all__ = ["Model"]


@dataclass(frozen=True)
class Model(Identified):
    """A dispersion model of one gas.

    Attributes
    ----------
    gas : Gas or str
        The species: a curated gas, or a custom species' name (a model a caller
        built), which has no curated virial ratio or equation of state.
    name : str
        A unique name, ``<source>/<formula>`` for a published formula.
    citation : str
        Where the model comes from.
    wavelength_range : tuple of float
        The vacuum-wavelength range over which the source supports the model, m
        (the upper end may be ``inf``): the envelope of :attr:`intervals`.
    kind : str
        The model's provenance: "published" (a formula as printed, corrected where it
        is misprinted), "legacy" (a recorded implementation formula without a
        publication of its own), "theory" (an ab initio or semi-empirical calculation)
        or "fitted" (fitted here). It is not ``identity.kind``, which names the
        quantity (``"gas_refractive"``).
    intervals : tuple of tuple of float
        The qualified intervals, vacuum wavelength, m: disjoint and ascending, with
        :attr:`wavelength_range` as their envelope. One interval, the range, unless
        the source qualifies several.
    status : str
        ``"validated"`` for the recommended model of its gas, ``"comparison"`` for the
        other curated models, ``"user"`` for a model a caller built.
    status_reason : str
        Why.
    """

    IDENTITY_KIND: ClassVar[str] = "gas_refractive"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.lorentz_lorenz",)

    absorption_warning_threshold: ClassVar[float] = 1e-4
    """``|Im alpha / Re alpha|`` above which the real-part methods warn."""
    INPUTS: ClassVar[tuple[str, ...]] = ("intervals",)
    DERIVED: ClassVar[tuple[str, ...]] = ("wavelength_range",)

    gas: Gas | str
    name: str
    citation: str
    wavelength_range: tuple[float, float]
    kind: str
    intervals: tuple[tuple[float, float], ...] = field(
        default=(), kw_only=True, repr=False
    )

    def __post_init__(self) -> None:
        """Take ownership of the fields, and check the intervals.

        Raises
        ------
        ValueError
            If the intervals are not disjoint and ascending with the range as their
            envelope.
        """
        super().__post_init__()
        intervals = _validity.qualified_intervals(
            self.intervals, self.wavelength_range, self.name
        )
        object.__setattr__(self, "intervals", intervals)

    def _override[T: "Model"](self: T, changes: Mapping[str, Any]) -> T:
        """Return the model with some inputs changed; the range follows the intervals.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values.

        Returns
        -------
        Model
            The rebuilt model.
        """
        changed = dict(changes)
        if "intervals" in changed:
            intervals = tuple(tuple(i) for i in changed["intervals"])
            changed["intervals"] = intervals
            changed["wavelength_range"] = (intervals[0][0], intervals[-1][1])
        return replace(self, **changed)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and source records.

        Returns
        -------
        tuple
            ``(name, gas, sources)``.
        """
        return self.name, str(self.gas), self._sources()

    def _sources(self) -> tuple[str, ...]:
        """Return the source records the model rests on.

        Returns
        -------
        tuple of str
            None known to the base class.
        """
        return ()

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the refractivity virial ratio, which the index methods read.

        Returns
        -------
        Mapping
            ``{"refractivity_virial": B_R / A_R}``, m^3/mol.
        """
        return {"refractivity_virial": self.refractivity_virial}

    @functools.cached_property
    def refractivity_virial(self) -> float:
        """The refractivity virial ratio ``B_R / A_R`` of the gas, resolved once.

        Returns
        -------
        float
            The ratio, m^3/mol; zero where none is curated.
        """
        return virial_ratio(self.gas) if isinstance(self.gas, Gas) else 0.0

    def _density(self, state: State, *, ideal: bool) -> float:
        """Return the number density at a state.

        A species without an equation of state (ozone, a custom species) is taken
        as an ideal gas.

        Parameters
        ----------
        state : State
            Pressure and temperature.
        ideal : bool
            Use the ideal-gas density.

        Returns
        -------
        float
            The number density, m^-3.
        """
        if not isinstance(self.gas, Gas) or info(self.gas).fluid is None:
            return ideal_number_density(state)
        return number_density(self.gas, state, ideal=ideal)

    @property
    def absorbing(self) -> bool:
        """Whether the model has an imaginary part (the gas absorbs in its range).

        Returns
        -------
        bool
            False for a transparent model; an absorbing model's real-part methods warn
            where the absorption is not negligible.
        """
        return False

    @abc.abstractmethod
    def _polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Return the polarizability at vacuum wavelengths, without a range check.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m (``inf`` for the static limit).

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V.
        """

    @abc.abstractmethod
    def _dcomplex_polarizability_domega(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d alpha / d omega, analytically, without a range check.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m (``inf`` for the static limit).

        Returns
        -------
        numpy.ndarray
            C m^2 V^-1 s, complex.
        """

    def _complex_polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex polarizability, without a range check.

        Transparent models return the real polarizability; an absorbing model
        overrides this with its imaginary part (``Im alpha > 0`` is absorption).

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m (``inf`` for the static limit).

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V (complex).
        """
        return np.asarray(self._polarizability(wavelength), dtype=complex)

    def _check_range(self, wavelength: Float[np.ndarray, "*shape"]) -> None:
        """Warn when any wavelength lies outside every one of :attr:`intervals`.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.
        """
        if np.any(_validity.outside(wavelength, self.intervals)):
            msg = (
                f"{self.name} ({self.gas}) evaluated outside "
                f"{_validity.describe(self.intervals)}, "
                f"the range its source supports{self._gap_note(wavelength)}"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)

    def _gap_note(self, wavelength: Float[np.ndarray, "*shape"]) -> str:  # noqa: ARG002 - a subclass reads it
        """Return what the message adds for wavelengths between two intervals.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        str
            Empty for a model without recorded gap reasons.
        """
        return ""

    def complex_polarizability(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex per-particle polarizability at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; ``inf`` is the static limit.
        warn : bool, optional
            Warn (:class:`~optmat.exceptions.ModelRangeWarning`) when a
            wavelength lies outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V (complex, ``Im alpha > 0`` is absorption; zero
            imaginary part for a transparent model), the shape of ``wavelength``.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.

        Warns
        -----
        ModelRangeWarning
            If ``warn`` and any wavelength lies outside the model's range.
        """
        wavelength = _validity.checked_wavelength(wavelength, self.name, static=True)
        if warn:
            self._check_range(wavelength)
        return np.asarray(self._complex_polarizability(wavelength))

    def dpolarizability_domega(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the frequency derivative of the complex polarizability.

        ``d alpha / d omega`` at fixed state, with ``omega = 2 pi c / lambda``: it
        enters the group index and every frame derivative of a propagation code. Each
        class gives its own analytic form (OEP 1, section 4.1), computed in ``omega``
        so that the static limit, where it is 0, needs no chain rule through
        ``lambda``.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; ``inf`` is the static limit.
        warn : bool, optional
            Warn (:class:`~optmat.exceptions.ModelRangeWarning`) when a wavelength
            lies outside every one of :attr:`intervals`; True by default.

        Returns
        -------
        numpy.ndarray
            ``d alpha / d omega``, C m^2 V^-1 s (complex), the shape of
            ``wavelength``.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.

        Warns
        -----
        ModelRangeWarning
            If ``warn`` and any wavelength lies outside the model's intervals.

        Examples
        --------
        >>> from optmat.gas.models.registry import model
        >>> float(model("fit/Ar").dpolarizability_domega(float("inf")).real)
        0.0
        """
        wavelength = _validity.checked_wavelength(wavelength, self.name, static=True)
        if warn:
            self._check_range(wavelength)
        return np.asarray(self._dcomplex_polarizability_domega(wavelength))

    def polarizability(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Float[np.ndarray, "*shape"]:
        """Return the real per-particle polarizability at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; ``inf`` is the static limit.
        warn : bool, optional
            Warn (:class:`~optmat.exceptions.ModelRangeWarning`) when a
            wavelength lies outside :attr:`wavelength_range`, and, for an absorbing
            model, (:class:`~optmat.exceptions.ModelAbsorptionWarning`) where the
            imaginary part is not negligible; True by default.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V (the real part), the shape of ``wavelength``.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.

        Warns
        -----
        ModelRangeWarning
            If ``warn`` and any wavelength lies outside the model's range.
        ModelAbsorptionWarning
            If ``warn``, the model absorbs and ``|Im alpha / Re alpha|`` exceeds
            :attr:`absorption_warning_threshold` at any wavelength.
        """
        wavelength = _validity.checked_wavelength(wavelength, self.name, static=True)
        if warn:
            self._check_range(wavelength)
        if not self.absorbing:
            return np.asarray(self._polarizability(wavelength))
        alpha = np.asarray(self._complex_polarizability(wavelength))
        if warn and np.any(
            np.abs(np.imag(alpha))
            > self.absorption_warning_threshold * np.abs(np.real(alpha))
        ):
            msg = (
                f"{self.name} ({self.gas}) absorbs at the requested wavelengths; "
                "the real part alone leaves the absorption out "
                "(use complex_polarizability or refractive_index)"
            )
            warnings.warn(msg, ModelAbsorptionWarning, stacklevel=2)
        return np.asarray(np.real(alpha))

    def refractivity(
        self,
        wavelength: ArrayLike | float,
        state: State = ONE_ATM_0C,
        *,
        ideal: bool = False,
        warn: bool = True,
    ) -> Float[np.ndarray, "*shape"]:
        """Return ``n - 1`` at a state, through Lorentz-Lorenz and the virial term.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State, optional
            Pressure and temperature; 0 degrees C and 1 atm by default.
        ideal : bool, optional
            Use the ideal-gas density instead of the real-gas one (always, for a
            species without an equation of state: ozone, a custom species).
        warn : bool, optional
            Warn outside the model's range (see :meth:`polarizability`).

        Returns
        -------
        numpy.ndarray
            ``n - 1``, the shape of ``wavelength``.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        """
        alpha = self.polarizability(wavelength, warn=warn)
        density = self._density(state, ideal=ideal)
        return refractivity(alpha, density, self.refractivity_virial)

    def refractive_index(
        self,
        wavelength: ArrayLike | float,
        state: State = ONE_ATM_0C,
        *,
        ideal: bool = False,
        warn: bool = True,
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex refractive index at a state (Lorentz-Lorenz, virial).

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State, optional
            Pressure and temperature; 0 degrees C and 1 atm by default.
        ideal : bool, optional
            Use the ideal-gas density instead of the real-gas one (always, for a
            species without an equation of state: ozone, a custom species).
        warn : bool, optional
            Warn outside the model's range (see :meth:`complex_polarizability`).

        Returns
        -------
        numpy.ndarray
            ``n`` (complex; ``Im n > 0`` is absorption, zero for a transparent model),
            the shape of ``wavelength``.

        Raises
        ------
        ArgumentError
            If a wavelength is zero or negative.
        """
        alpha = self.complex_polarizability(wavelength, warn=warn)
        density = self._density(state, ideal=ideal)
        return refractive_index(alpha, density, self.refractivity_virial)

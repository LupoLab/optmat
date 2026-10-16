# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The solid model interface: the complex dielectric function against vacuum wavelength.

A solid's model gives its complex relative permittivity epsilon = (n + ik)^2, with
k >= 0 for absorption, and its exact derivative with respect to angular frequency. The
refractive index, the group index and the group-velocity dispersion follow. The
interface matches what a propagation code needs of a cladding or a bulk medium: n(omega)
and epsilon(omega) with their omega derivatives, a name, a citation and a stated range.
"""

import abc
import math
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any, ClassVar

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat import _validity
from optmat.constants import C
from optmat.exceptions import ModelRangeWarning
from optmat.identity import Identified

__all__ = ["Model", "complex_index"]


def complex_index(
    permittivity: Complex[np.ndarray, "*shape"],
) -> Complex[np.ndarray, "*shape"]:
    """Return n + ik = sqrt(epsilon) with k >= 0, to full relative precision.

    numpy's complex square root returns k = 0 when k is below about 1e-16 n, which
    loses the weak absorption of a transparent window (TC-321); here the smaller of the
    two parts is formed from the larger as epsilon_i / (2 x).

    Parameters
    ----------
    permittivity : numpy.ndarray
        Complex relative permittivity, epsilon_i >= 0.

    Returns
    -------
    numpy.ndarray
        The complex index with Re n >= 0 and Im n >= 0.

    Examples
    --------
    >>> complex_index(np.array([2.25 + 3e-30j]))
    array([1.5+1.e-30j])
    """
    eps = np.asarray(permittivity, dtype=complex)
    modulus = np.abs(eps)
    re, im = np.real(eps), np.imag(eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        n_first = np.sqrt((modulus + re) / 2)
        k_first = np.sqrt(np.maximum(modulus - re, 0.0) / 2)
        n = np.where(re >= 0, n_first, np.where(k_first > 0, im / (2 * k_first), 0.0))
        k = np.where(re >= 0, np.where(n_first > 0, im / (2 * n_first), 0.0), k_first)
    return np.asarray(n + 1j * k)


@dataclass(frozen=True)
class Model(Identified):
    """A model of a solid's complex dielectric function.

    Attributes
    ----------
    material : str
        The material's id (:mod:`optmat.materials`).
    name : str
        A unique name, ``<source>/<formula>`` for a published formula.
    citation : str
        Where the model comes from.
    wavelength_range : tuple of float
        The vacuum-wavelength range over which the source supports the model, m: the
        envelope of :attr:`intervals`.
    kind : str
        The model's provenance: "published" (a formula or model as its authors give
        it), "fitted" (fitted here), "tabulated" (curated tables joined by range and
        interpolated) or "comparison" (a database entry, compared only). It is not
        ``identity.kind``, which names the quantity (``"solid"``).
    axis : str
        For an anisotropic crystal, the axis the model describes; empty otherwise.
    intervals : tuple of tuple of float
        The qualified intervals, vacuum wavelength, m; one, the range, unless the
        source qualifies several.
    status : str
        ``validated`` for a material's recommended model, ``comparison`` for the
        other curated models, ``user`` for a model a caller built.
    """

    IDENTITY_KIND: ClassVar[str] = "solid"
    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = ("intervals",)
    DERIVED: ClassVar[tuple[str, ...]] = ("wavelength_range",)

    material: str
    name: str
    citation: str
    wavelength_range: tuple[float, float]
    kind: str
    axis: str = ""
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
        """Return the model's name, material and source records.

        Returns
        -------
        tuple
            ``(name, material, sources)``.
        """
        return self.name, self.material, self._sources()

    def _sources(self) -> tuple[str, ...]:
        """Return the source records the model rests on.

        Returns
        -------
        tuple of str
            The first part of a ``<source>/<formula>`` name, where it is one.
        """
        head = self.name.split("/", 1)[0]
        return (head,) if "/" in self.name else ()

    @abc.abstractmethod
    def _permittivity(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return epsilon at angular frequencies, without a range check.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s, > 0.

        Returns
        -------
        numpy.ndarray
            The complex relative permittivity.
        """

    @abc.abstractmethod
    def _dpermittivity_domega(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d epsilon / d omega, without a range check.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s, > 0.

        Returns
        -------
        numpy.ndarray
            The derivative, s/rad.
        """

    def _omega(self, wavelength: ArrayLike | float, *, warn: bool) -> np.ndarray:
        """Convert vacuum wavelengths to angular frequencies, warning out of range.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool
            Warn when a wavelength lies outside :attr:`wavelength_range`.

        Returns
        -------
        numpy.ndarray
            Angular frequency, rad/s.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        wavelength = _validity.checked_wavelength(wavelength, self.name, static=False)
        if warn and np.any(_validity.outside(wavelength, self.intervals)):
            msg = (
                f"{self.name} ({self.material}) evaluated outside "
                f"{_validity.describe(self.intervals)}, the range its source supports"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)
        return np.asarray(2 * math.pi * C / wavelength)

    def permittivity(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex relative permittivity at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn (:class:`~optmat.exceptions.ModelRangeWarning`) outside
            :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            The permittivity epsilon = (n + ik)^2, with Im epsilon >= 0 for absorption.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        return np.asarray(self._permittivity(self._omega(wavelength, warn=warn)))

    def dpermittivity_domega(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d epsilon / d omega at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            The derivative with respect to angular frequency, s/rad.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        return np.asarray(
            self._dpermittivity_domega(self._omega(wavelength, warn=warn))
        )

    def refractive_index(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex refractive index n + ik at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            The index n + ik, with k >= 0 for absorption (the intensity falls as
            exp(-4 pi k z / lambda)).

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        return complex_index(self.permittivity(wavelength, warn=warn))

    def dindex_domega(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d(n + ik) / d omega at vacuum wavelengths.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            The derivative, s/rad.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        omega = self._omega(wavelength, warn=warn)
        index = complex_index(self._permittivity(omega))
        return np.asarray(self._dpermittivity_domega(omega) / (2 * index))

    def group_index(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Float[np.ndarray, "*shape"]:
        """Return the group index n + omega dn/domega (its real part).

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            The group index.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        omega = self._omega(wavelength, warn=warn)
        index = complex_index(self._permittivity(omega))
        slope = self._dpermittivity_domega(omega) / (2 * index)
        return np.asarray(np.real(index + omega * slope))

    def gvd(
        self, wavelength: ArrayLike | float, *, warn: bool = True
    ) -> Float[np.ndarray, "*shape"]:
        """Return the group-velocity dispersion beta_2 of the bulk medium.

        beta_2 = (2 dn/domega + omega d^2n/domega^2) / c, with the second derivative
        from a fourth-order central difference of the exact first derivative (step
        1e-3 omega: truncation below 1e-11 relative where the model is smooth).

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        warn : bool, optional
            Warn outside :attr:`wavelength_range`; True by default.

        Returns
        -------
        numpy.ndarray
            Re beta_2, s^2/m.

        Raises
        ------
        ArgumentError
            If a wavelength is not positive and finite.
        """
        omega = self._omega(wavelength, warn=warn)

        def slope(w: np.ndarray) -> np.ndarray:
            """Return dn/domega at w.

            Parameters
            ----------
            w : numpy.ndarray
                Angular frequency, rad/s.

            Returns
            -------
            numpy.ndarray
                The complex derivative, s/rad.
            """
            w = np.asarray(w)
            return np.asarray(
                self._dpermittivity_domega(w)
                / (2 * complex_index(self._permittivity(w)))
            )

        step = 1e-3 * omega
        curvature = (
            -slope(np.asarray(omega + 2 * step))
            + 8 * slope(np.asarray(omega + step))
            - 8 * slope(np.asarray(omega - step))
            + slope(np.asarray(omega - 2 * step))
        ) / (12 * step)
        return np.asarray(np.real(2 * slope(omega) + omega * curvature) / C)

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Oscillator-strength (pseudo-spectral) polarizability models.

The dynamic dipole polarizability of an atom or molecule below its first resonance is a
sum over its dipole oscillator-strength distribution. Represented by a finite set of
lines, or *pseudo-states*, of energies :math:`E_i` and oscillator strengths :math:`f_i`,

.. math::

   \alpha(\omega) = \sum_i \frac{f_i}{E_i^2 - \omega^2}

in atomic units (hartree for :math:`E_i` and :math:`\omega`, :math:`a_0^3` for
:math:`\alpha`). This is the form of Kumar and Meath's constrained distributions and,
in the dilute limit, of every Sellmeier formula; its Cauchy moments are
:math:`S(-2k) = \sum_i f_i E_i^{-2k}`, with :math:`\alpha(\omega) = \sum_{k \geq 1}
S(-2k)\,\omega^{2k-2}` below the lowest :math:`E_i`, and the Thomas-Reiche-Kuhn sum
:math:`S(0) = \sum_i f_i` equals the number of electrons for a complete distribution.

A model whose lowest pseudo-states sit at the measured resonance lines has physical
poles: it can be evaluated up to the first line and down to zero frequency without the
spurious structure a free Sellmeier fit may place in between.

The refractive index measures the *refractive* polarizability: the electric
polarizability plus the atom's (diamagnetic) magnetic susceptibility :math:`\chi_0`,
because :math:`n^2 = \varepsilon\mu`. The model carries :math:`\chi_0` as a
frequency-independent ``offset_au``; its sums :math:`S(k)` are those of the electric
distribution alone. Where it is known (helium), the retardation of the
forward-scattering amplitude adds a term :math:`c\,\omega^2` to the refractive
polarizability (``retardation_au``; Pachucki & Puchalski 2019, TC-705), which the sums
do not contain either.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_ENERGY, AU_POLARIZABILITY, C
from optmat.gas.absorption import MeasuredAbsorption
from optmat.gas.gases import Gas
from optmat.gas.models.base import Model

__all__ = ["OscillatorModel", "angular_frequency_au"]


def _species(name: Gas | str) -> Gas | str:
    """Return a curated gas by its name, or a custom species' name as given.

    Parameters
    ----------
    name : Gas or str
        The species.

    Returns
    -------
    Gas or str
        The :class:`Gas` member it names, else the name.
    """
    try:
        return Gas(name)
    except ValueError:
        return str(name)


def angular_frequency_au(
    wavelength: Float[np.ndarray, "*shape"],
) -> Float[np.ndarray, "*shape"]:
    """Return the photon energy in hartree of a vacuum wavelength.

    Parameters
    ----------
    wavelength : numpy.ndarray
        Vacuum wavelength, m (``inf`` gives 0).

    Returns
    -------
    numpy.ndarray
        ``omega`` in atomic units.

    Examples
    --------
    >>> f"{float(angular_frequency_au(np.array(45.5634e-9))):.5f}"
    '1.00000'
    """
    return np.asarray(2 * np.pi * C / wavelength / AU_ANGULAR_FREQUENCY)


@dataclass(frozen=True)
class OscillatorModel(Model):
    """A polarizability model made of oscillator strengths at fixed energies.

    :meth:`from_parameters` builds one from a caller's poles, in SI; the overrides
    (:meth:`~optmat.identity.Identified.with_parameters`) name the fields, in atomic
    units: :attr:`INPUTS`.

    Attributes
    ----------
    strengths : tuple of float
        Oscillator strengths ``f_i`` (dimensionless).
    energies : tuple of float
        Energies ``E_i``, hartree, > 0.
    labels : tuple of str
        What each term represents ("He 1s2p 1P", "effective continuum", ...).
    offset_au : float
        A frequency-independent term added to the sum, atomic units: the magnetic
        susceptibility ``chi_0`` that makes the sum's electric polarizability the
        refractive one.
    retardation_au : float
        The coefficient ``c`` of a term ``c omega^2`` (omega in hartree) added to the
        refractive polarizability, atomic units: the retardation (quadrupole and
        magnetic-dipole) part of the forward-scattering amplitude, known for helium
        (Pachucki & Puchalski 2019; TC-705). Zero for the other gases.
    absorption : MeasuredAbsorption or None
        A measured absorption spectrum whose complex polarizability (Kramers-Kronig)
        is added to the sum: the model then absorbs (ozone). None for a transparent
        model.
    """

    strengths: tuple[float, ...] = field(default=(), repr=False)
    energies: tuple[float, ...] = field(default=(), repr=False)
    labels: tuple[str, ...] = field(default=(), repr=False)
    offset_au: float = 0.0
    retardation_au: float = 0.0
    absorption: MeasuredAbsorption | None = field(default=None, repr=False)

    # 2: the first frequency derivative of a measured absorption term is the spline's
    # own, not a gradient of the grid (OEP 1, section 4.1; TC-732).
    EVALUATION_VERSION: ClassVar[int] = 2
    INPUTS: ClassVar[tuple[str, ...]] = (
        "energies",
        "intervals",
        "labels",
        "offset_au",
        "retardation_au",
        "strengths",
    )

    @classmethod
    def from_parameters(  # noqa: PLR0913 - the model's parameters
        cls,
        species: Gas | str,
        strengths: ArrayLike,
        energies: ArrayLike,
        *,
        intervals: Sequence[tuple[float, float]],
        offset: float = 0.0,
        labels: Sequence[str] | None = None,
        citation: str = "",
        name: str = "",
    ) -> "OscillatorModel":
        """Return a gas polarizability as a sum of poles, from SI parameters.

        ``alpha(omega) = (e^2/m_e) sum_k f_k / (omega_k^2 - omega^2) + offset``, the
        form of OEP 2's compact models, with ``E_k = hbar omega_k``.

        Parameters
        ----------
        species : Gas or str
            A curated gas, or a custom species' name, the identity's material.
        strengths : array_like
            The oscillator strengths ``f_k``, dimensionless.
        energies : array_like
            The poles' photon energies ``E_k``, J, > 0.
        intervals : sequence of (float, float)
            The qualified intervals, vacuum wavelength, m, disjoint and ascending.
        offset : float, optional
            A frequency-independent term, C m^2/V (a magnetic susceptibility).
        labels : sequence of str or None, optional
            What each pole is; ``pole 1``, ... by default.
        citation : str, optional
            Where the parameters come from.
        name : str, optional
            The model's name; ``user/<species>`` by default.

        Returns
        -------
        OscillatorModel
            The model, status ``user``.

        Raises
        ------
        ValueError
            If the lengths differ, an energy is not positive or the intervals are not
            disjoint and ascending.

        Examples
        --------
        >>> from optmat.constants import E_CHARGE
        >>> model = OscillatorModel.from_parameters(
        ...     "Ar", [1.2, 4.0], [11.8 * E_CHARGE, 20.0 * E_CHARGE],
        ...     intervals=[(150e-9, float("inf"))])
        >>> model.identity.status, model.name
        ('user', 'user/Ar')
        """
        spans = tuple((float(a), float(b)) for a, b in intervals)
        return cls(
            gas=_species(species),
            name=name or f"user/{species}",
            citation=citation,
            wavelength_range=(spans[0][0], spans[-1][1]),
            kind="user",
            intervals=spans,
            strengths=tuple(float(f) for f in np.ravel(strengths)),
            energies=tuple(float(e) / AU_ENERGY for e in np.ravel(energies)),
            labels=tuple(
                labels or (f"pole {k + 1}" for k in range(np.size(strengths)))
            ),
            offset_au=offset / AU_POLARIZABILITY,
        )

    @property
    def absorbing(self) -> bool:
        """Whether the model has a measured absorption term.

        Returns
        -------
        bool
            True when :attr:`absorption` is set.
        """
        return self.absorption is not None

    def __post_init__(self) -> None:
        """Check that the terms agree in number and the energies are positive.

        Raises
        ------
        ValueError
            If the lengths differ or an energy is not positive.
        """
        super().__post_init__()
        if not len(self.strengths) == len(self.energies) == len(self.labels):
            msg = f"{self.name}: strengths, energies and labels differ in length"
            raise ValueError(msg)
        if any(e <= 0 for e in self.energies):
            msg = f"{self.name}: every energy must be positive"
            raise ValueError(msg)

    def polarizability_au(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Return the refractive polarizability in atomic units.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energy, hartree.

        Returns
        -------
        numpy.ndarray
            ``sum_i f_i / (E_i^2 - omega^2) + offset_au + retardation_au omega^2``,
            atomic units.
        """
        omega2 = np.asarray(omega, dtype=float) ** 2
        measured = (
            0.0
            if self.absorption is None
            else self.absorption.complex_polarizability_au(np.asarray(np.sqrt(omega2)))[
                0
            ]
        )
        return np.asarray(
            sum(
                f / (e**2 - omega2)
                for f, e in zip(self.strengths, self.energies, strict=True)
            )
            + self.offset_au
            + self.retardation_au * omega2
            + measured
        )

    def dpolarizability_domega_au(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Return the derivative of the polarizability with frequency, analytically.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energy, hartree.

        Returns
        -------
        numpy.ndarray
            ``sum_i 2 f_i omega / (E_i^2 - omega^2)^2 + 2 retardation_au omega``,
            atomic units (a0^3 per hartree).

        Examples
        --------
        >>> from optmat.gas.gases import Gas
        >>> model = OscillatorModel(
        ...     gas=Gas.He, name="demo", citation="", kind="fitted",
        ...     wavelength_range=(60e-9, float("inf")),
        ...     strengths=(0.5,), energies=(0.8,), labels=("a",),
        ... )
        >>> omega, h = np.array([0.3]), 1e-6
        >>> numeric = (model.polarizability_au(omega + h)
        ...            - model.polarizability_au(omega - h)) / (2 * h)
        >>> bool(abs(model.dpolarizability_domega_au(omega)[0] / numeric[0] - 1) < 1e-8)
        True
        """
        omega = np.asarray(omega, dtype=float)
        measured = (
            0.0
            if self.absorption is None
            else self.absorption.complex_derivative_au(omega)[0]
        )
        return np.asarray(
            sum(
                2 * f * omega / (e**2 - omega**2) ** 2
                for f, e in zip(self.strengths, self.energies, strict=True)
            )
            + 2 * self.retardation_au * omega
            + measured
        )

    def dispersion_au(
        self, omega: Float[np.ndarray, "*shape"], order: int
    ) -> Float[np.ndarray, "*shape"]:
        """Return the ``order``-th derivative of ``omega alpha(omega)``, analytically.

        In the dilute limit ``beta = (omega / c) [1 + N alpha / (2 epsilon_0)]``, so
        for ``order >= 2`` the dispersion coefficient ``beta_k`` is proportional to
        this derivative with a factor that depends only on the density: its ratio
        between models, or its relative change, needs no state. The magnetic term
        contributes only to orders 0 and 1, the retardation term ``c omega^3`` to orders
        0 to 3.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energy, hartree, below the lowest term's energy.
        order : int
            The order ``k >= 0``.

        Returns
        -------
        numpy.ndarray
            ``d^k (omega alpha) / d omega^k``, atomic units.

        Raises
        ------
        ValueError
            If ``order`` is negative.

        Examples
        --------
        >>> from optmat.gas.gases import Gas
        >>> model = OscillatorModel(
        ...     gas=Gas.He, name="demo", citation="", kind="fitted",
        ...     wavelength_range=(60e-9, float("inf")),
        ...     strengths=(0.5,), energies=(0.8,), labels=("a",),
        ... )
        >>> omega = np.array([0.3])
        >>> first = model.dispersion_au(omega, 1)
        >>> direct = (model.polarizability_au(omega)
        ...           + omega * model.dpolarizability_domega_au(omega))
        >>> bool(abs(first[0] / direct[0] - 1) < 1e-12)
        True
        """
        if order < 0:
            msg = f"order must be non-negative, got {order}"
            raise ValueError(msg)
        omega = np.asarray(omega, dtype=float)
        # omega f / (E^2 - omega^2) = (f / 2) [1 / (E - omega) - 1 / (E + omega)].
        factor, sign = 0.5 * math.factorial(order), (-1) ** order
        total = sum(
            f
            * factor
            * ((e - omega) ** -(order + 1) - sign * (e + omega) ** -(order + 1))
            for f, e in zip(self.strengths, self.energies, strict=True)
        )
        magnetic = {0: self.offset_au * omega, 1: self.offset_au}.get(order, 0.0)
        # d^k (c omega^3) / d omega^k.
        c = self.retardation_au
        retardation = {0: c * omega**3, 1: 3 * c * omega**2, 2: 6 * c * omega}.get(
            order, 6 * c if order == 3 else 0.0
        )
        measured = 0.0
        if self.absorption is not None:
            # d^k (omega a) / d omega^k = omega a^(k) + k a^(k-1). The first order
            # uses the spline's own derivative (OEP 1, section 4.1); the higher ones
            # stay gradients of the grid, since a cubic spline's third derivative is
            # piecewise constant.
            if order == 1:
                derivative = self.absorption.complex_derivative_au(omega)[0]
            else:
                derivative = self.absorption.real_derivative_au(omega, order)
            measured = omega * derivative
            if order > 0:
                measured = measured + order * self.absorption.real_derivative_au(
                    omega, order - 1
                )
        return np.asarray(total + magnetic + retardation + measured)

    def _polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Float[np.ndarray, "*shape"]:
        """Evaluate the model at vacuum wavelengths.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V.
        """
        return np.asarray(
            self.polarizability_au(angular_frequency_au(wavelength)) * AU_POLARIZABILITY
        )

    def _complex_polarizability(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Evaluate the complex polarizability at vacuum wavelengths.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            Polarizability, C m^2/V (complex; the imaginary part is the measured
            absorption's, zero without it).
        """
        real = self._polarizability(wavelength)
        if self.absorption is None:
            return np.asarray(real, dtype=complex)
        omega = np.asarray(angular_frequency_au(wavelength))
        imaginary = self.absorption.complex_polarizability_au(omega)[1]
        return np.asarray(real + 1j * imaginary * AU_POLARIZABILITY)

    def _dcomplex_polarizability_domega(
        self, wavelength: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d alpha / d omega at vacuum wavelengths, analytically.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelength, m (``inf`` for the static limit, where it is 0).

        Returns
        -------
        numpy.ndarray
            C m^2 V^-1 s (complex; the imaginary part is the measured absorption's).
        """
        omega = np.asarray(angular_frequency_au(wavelength))
        real = self.dpolarizability_domega_au(omega)
        imaginary = (
            0.0
            if self.absorption is None
            else self.absorption.complex_derivative_au(omega)[1]
        )
        scale = AU_POLARIZABILITY / AU_ANGULAR_FREQUENCY
        return np.asarray((real + 1j * imaginary) * scale, dtype=complex)

    def cauchy_moment(self, k: int) -> float:
        """Return the electric distribution's dipole sum ``S(k) = sum_i f_i E_i^k``.

        Parameters
        ----------
        k : int
            The order: ``S(-2)`` is the static polarizability, ``S(0)`` the
            Thomas-Reiche-Kuhn sum, ``S(-4)`` the first dispersion coefficient.

        Returns
        -------
        float
            The sum, atomic units.

        Examples
        --------
        >>> from optmat.gas.gases import Gas
        >>> model = OscillatorModel(
        ...     gas=Gas.He, name="demo", citation="", kind="fitted",
        ...     wavelength_range=(60e-9, float("inf")),
        ...     strengths=(0.5, 1.5), energies=(0.8, 1.2), labels=("a", "b"),
        ... )
        >>> model.cauchy_moment(-2) == float(model.polarizability_au(np.array(0.0)))
        True

        (``offset_au`` is 0 here; otherwise the static refractive polarizability is
        ``S(-2) + offset_au``. The refractive polarizability's omega^2 coefficient is
        ``S(-4) + retardation_au``.)
        """
        measured = 0.0 if self.absorption is None else self.absorption.sum_rule(k)
        return float(
            sum(f * e**k for f, e in zip(self.strengths, self.energies, strict=True))
            + measured
        )

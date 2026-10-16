# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The refractive index of a gas mixture at finite density (OEP 4, sections 4.1, 4.3).

In the low-density cluster expansion the Lorentz-Lorenz function of a mixture is

.. math::

   y = \frac{n^2 - 1}{n^2 + 2} = \sum_s A_s \rho_s + \sum_{s,t} B_{st}\,\rho_s \rho_t,
   \qquad A_s = \frac{N_A \alpha_s}{3\varepsilon_0},

with :math:`\rho_s = N_s/N_A` the molar densities and :math:`B_{ss}` each gas's own
refractivity virial coefficient, so that one species is optmat's pure law exactly. The
whole mixture's :math:`y` is inverted once, :math:`\chi = 3y/(1-y)`.

- **Like terms** disperse with their gas, :math:`B_{ss} = b_s A_s`, as in the pure
  models. O2's is a published estimate, ``provisional`` with its 50 % uncertainty. A
  gas whose :math:`B_R` optmat has not curated (D2, N2O, CO, O3) is
  ``unknown``: it enters with :math:`b = 0` and the declared scale
  :math:`|b| \le 2.5` cm^3/mol.
- **Unlike terms** are not curated for any pair. Each is estimated as the mean of its
  like terms, :math:`B_{st} = \tfrac12(B_{ss} + B_{tt})`, each like term with its own
  :math:`A`, with standard uncertainty :math:`\max(|B_{ss}|, |B_{tt}|)`, taking an
  unknown like term at its scale. The status is ``estimate``, or ``unknown`` where a
  like term is unknown (section 7.2).
- **The uncertainty** of the estimated and unknown terms is combined linearly, as fully
  correlated (:meth:`MixtureOptics.refractivity_uncertainty`).

**The optical domain** (:meth:`MixtureOptics.density_bound`) is the lower of two
densities. One is where the neglected third term reaches 1e-3 of :math:`y`, from the
largest constituent :math:`|C_R/A_R|`: measured, a measured bound, or the declared scale
:math:`600 B_R^2/A_R`. The other is where the estimated and unknown terms' uncertainty
reaches 1e-3 of :math:`y`. Beyond the bound the optics warn. In wavelength, a mixture
holds on the intersection of its constituents' qualified intervals.
"""

import functools
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat import _validity
from optmat.constants import EPS0, N_A
from optmat.exceptions import ModelRangeWarning, OptmatError
from optmat.gas.air import _intersection  # pyright: ignore[reportPrivateUsage]
from optmat.gas.gases import Gas, as_gas
from optmat.gas.models.base import Model
from optmat.gas.models.registry import default_model
from optmat.gas.virial import like_term, pair_rules, third_term
from optmat.identity import Identified

__all__ = ["MixtureOptics", "PairTerm", "mixture_optics"]

_MOLAR: Final = N_A / (3 * EPS0)
"""``A / alpha = N_A / (3 eps0)``, m^3/mol per C m^2/V."""

_BOUND_POINTS: Final = 48
"""Wavelengths per qualified interval at which the optical bound is taken."""


@dataclass(frozen=True)
class PairTerm:
    """One pair's refractivity virial term, like or unlike (OEP 4, section 6).

    Attributes
    ----------
    species, partner : str
        The pair, in sorted order; the same species twice for a like term.
    b : float
        ``B_st / mean(A_s, A_t)``, m^3/mol: a like term's curated ratio, or 0 for an
        unknown one. NaN for an unlike pair, whose rule is named in ``source``.
    sigma : float
        The standard uncertainty of ``b``, m^3/mol: the curated one, or the declared
        scale for an unknown like term. NaN for an unlike pair (see ``reason``).
    status : str
        ``curated`` (a like term optmat curates), ``estimate`` (the mean rule) or
        ``unknown`` (a like term optmat has not curated, and its pairs); a like
        term that is a published estimate is ``provisional``.
    source : str
        The curated dataset, or ``rule: mean``.
    reason : str
        Why the term is what it is.
    frequency : str
        ``static`` or the source's vacuum wavelength, nm; empty where there is none.
    temperature : float
        The source's temperature, K; NaN where there is none.
    contribution : str
        ``induced``, or ``total`` for a static value that includes orientation.
    """

    species: str
    partner: str
    b: float
    sigma: float
    status: str
    source: str
    reason: str
    frequency: str
    temperature: float
    contribution: str


def _pairs(species: tuple[str, ...]) -> tuple[PairTerm, ...]:
    """Return every pair of a mixture's species, like terms first.

    Parameters
    ----------
    species : tuple of str
        The species, sorted.

    Returns
    -------
    tuple of PairTerm
        One per unordered pair.
    """
    rules = pair_rules()
    like = {s: like_term(s) for s in species}
    out = [
        PairTerm(
            species=s,
            partner=s,
            b=like[s].b,
            sigma=like[s].sigma,
            status=like[s].status,
            source=like[s].source,
            reason={
                "curated": "the B_R chosen in refractivity_virial.toml",
                "provisional": "a published estimate with its stated uncertainty "
                "(refractivity_virial.toml)",
                "unknown": "not curated: b = 0 with the declared scale "
                "(OEP 4, section 7.2)",
            }[like[s].status],
            frequency=like[s].frequency,
            temperature=like[s].temperature,
            contribution=like[s].contribution,
        )
        for s in species
    ]
    for i, s in enumerate(species):
        for t in species[i + 1 :]:
            known = "unknown" not in (like[s].status, like[t].status)
            out.append(
                PairTerm(
                    species=s,
                    partner=t,
                    b=math.nan,
                    sigma=math.nan,
                    status="estimate" if known else "unknown",
                    source=f"rule: {rules['rule']}",
                    reason=str(rules["reason"]),
                    frequency="",
                    temperature=math.nan,
                    contribution="induced",
                )
            )
    return tuple(out)


def _bound_wavelengths(intervals: tuple[tuple[float, float], ...]) -> np.ndarray:
    """Return the wavelengths at which the optical bound is taken.

    Parameters
    ----------
    intervals : tuple of tuple of float
        The mixture's qualified intervals, m.

    Returns
    -------
    numpy.ndarray
        :data:`_BOUND_POINTS` geometric points across each interval (to 1 mm where it
        is open), and the static limit where an interval is open.
    """
    points = []
    for low, high in intervals:
        points.append(np.geomspace(low, min(high, 1e-3), _BOUND_POINTS))
        if math.isinf(high):
            points.append(np.array([math.inf]))
    return np.concatenate(points)


@dataclass(frozen=True, eq=False)
class MixtureOptics(Identified):
    """The optics of a gas mixture at finite density (OEP 4).

    Attributes
    ----------
    species : tuple of str
        The species, sorted.
    models : tuple of Model
        Each species' recommended dispersion model, in the same order.
    pairs : tuple of PairTerm
        Every like and unlike pair, like terms first.
    intervals : tuple of tuple of float
        The intersection of the models' qualified intervals, m.
    status : str
        ``provisional`` where any like term is ``unknown`` or ``provisional``, else
        ``validated``.
    """

    IDENTITY_KIND: ClassVar[str] = "gas_refractive_mixture"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ()

    species: tuple[str, ...]
    models: tuple[Model, ...]
    pairs: tuple[PairTerm, ...]
    intervals: tuple[tuple[float, float], ...]

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the mixture's name, material and source records.

        Returns
        -------
        tuple
            ``("mixture/<species joined by +>", the same, sources)``.
        """
        material = "+".join(self.species)
        found = {p.source.split("/", 1)[0] for p in self.pairs if "/" in p.source}
        return f"mixture/{material}", material, tuple(sorted(found))

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the pair rules and the third terms, which the bound reads.

        Returns
        -------
        Mapping
            The ``[pairs]`` rules and each species' third-term status and ratio.
        """
        return {
            "pair_rules": dict(pair_rules()),
            "third_terms": {
                s: [third_term(s).status, third_term(s).ratio] for s in self.species
            },
        }

    @functools.cached_property
    def _b(self) -> dict[str, float]:
        """Each species' like ratio ``b``, m^3/mol (0 for an unknown one).

        Returns
        -------
        dict
            Species to ``b``.
        """
        return {p.species: p.b for p in self.pairs if p.species == p.partner}

    def _molar(
        self,
        wavelength: np.ndarray,
        *,
        complex_: bool,
        derivative: bool = False,
        absorption: bool = False,
    ) -> dict[str, Any]:
        """Return each species' ``A`` (or ``dA/d omega``) at wavelengths.

        The constituents' own range warnings are silenced: the mixture warns once, for
        its own intervals.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelengths, m.
        complex_ : bool
            From the complex polarizability, or the real.
        derivative : bool, optional
            ``dA/d omega`` (always complex) in place of ``A``.
        absorption : bool, optional
            For the real ``A``, warn where a constituent absorbs, as its model does.

        Returns
        -------
        dict
            Species to m^3/mol (or m^3 s/mol).
        """
        out = {}
        for s, model in zip(self.species, self.models, strict=True):
            if derivative:
                alpha = model.dpolarizability_domega(wavelength, warn=False)
            elif complex_:
                alpha = model.complex_polarizability(wavelength, warn=False)
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ModelRangeWarning)
                    alpha = model.polarizability(wavelength, warn=absorption)
            out[s] = _MOLAR * np.asarray(alpha)
        return out

    def _unlike(self, molar: Mapping[str, Any], s: str, t: str) -> np.ndarray:
        """Return the estimate ``B_st``: the mean of the like terms, each with its A.

        Parameters
        ----------
        molar : Mapping
            Species to ``A`` (or ``dA/d omega``).
        s, t : str
            The pair.

        Returns
        -------
        numpy.ndarray
            ``(b_s A_s + b_t A_t) / 2``.
        """
        return np.asarray(0.5 * (self._b[s] * molar[s] + self._b[t] * molar[t]))

    def _y(self, molar: Mapping[str, Any], rho: Mapping[str, Any]) -> np.ndarray:
        """Return the Lorentz-Lorenz function of the mixture.

        It is linear in the ``A`` at fixed ``b``, so ``dA/d omega`` in place of ``A``
        gives ``dy/d omega``.

        Parameters
        ----------
        molar : Mapping
            Species to ``A``.
        rho : Mapping
            Species to molar density.

        Returns
        -------
        numpy.ndarray
            ``y``, broadcast over wavelength and density.
        """
        y: Any = 0.0
        for i, s in enumerate(self.species):
            # A rho (1 + b rho), as the pure law writes it.
            y = y + molar[s] * rho[s] * (1 + self._b[s] * rho[s])
            for t in self.species[i + 1 :]:
                y = y + 2 * self._unlike(molar, s, t) * rho[s] * rho[t]
        return np.asarray(y)

    def _spread(self, molar: Mapping[str, Any]) -> dict[tuple[str, str], Any]:
        """Return each estimated or unknown term's standard uncertainty, as a ``B``.

        Parameters
        ----------
        molar : Mapping
            Species to real ``A``.

        Returns
        -------
        dict
            ``(s, t)``, s <= t, to its uncertainty, m^6/mol^2. Curated like terms are
            left out: their uncertainty is the source's, not a sensitivity. A
            provisional like term enters with its stated uncertainty, an unknown one
            with its scale.
        """
        scale = {
            p.species: p.sigma if p.status == "unknown" else abs(p.b)
            for p in self.pairs
            if p.species == p.partner
        }
        out: dict[tuple[str, str], Any] = {}
        for p in self.pairs:
            s, t = p.species, p.partner
            if s == t:
                if p.status != "curated":
                    out[s, s] = p.sigma * np.abs(molar[s])
            else:
                out[s, t] = np.maximum(
                    scale[s] * np.abs(molar[s]), scale[t] * np.abs(molar[t])
                )
        return out

    def _densities(self, densities: Mapping[Any, ArrayLike]) -> dict[str, np.ndarray]:
        """Return the molar densities by species, checking the keys.

        Parameters
        ----------
        densities : Mapping
            Species (a :class:`~optmat.gas.gases.Gas` or its name) to number density,
            m^-3.

        Returns
        -------
        dict
            Species to molar density, mol/m^3.

        Raises
        ------
        OptmatError
            If a species is missing or extra, or a density is negative.
        """
        # A Gas and its name are one key (Gas is a str enum), so no mapping repeats a
        # species; two labelled populations are combined by the caller (section 4.1).
        given = {
            as_gas(k).value: np.asarray(v, dtype=float) for k, v in densities.items()
        }
        missing = sorted(set(self.species) - set(given))
        extra = sorted(set(given) - set(self.species))
        if missing or extra:
            msg = (
                f"the mixture of {', '.join(self.species)} takes exactly their "
                f"densities: missing {missing}, extra {extra}"
            )
            raise OptmatError(msg)
        if any(np.any(v < 0) for v in given.values()):
            msg = "number densities must not be negative"
            raise OptmatError(msg)
        return {s: np.asarray(given[s] / N_A) for s in self.species}

    @functools.cached_property
    def _bound_terms(
        self,
    ) -> tuple[dict[str, Any], dict[tuple[str, str], Any]]:
        """The real ``A`` and the uncertainties at the bound's wavelengths.

        Returns
        -------
        tuple
            ``(A by species, uncertainty by pair)``, each over the wavelengths of
            :func:`_bound_wavelengths`.
        """
        molar = self._molar(_bound_wavelengths(self.intervals), complex_=False)
        return molar, self._spread(molar)

    def _molar_bound(self, x: Mapping[str, Any]) -> np.ndarray:
        """Return the optical bound on the total molar density at compositions.

        Parameters
        ----------
        x : Mapping
            Species to mole fraction, each broadcastable to the others.

        Returns
        -------
        numpy.ndarray
            The bound, mol/m^3; inf where nothing bounds it.
        """
        fraction = float(pair_rules()["density_fraction"])
        ratio: Any = 0.0
        for s in self.species:
            third = third_term(s).ratio
            if not math.isnan(third):
                ratio = np.maximum(ratio, np.where(x[s] > 0, third, 0.0))
        with np.errstate(divide="ignore"):
            third_bound = np.sqrt(fraction / np.asarray(ratio))
        molar, spread = self._bound_terms
        # The bound's wavelengths are the leading axis; the smallest is taken.
        linear: Any = 0.0
        for s in self.species:
            linear = linear + np.multiply.outer(molar[s], x[s])
        quadratic: Any = 0.0
        for (s, t), value in spread.items():
            weight = 1 if s == t else 2
            quadratic = quadratic + weight * np.multiply.outer(value, x[s] * x[t])
        quadratic = np.broadcast_to(quadratic, np.shape(linear))
        with np.errstate(divide="ignore", invalid="ignore"):
            pairs = np.where(quadratic > 0, fraction * linear / quadratic, np.inf)
        return np.asarray(np.minimum(third_bound, np.min(pairs, axis=0)))

    def density_bound(self, composition: Mapping[Any, float]) -> float:
        """Return the optical domain's bound on the total number density.

        The lower of the density at which the third term reaches 1e-3 of ``y`` (the
        largest ``|C_R/A_R|`` of the species present) and that at which the estimated
        and unknown terms' uncertainty does, the smallest over the qualified intervals
        (section 4.3).

        Parameters
        ----------
        composition : Mapping
            Species to mole fraction, summing to one; a species left out is absent.

        Returns
        -------
        float
            The bound, m^-3; inf where nothing bounds it.

        Raises
        ------
        OptmatError
            For a species not in the mixture, a negative fraction, or fractions that do
            not sum to one.

        Examples
        --------
        >>> f"{mixture_optics(['Ar']).density_bound({'Ar': 1.0}) / 6.02214076e23:.0f}"
        '6961'
        """
        x = {as_gas(k).value: float(v) for k, v in composition.items()}
        extra = sorted(set(x) - set(self.species))
        if extra or any(v < 0 for v in x.values()):
            msg = f"not a composition of {', '.join(self.species)}: {composition}"
            raise OptmatError(msg)
        if not math.isclose(sum(x.values()), 1.0, rel_tol=1e-9):
            msg = f"mole fractions sum to {sum(x.values())}, not one"
            raise OptmatError(msg)
        fractions = {s: np.asarray(x.get(s, 0.0)) for s in self.species}
        return float(self._molar_bound(fractions)) * N_A

    def _check_range(self, wavelength: np.ndarray, rho: Mapping[str, Any]) -> None:
        """Warn outside the mixture's intervals, or beyond its optical domain.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Vacuum wavelengths, m.
        rho : Mapping
            Species to molar density.
        """
        if np.any(_validity.outside(wavelength, self.intervals)):
            msg = (
                f"mixture of {', '.join(self.species)} evaluated outside "
                f"{_validity.describe(self.intervals)}, its constituents' common range"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=4)
        total = np.asarray(sum(rho.values()))
        if not np.any(total > 0):
            return
        safe = np.where(total > 0, total, 1.0)
        bound = self._molar_bound({s: rho[s] / safe for s in self.species})
        if np.any(total > bound):
            msg = (
                f"mixture of {', '.join(self.species)} beyond its optical domain, "
                f"{float(np.min(bound)):.4g} mol/m^3 at its composition (OEP 4, "
                "section 4.3)"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=4)

    def _prepare(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        """Check the inputs, warn, and return the wavelengths and molar densities.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        tuple
            The wavelengths and the molar densities by species.
        """
        wl = np.asarray(wavelength, dtype=float)
        rho = self._densities(densities)
        self._check_range(wl, rho)
        return wl, rho

    def susceptibility(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex susceptibility ``chi = n^2 - 1`` of the mixture.

        For one species, ``sqrt(1 + chi)`` is the pure model's ``refractive_index``.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; broadcast with the densities.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        numpy.ndarray
            ``chi``, from each constituent's complex polarizability.

        Warns
        -----
        ModelRangeWarning
            Outside the constituents' common intervals or the optical domain.
        """
        wl, rho = self._prepare(wavelength, densities)
        y = self._y(self._molar(wl, complex_=True), rho)
        return np.asarray(3 * y / (1 - y))

    def refractivity(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> Float[np.ndarray, "*shape"]:
        """Return ``n - 1`` of the mixture, in the real-polarizability convention.

        For one species it is the pure model's ``refractivity`` at that density.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m; broadcast with the densities.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        numpy.ndarray
            ``n - 1 = chi / (1 + sqrt(1 + chi))``.

        Warns
        -----
        ModelRangeWarning
            Outside the constituents' common intervals or the optical domain.
        ModelAbsorptionWarning
            Where a constituent absorbs, as its model's ``polarizability`` does.

        Examples
        --------
        >>> he_ar = mixture_optics(["He", "Ar"])
        >>> n = 1.25e25  # m^-3 each: about 1 bar together at 294 K
        >>> f"{float(he_ar.refractivity(800e-9, {'He': n, 'Ar': n})):.5e}"
        '1.46154e-04'
        """
        wl, rho = self._prepare(wavelength, densities)
        y = self._y(self._molar(wl, complex_=False, absorption=True), rho)
        chi = 3 * y / (1 - y)
        return np.asarray(chi / (1 + np.sqrt(1 + chi)))

    def refractivity_uncertainty(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> Float[np.ndarray, "*shape"]:
        """Return the estimated and unknown terms' effect on ``n - 1``.

        Their uncertainties are combined linearly, as fully correlated (section 7.2):
        a sensitivity, not a statistical uncertainty.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        numpy.ndarray
            The change in ``n - 1`` they span.

        Warns
        -----
        ModelRangeWarning
            Outside the constituents' common intervals or the optical domain.
        """
        wl, rho = self._prepare(wavelength, densities)
        molar = self._molar(wl, complex_=False, absorption=True)
        y = self._y(molar, rho)
        dy: Any = 0.0
        for (s, t), value in self._spread(molar).items():
            dy = dy + (1 if s == t else 2) * value * rho[s] * rho[t]
        # d(n - 1)/dy, with n^2 = (1 + 2y)/(1 - y).
        slope = 1.5 / ((1 - y) ** 2 * np.sqrt((1 + 2 * y) / (1 - y)))
        return np.asarray(slope * dy)

    def dsusceptibility_domega(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return ``d chi / d omega`` of the mixture.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        numpy.ndarray
            The derivative, s.

        Warns
        -----
        ModelRangeWarning
            Outside the constituents' common intervals or the optical domain.
        """
        wl, rho = self._prepare(wavelength, densities)
        y = self._y(self._molar(wl, complex_=True), rho)
        dy = self._y(self._molar(wl, complex_=True, derivative=True), rho)
        return np.asarray(3 / (1 - y) ** 2 * dy)

    def dsusceptibility_ddensity(
        self, wavelength: ArrayLike | float, densities: Mapping[Any, ArrayLike]
    ) -> dict[str, Complex[np.ndarray, "*shape"]]:
        """Return ``d chi / d N_s`` for every species.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        densities : Mapping
            Each species' number density, m^-3.

        Returns
        -------
        dict
            Species to the derivative with respect to its number density, m^3.

        Warns
        -----
        ModelRangeWarning
            Outside the constituents' common intervals or the optical domain.
        """
        wl, rho = self._prepare(wavelength, densities)
        molar = self._molar(wl, complex_=True)
        factor = 3 / (1 - self._y(molar, rho)) ** 2 / N_A
        out = {}
        for s in self.species:
            pair: Any = 2 * self._b[s] * molar[s] * rho[s]
            for t in self.species:
                if t != s:
                    pair = pair + 2 * self._unlike(molar, s, t) * rho[t]
            out[s] = np.asarray(factor * (molar[s] + pair))
        return out


@functools.cache
def _optics(species: tuple[str, ...]) -> MixtureOptics:
    """Return the optics of sorted species, built once.

    Parameters
    ----------
    species : tuple of str
        The species, sorted and distinct.

    Returns
    -------
    MixtureOptics
        With each species' recommended model.
    """
    models = tuple(default_model(s) for s in species)
    pairs = _pairs(species)
    like = [p for p in pairs if p.species == p.partner]
    unknown = [p.species for p in like if p.status == "unknown"]
    estimated = [p.species for p in like if p.status == "provisional"]
    reasons = []
    if unknown:
        reasons.append(
            f"B_R not curated for {', '.join(unknown)}: its terms are unknown"
        )
    if estimated:
        reasons.append(f"B_R of {', '.join(estimated)} a provisional estimate")
    return MixtureOptics(
        species=species,
        models=models,
        pairs=pairs,
        intervals=_intersection([m.intervals for m in models]),
        status="provisional" if reasons else "validated",
        status_reason="; ".join(reasons),
    )


def mixture_optics(species: Sequence[Gas | str]) -> MixtureOptics:
    """Return the optics of a mixture of gases (OEP 4).

    Parameters
    ----------
    species : sequence of Gas or str
        The mixture's species, each once.

    Returns
    -------
    MixtureOptics
        With each species' recommended dispersion model; one object per set of
        species, whatever their order.

    Raises
    ------
    OptmatError
        If the list is empty or a species repeats.

    Examples
    --------
    >>> optics = mixture_optics(["Ar", "He"])
    >>> optics.species, [p.status for p in optics.pairs]
    (('Ar', 'He'), ['curated', 'curated', 'estimate'])
    """
    names = [as_gas(s).value for s in species]
    if not names or len(set(names)) != len(names):
        msg = f"a mixture takes distinct species, not {names}"
        raise OptmatError(msg)
    return _optics(tuple(sorted(names)))

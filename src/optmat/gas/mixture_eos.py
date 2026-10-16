# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The density of a gas mixture at a pressure and temperature (OEP 4, section 7.4).

A mixture's equation of state is chosen per mixture, by the evidence for its pairs:

- **multi-fluid**: CoolProp's multi-fluid Helmholtz model of the mixture, where every
  pair has a binary model fitted to data of that pair. It reduces to each pure fluid's
  equation at :math:`x_s = 1`. It is ``qualified`` only where every pair passed the
  comparison with the measured cross second virials and mixture densities
  (``data/gas/mixture_pairs.toml``).
- **virial**: the fallback of section 4.2, each species' own equation at its partial
  density plus the curated cross second virial coefficients,

  .. math::

     p(\rho, T) = \sum_s p_s(x_s\rho, T) + RT\rho^2\sum_{s \ne t} x_s x_t B_{st}(T).

  It preserves the composition, :math:`N_s = x_s N` exactly, and holds while
  :math:`|\rho\sum_{s\ne t}x_s x_t B_{st}| \le 10^{-2}`. It is ``qualified`` where
  every pair's fallback passed the same comparison with the measured densities.
- **partial-density**: the same with every unknown cross coefficient zero: the
  composition-preserving form of the partial-pressure rule, ``provisional``, which
  warns when it is built.

The densities follow the CODATA convention of :mod:`optmat.gas.eos`: ``N = p / (Z k_B
T)`` with the equation's dimensionless ``Z``.
"""

import functools
import math
import tomllib
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike
from scipy import optimize

from optmat.conditions import State
from optmat.constants import K_B, N_A
from optmat.datasets import datasets
from optmat.exceptions import DataError, ModelRangeWarning, OptmatError, StateError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.composition import REDUCED, composition
from optmat.gas.eos import EquationOfState, _coolprop, equation_of_state
from optmat.gas.gases import as_gas
from optmat.identity import Identified

__all__ = [
    "METHODS",
    "PAIRS_FILE",
    "TRUNCATION_BOUND",
    "MixtureEquationOfState",
    "PairDensity",
    "mixture_equation_of_state",
]

PAIRS_FILE: Final = GAS_DATA_DIRECTORY / "mixture_pairs.toml"
"""The curated pairs of the mixture equations of state."""

METHODS: Final = ("multi-fluid", "virial", "partial-density")
"""The density methods, in the order the default is chosen."""

TRUNCATION_BOUND: Final = 1e-2
"""The bound on ``|rho sum_{s != t} x_s x_t B_st|`` within which the fallback holds."""

_R_CODATA: Final = N_A * K_B

_GAS_PHASES: Final = ("phase_gas", "phase_supercritical", "phase_supercritical_gas")

_CROSS_QUANTITIES: Final = (
    "cross_second_virial_coefficient",
    "excess_second_virial_coefficient",
)
"""What a pair's ``cross_virial`` may name: B_12, or the excess E it is built from."""


@dataclass(frozen=True)
class PairDensity:
    """One unlike pair of a mixture equation of state.

    Attributes
    ----------
    species, partner : str
        The pair, in sorted order.
    model : str
        CoolProp's reference key of the pair's binary model, or empty where it has none.
    kind : str
        ``fitted``, ``placeholder`` (every reducing parameter one and no departure
        function: not fitted to data of the pair) or ``none``.
    curated : bool
        Whether the pair is in ``mixture_pairs.toml``.
    qualified : bool
        Whether its multi-fluid model passed the comparison with the measured data.
    cross_virial : str
        The curated cross (or excess) second virial dataset, or empty.
    reason : str
        The comparison's result, or why the pair is what it is.
    fallback_qualified : bool
        Whether the fallback (method ``virial``) passed the same comparison with the
        measured mixture densities, within its truncation bound.
    fallback_reason : str
        That comparison's result.
    prefer_multi_fluid : bool
        Whether the author chose the multi-fluid model as this pair's default even
        though it is not qualified (dry air's pairs, TC-796).
    """

    species: str
    partner: str
    model: str
    kind: str
    curated: bool
    qualified: bool
    cross_virial: str
    reason: str
    fallback_qualified: bool = False
    fallback_reason: str = ""
    prefer_multi_fluid: bool = False


@functools.cache
def _pair_records() -> dict[tuple[str, str], Mapping[str, Any]]:
    """Return the curated pairs of ``mixture_pairs.toml``.

    Returns
    -------
    dict
        ``(species, partner)``, sorted, to the pair's table.

    Raises
    ------
    DataError
        For a pair not in sorted order, or a cross virial that is not a curated
        ``cross_second_virial_coefficient`` of the pair.
    """
    with PAIRS_FILE.open("rb") as handle:
        record = tomllib.load(handle)
    by_key = {d.key: d for d in datasets()}
    out = {}
    for name, entry in record["pairs"].items():
        species, partner = (as_gas(s).value for s in name.split("-"))
        if species >= partner:
            msg = f"mixture_pairs.toml [{name}]: name a pair in sorted order"
            raise DataError(msg)
        key = entry.get("cross_virial", "")
        if key:
            dataset = by_key.get(key)
            if (
                dataset is None
                or dataset.quantity not in _CROSS_QUANTITIES
                or {str(dataset.material), str(dataset.parameters.get("partner"))}
                != {species, partner}
            ):
                msg = (
                    f"mixture_pairs.toml [{name}]: {key!r} is not a cross virial of it"
                )
                raise DataError(msg)
        out[species, partner] = entry
    return out


def _cas(species: str) -> str:
    """Return a species' CAS number, CoolProp's key for its pairs.

    Parameters
    ----------
    species : str
        The species.

    Returns
    -------
    str
        The CAS number.
    """
    return str(
        _coolprop().get_fluid_param_string(equation_of_state(species).fluid, "CAS")
    )


def _binary(species: str, partner: str) -> tuple[str, str]:
    """Return CoolProp's binary model of a pair: its reference key and kind.

    Parameters
    ----------
    species, partner : str
        The pair.

    Returns
    -------
    tuple of str
        ``(model, kind)``: the reference key and ``fitted`` or ``placeholder``, or
        ``("", "none")`` where CoolProp has no model of the pair.
    """
    coolprop = _coolprop()
    first, second = _cas(species), _cas(partner)
    for a, b in ((first, second), (second, first)):
        try:
            model = str(coolprop.get_mixture_binary_pair_data(a, b, "BibTeX"))
        except ValueError:
            continue
        values = [
            float(coolprop.get_mixture_binary_pair_data(a, b, key))
            for key in ("betaT", "gammaT", "betaV", "gammaV", "F")
        ]
        unfitted = all(v == 1.0 for v in values[:4]) and values[4] == 0.0
        return model, "placeholder" if unfitted else "fitted"
    return "", "none"


def _pair(species: str, partner: str) -> PairDensity:
    """Return one pair, checked against CoolProp's model of it.

    Parameters
    ----------
    species, partner : str
        The pair, sorted.

    Returns
    -------
    PairDensity
        Its model, kind, curation and qualification.

    Raises
    ------
    DataError
        If a curated pair's model is not CoolProp's.
    """
    model, kind = _binary(species, partner)
    entry = _pair_records().get((species, partner))
    if entry is None:
        reason = {
            "fitted": "CoolProp's binary model is not curated or compared",
            "placeholder": "CoolProp's binary model is a placeholder, not fitted",
            "none": "CoolProp has no binary model of the pair",
        }[kind]
        return PairDensity(species, partner, model, kind, False, False, "", reason)
    if entry["model"] != model:
        msg = (
            f"mixture_pairs.toml [{species}-{partner}]: model {entry['model']!r}, "
            f"but CoolProp's is {model!r}"
        )
        raise DataError(msg)
    return PairDensity(
        species=species,
        partner=partner,
        model=model,
        kind=kind,
        curated=True,
        qualified=bool(entry["qualified"]) and kind == "fitted",
        cross_virial=str(entry.get("cross_virial", "")),
        reason=str(entry["reason"]),
        fallback_qualified=bool(entry.get("fallback_qualified", False)),
        fallback_reason=str(entry.get("fallback_reason", "")),
        prefer_multi_fluid=bool(entry.get("prefer_multi_fluid", False)),
    )


@functools.cache
def _cross_table(key: str) -> tuple[np.ndarray, np.ndarray, bool]:
    """Return a curated cross or excess virial's temperatures and values.

    Parameters
    ----------
    key : str
        The dataset, ``<source>/<id>``.

    Returns
    -------
    tuple
        Temperatures (K, ascending), the values (m^3/mol), and whether they are the
        excess E rather than B_12.
    """
    dataset = {d.key: d for d in datasets()}[key]
    temperatures = dataset.abscissa_values()
    order = np.argsort(temperatures)
    temperatures = temperatures[order]
    values = np.asarray(dataset.scaled_values(), dtype=float)[order] * 1e-6
    return temperatures, values, dataset.quantity == "excess_second_virial_coefficient"


@functools.cache
def _pure_second_virial(fluid: str, temperature: float) -> float:
    """Return a pure fluid's second virial coefficient from its reference equation.

    Parameters
    ----------
    fluid : str
        CoolProp's name of the fluid.
    temperature : float
        Temperature, K.

    Returns
    -------
    float
        The coefficient, m^3/mol.
    """
    backend = _coolprop().AbstractState("HEOS", fluid)
    backend.update(_coolprop().DmolarT_INPUTS, 1e-6, temperature)
    return float(backend.Bvirial())


@functools.cache
def _mixture_state(fluids: tuple[str, ...], fractions: tuple[float, ...]) -> Any:
    """Return CoolProp's multi-fluid state of a composition (one per composition).

    Parameters
    ----------
    fluids : tuple of str
        CoolProp's fluid names.
    fractions : tuple of float
        Their mole fractions.

    Returns
    -------
    CoolProp.AbstractState
        The ``HEOS`` mixture, its fractions set.
    """
    backend = _coolprop().AbstractState("HEOS", "&".join(fluids))
    backend.set_mole_fractions(list(fractions))
    return backend


@dataclass(frozen=True, eq=False)
class MixtureEquationOfState(Identified):
    """The equation of state of a gas mixture of fixed composition (OEP 4).

    Attributes
    ----------
    material : str
        The composition's name, or its species and fractions.
    composition : tuple of tuple
        ``(species, mole fraction)``, in the order given, summing to one.
    method : str
        One of :data:`METHODS`.
    qualified : bool
        Whether the method holds for this mixture by the evidence of section 7.4.
    pairs : tuple of PairDensity
        Every unlike pair.
    constituents : tuple of EquationOfState
        Each species' own equation of state, in the order of :attr:`composition`.
    backend : str
        CoolProp's version.
    temperature_range : tuple of float
        The temperatures every constituent's equation covers, K, and for the
        fallback every cross coefficient's measured temperatures.
    status : str
        ``validated`` where :attr:`qualified`, else ``provisional``.
    """

    IDENTITY_KIND: ClassVar[str] = "mixture_equation_of_state"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.eos",)

    material: str
    composition: tuple[tuple[str, float], ...]
    method: str
    qualified: bool
    pairs: tuple[PairDensity, ...]
    constituents: tuple[EquationOfState, ...]
    backend: str
    temperature_range: tuple[float, float]

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the equation's name, mixture and source records.

        Returns
        -------
        tuple
            ``("mixture_eos/<material>", material, sources)``.
        """
        found = {
            str(_pair_records()[p.species, p.partner]["source"])
            for p in self.pairs
            if p.curated
        } | {p.cross_virial.split("/", 1)[0] for p in self.pairs if p.cross_virial}
        return f"mixture_eos/{self.material}", self.material, tuple(sorted(found))

    @property
    def species(self) -> tuple[str, ...]:
        """The species, in the order of :attr:`composition`.

        Returns
        -------
        tuple of str
            The species.
        """
        return tuple(s for s, _ in self.composition)

    @property
    def molar_mass(self) -> float:
        """The composition's molar mass, kg/mol.

        The constituents' equation-of-state molar masses
        (:attr:`~optmat.gas.eos.EquationOfState.molar_mass`) weighted by their mole
        fractions, for every method. For dry air it is the mean of its composition,
        not the pseudo-pure air of
        :func:`~optmat.gas.eos.air_equation_of_state` (TC-832).

        Returns
        -------
        float
            The molar mass, kg/mol.

        Examples
        --------
        >>> mixture_equation_of_state({"Ar": 1.0}).molar_mass
        0.039948
        """
        return math.fsum(
            x * e.molar_mass
            for (_, x), e in zip(self.composition, self.constituents, strict=True)
        )

    def _fractions(self) -> np.ndarray:
        """Return the mole fractions.

        Returns
        -------
        numpy.ndarray
            In the order of :attr:`composition`.
        """
        return np.array([x for _, x in self.composition])

    def _check_temperature(self, temperature: float) -> None:
        """Refuse a temperature outside the constituents' equations.

        Parameters
        ----------
        temperature : float
            K.

        Raises
        ------
        StateError
            Outside :attr:`temperature_range`.
        """
        low, high = self.temperature_range
        if not low <= temperature <= high:
            msg = (
                f"no equation-of-state solution for {self.material} at {temperature} "
                f"K: outside its {low:g}-{high:g} K"
            )
            raise StateError(msg)

    # The multi-fluid model.

    def _backend(self) -> Any:
        """Return CoolProp's state of the mixture.

        Returns
        -------
        CoolProp.AbstractState
            The multi-fluid mixture.
        """
        fluids = tuple(e.fluid for e in self.constituents)
        return _mixture_state(fluids, tuple(float(x) for x in self._fractions()))

    def _multi_fluid_z(self, pressure: float, temperature: float) -> float:
        """Return the multi-fluid compressibility factor, checking the phase.

        Parameters
        ----------
        pressure : float
            Pa.
        temperature : float
            K.

        Returns
        -------
        float
            ``Z``.

        Raises
        ------
        StateError
            Where the mixture's flash fails or finds it condensed.
        """
        backend = self._backend()
        coolprop = _coolprop()
        try:
            backend.update(coolprop.PT_INPUTS, pressure, temperature)
        except ValueError as exc:
            msg = (
                f"no equation-of-state solution for {self.material} at {pressure} Pa "
                f"and {temperature} K: {exc}"
            )
            raise StateError(msg) from exc
        phase = backend.phase()
        if phase not in {coolprop.get_phase_index(p) for p in _GAS_PHASES}:
            msg = (
                f"no equation-of-state solution for {self.material} at {pressure} Pa "
                f"and {temperature} K: not a single gas phase"
            )
            raise StateError(msg)
        return float(backend.compressibility_factor())

    # The fallback.

    def _cross(self, temperature: float) -> dict[tuple[str, str], float]:
        """Return the cross coefficients the fallback takes, m^3/mol.

        Parameters
        ----------
        temperature : float
            Temperature, K.

        Returns
        -------
        dict
            ``(species, partner)`` to B_12; a pair without one is left out (zero).
        """
        out = {}
        fluids = {e.material: e.fluid for e in self.constituents}
        for p in self.pairs:
            if not p.cross_virial:
                continue
            # The fallback's temperature range is within every table's
            # (_check_temperature), so this interpolates and never extrapolates.
            temperatures, values, excess = _cross_table(p.cross_virial)
            value = float(np.interp(temperature, temperatures, values))
            if excess:
                # B_12 = E + (B_11 + B_22)/2, with the reference equations' pure B's.
                value += 0.5 * sum(
                    _pure_second_virial(fluids[g], temperature)
                    for g in (p.species, p.partner)
                )
            out[p.species, p.partner] = value
        return out

    def _cross_sum(self, temperature: float) -> float:
        """Return ``sum_{s != t} x_s x_t B_st``, m^3/mol.

        Parameters
        ----------
        temperature : float
            K.

        Returns
        -------
        float
            Over ordered pairs, so twice the sum over unordered ones.
        """
        x = dict(self.composition)
        return sum(
            2 * x[s] * x[t] * b for (s, t), b in self._cross(temperature).items()
        )

    def _partial_pressure(
        self, eos: EquationOfState, density: float, temperature: float
    ) -> float:
        """Return a constituent's own pressure at its partial density.

        Parameters
        ----------
        eos : EquationOfState
            The constituent's equation.
        density : float
            Its number density, m^-3.
        temperature : float
            K.

        Returns
        -------
        float
            Pa, in the CODATA convention.

        Raises
        ------
        StateError
            Where the partial density is one the pure gas reaches only condensed.
        """
        _, top = eos.phase_domain(temperature)
        # Below the critical temperature the gas ends at its saturated vapour (the
        # upper end excluded); above it, the fluid does not condense.
        below = temperature < eos.critical_temperature
        if below and density >= eos._density_at(top * (1 - 1e-12), temperature):  # pyright: ignore[reportPrivateUsage]
            msg = (
                f"no equation-of-state solution for {self.material}: {eos.material} at "
                f"its partial density {density:.6g} m^-3 and {temperature} K is beyond "
                "its gas phase"
            )
            raise StateError(msg)
        return eos._start(density, temperature)  # pyright: ignore[reportPrivateUsage]

    def _fallback_pressure(self, density: float, temperature: float) -> float:
        """Return the fallback's pressure at a total number density.

        Parameters
        ----------
        density : float
            Total number density, m^-3, > 0.
        temperature : float
            K.

        Returns
        -------
        float
            The pressure, Pa.
        """
        pure = sum(
            self._partial_pressure(eos, x * density, temperature)
            for eos, (_, x) in zip(self.constituents, self.composition, strict=True)
        )
        rho = density / N_A
        return pure + _R_CODATA * temperature * rho**2 * self._cross_sum(temperature)

    def _truncation(self, density: float, temperature: float) -> None:
        """Warn where the fallback is beyond its truncation bound.

        Parameters
        ----------
        density : float
            Total number density, m^-3.
        temperature : float
            K.
        """
        if self.method == "multi-fluid":
            return
        term = abs(density / N_A * self._cross_sum(temperature))
        if term > TRUNCATION_BOUND:
            msg = (
                f"{self.material}: the fallback's cross term is {term:.3g}, beyond its "
                f"truncation bound {TRUNCATION_BOUND:g} (OEP 4, section 4.2)"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=4)

    # The interface.

    def _number_density(self, pressure: float, temperature: float) -> float:
        """Return the total number density at one state.

        Parameters
        ----------
        pressure : float
            Pa, > 0.
        temperature : float
            K.

        Returns
        -------
        float
            The total number density, m^-3.
        """
        self._check_temperature(temperature)
        if self.method == "multi-fluid":
            z = self._multi_fluid_z(pressure, temperature)
            return pressure / (z * K_B * temperature)
        ideal = pressure / (K_B * temperature)
        low, high = ideal / 2, ideal * 2
        # p(N) rises on the gas phase: widen the bracket until it holds the root.
        while self._fallback_pressure(low, temperature) > pressure:
            low /= 2
        while self._fallback_pressure(high, temperature) < pressure:
            high *= 2
        root = optimize.brentq(
            lambda n: self._fallback_pressure(n, temperature) - pressure,
            low,
            high,
            xtol=np.finfo(float).tiny,
            rtol=4 * np.finfo(float).eps,
        )
        return float(root)

    def number_density(self, state: State) -> float:
        """Return the mixture's total number density at a state.

        Parameters
        ----------
        state : State
            Pressure and temperature.

        Returns
        -------
        float
            The total number density, m^-3.

        Raises
        ------
        StateError
            Outside the constituents' temperatures, or where the mixture (or, for the
            fallback, a constituent at its partial density) is not a gas.

        Warns
        -----
        ModelRangeWarning
            Where the fallback is beyond its truncation bound.

        Examples
        --------
        >>> from optmat.conditions import State
        >>> eos = mixture_equation_of_state({"He": 0.5, "Ne": 0.5})
        >>> eos.method
        'virial'
        >>> f"{eos.number_density(State(pressure=1e6, temperature=294.0)):.6e}"
        '2.451570e+26'
        """
        density = self._number_density(state.pressure, state.temperature)
        self._truncation(density, state.temperature)
        return density

    def partial_densities(self, state: State) -> dict[str, float]:
        """Return each species' number density, ``x_s N``.

        Parameters
        ----------
        state : State
            Pressure and temperature.

        Returns
        -------
        dict
            Species to m^-3, the keys :class:`~optmat.gas.mixture.MixtureOptics`
            takes.
        """
        total = self.number_density(state)
        return {s: x * total for s, x in self.composition}

    def compressibility(self, state: State) -> float:
        """Return the mixture's compressibility factor ``Z = p / (N k_B T)``.

        Parameters
        ----------
        state : State
            Pressure and temperature.

        Returns
        -------
        float
            Dimensionless.
        """
        return state.pressure / (self.number_density(state) * K_B * state.temperature)

    def _pressure(self, density: float, temperature: float) -> float:
        """Return the pressure at one total number density and temperature.

        Parameters
        ----------
        density : float
            Total number density, m^-3, > 0.
        temperature : float
            K.

        Returns
        -------
        float
            The pressure, Pa.

        Raises
        ------
        StateError
            At a non-positive density, or where the mixture is not a gas.
        """
        if not density > 0:
            msg = f"the pressure of {self.material} needs a density > 0, not {density}"
            raise StateError(msg)
        self._check_temperature(temperature)
        if self.method != "multi-fluid":
            self._truncation(density, temperature)
            return self._fallback_pressure(density, temperature)
        backend = self._backend()
        molar = density / N_A * _R_CODATA / float(backend.gas_constant())
        backend.update(_coolprop().DmolarT_INPUTS, molar, temperature)
        start = float(backend.p())

        def forward(p: float) -> float:
            """Return the forward law's density minus the target.

            Parameters
            ----------
            p : float
                Pressure, Pa.

            Returns
            -------
            float
                The difference, m^-3.
            """
            z = self._multi_fluid_z(p, temperature)
            return p / (z * K_B * temperature) - density

        # CoolProp's (rho, T) pressure is within the equation's iteration tolerance
        # of the root, far inside 1e-6.
        return float(
            optimize.brentq(
                forward,
                start * (1 - 1e-6),
                start * (1 + 1e-6),
                xtol=np.finfo(float).tiny,
                rtol=4 * np.finfo(float).eps,
            )
        )

    def pressure(
        self, number_density: ArrayLike | float, temperature: ArrayLike | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the pressure at total number densities and temperatures.

        Parameters
        ----------
        number_density : array_like or float
            Total number density, m^-3, > 0.
        temperature : array_like or float
            K; broadcast against ``number_density``.

        Returns
        -------
        numpy.ndarray
            Pa.

        Raises
        ------
        StateError
            At a non-positive density, or where the mixture is not a gas.

        Examples
        --------
        >>> eos = mixture_equation_of_state({"He": 0.5, "Ne": 0.5})
        >>> f"{float(eos.pressure(2.451570e26, 294.0)) / 1e6:.6f}"
        '1.000000'
        """
        densities, temperatures = np.broadcast_arrays(
            np.asarray(number_density, dtype=float),
            np.asarray(temperature, dtype=float),
        )
        out = np.array(
            [
                self._pressure(float(n), float(t))
                for n, t in zip(densities.ravel(), temperatures.ravel(), strict=True)
            ]
        )
        return out.reshape(densities.shape)

    def second_virial(self, temperature: float) -> float:
        """Return the mixture's second virial coefficient, the method's own.

        Parameters
        ----------
        temperature : float
            K.

        Returns
        -------
        float
            ``B_mix``, m^3/mol: the multi-fluid model's at zero density, or the
            fallback's ``sum_s x_s^2 B_ss + sum_{s != t} x_s x_t B_st``.
        """
        self._check_temperature(temperature)
        if self.method == "multi-fluid":
            backend = self._backend()
            backend.update(_coolprop().DmolarT_INPUTS, 1e-6, temperature)
            return float(backend.Bvirial())
        like = sum(
            x**2 * _pure_second_virial(eos.fluid, temperature)
            for eos, (_, x) in zip(self.constituents, self.composition, strict=True)
        )
        return like + self._cross_sum(temperature)


def _named(name: str) -> tuple[tuple[str, float], ...]:
    """Return a named composition's species and fractions.

    Parameters
    ----------
    name : str
        A reduced composition (:data:`~optmat.gas.composition.REDUCED`).

    Returns
    -------
    tuple of tuple
        ``(species, fraction)``, renormalised.

    Raises
    ------
    OptmatError
        For another name.
    """
    if name not in REDUCED:
        msg = (
            f"no mixture composition named {name!r}; the named ones are "
            f"{tuple(REDUCED)}"
        )
        raise OptmatError(msg)
    return tuple((c.material, c.mole_fraction) for c in composition(name))


def _composition(
    spec: Mapping[Any, float] | str,
) -> tuple[str, tuple[tuple[str, float], ...]]:
    """Return a composition's name and its checked species and fractions.

    Parameters
    ----------
    spec : Mapping or str
        Species to mole fraction, or a named composition.

    Returns
    -------
    tuple
        The name and ``(species, fraction)`` pairs.

    Raises
    ------
    OptmatError
        For a repeated or unknown species, a fraction not in (0, 1], or fractions
        that do not sum to one.
    """
    if isinstance(spec, str):
        return spec, _named(spec)
    pairs = tuple((as_gas(k).value, float(v)) for k, v in spec.items())
    names = [s for s, _ in pairs]
    if not pairs or len(set(names)) != len(names):
        msg = f"a mixture takes distinct species, not {names}"
        raise OptmatError(msg)
    if any(not 0 < x <= 1 for _, x in pairs):
        msg = f"mole fractions must be in (0, 1]: {dict(pairs)}"
        raise OptmatError(msg)
    if not math.isclose(sum(x for _, x in pairs), 1.0, rel_tol=1e-9):
        msg = f"mole fractions sum to {sum(x for _, x in pairs)}, not one"
        raise OptmatError(msg)
    name = " + ".join(f"{s} {x:g}" for s, x in pairs)
    return name, pairs


def _choose(pairs: tuple[PairDensity, ...], method: str | None) -> tuple[str, bool]:
    """Return the method and whether it is qualified (section 7.4).

    Parameters
    ----------
    pairs : tuple of PairDensity
        The mixture's unlike pairs.
    method : str or None
        The method asked for, or None for the default.

    Returns
    -------
    tuple
        The method and its qualification.

    Raises
    ------
    OptmatError
        For an unknown method, or one the pairs cannot support.
    """
    fitted = all(p.kind == "fitted" for p in pairs)
    crossed = all(p.cross_virial for p in pairs)
    if method is None:
        if all(p.qualified for p in pairs):
            return "multi-fluid", True
        if pairs and fitted and all(p.prefer_multi_fluid for p in pairs):
            # The author's choice over the fallback (TC-796): unqualified, provisional.
            return "multi-fluid", False
        if crossed:
            return "virial", all(p.fallback_qualified for p in pairs)
        return "partial-density", False
    if method not in METHODS:
        msg = f"unknown density method {method!r}; the methods are {METHODS}"
        raise OptmatError(msg)
    if method == "multi-fluid" and not fitted:
        missing = [
            f"{p.species}-{p.partner} ({p.kind})" for p in pairs if p.kind != "fitted"
        ]
        msg = f"no fitted multi-fluid model for {', '.join(missing)}"
        raise OptmatError(msg)
    if method == "virial" and not crossed:
        missing = [f"{p.species}-{p.partner}" for p in pairs if not p.cross_virial]
        msg = f"no curated cross virial for {', '.join(missing)}"
        raise OptmatError(msg)
    qualified = {
        "multi-fluid": all(p.qualified for p in pairs),
        "virial": all(p.fallback_qualified for p in pairs),
        "partial-density": not pairs,
    }[method]
    return method, qualified


def mixture_equation_of_state(
    composition: Mapping[Any, float] | str, *, method: str | None = None
) -> MixtureEquationOfState:
    """Return the equation of state of a gas mixture (OEP 4, sections 4.2 and 7.4).

    Parameters
    ----------
    composition : Mapping or str
        Species to mole fraction, summing to one, or a named composition
        (``"dry air (N2, O2, Ar)"``).
    method : str, optional
        One of :data:`METHODS`; by default the first the mixture's pairs qualify for,
        else ``partial-density``.

    Returns
    -------
    MixtureEquationOfState
        The equation, with its method, qualification and pairs.

    Raises
    ------
    OptmatError
        For an invalid composition or a method the pairs cannot support.
    StateError
        For a species without an equation of state (ozone).

    Warns
    -----
    ModelRangeWarning
        When the equation is provisional.

    Examples
    --------
    >>> eos = mixture_equation_of_state({"Ar": 1.0})
    >>> eos.method, eos.qualified
    ('multi-fluid', True)
    """
    name, pairs = _composition(composition)
    constituents = tuple(equation_of_state(s) for s, _ in pairs)
    species = sorted(s for s, _ in pairs)
    unlike = tuple(_pair(s, t) for i, s in enumerate(species) for t in species[i + 1 :])
    chosen, qualified = _choose(unlike, method)
    low = max(e.temperature_range[0] for e in constituents)
    high = min(e.temperature_range[1] for e in constituents)
    if chosen != "multi-fluid":
        # The fallback holds where its cross coefficients are measured.
        for p in unlike:
            if p.cross_virial:
                temperatures = _cross_table(p.cross_virial)[0]
                low = max(low, float(temperatures[0]))
                high = min(high, float(temperatures[-1]))
    unqualified = [
        f"{p.species}-{p.partner}"
        for p in unlike
        if not (p.fallback_qualified if chosen == "virial" else p.qualified)
    ]
    eos = MixtureEquationOfState(
        material=name,
        composition=pairs,
        method=chosen,
        qualified=qualified,
        pairs=unlike,
        constituents=constituents,
        backend=constituents[0].backend,
        temperature_range=(low, high),
        status="validated" if qualified else "provisional",
        status_reason=(
            ""
            if qualified
            else f"{chosen}: pairs not qualified: {', '.join(unqualified)}"
        ),
    )
    if not qualified:
        msg = f"{name}: a provisional mixture equation of state ({eos.status_reason})"
        warnings.warn(msg, ModelRangeWarning, stacklevel=2)
    return eos

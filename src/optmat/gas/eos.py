# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Number densities of the gases: the ideal gas, the real gas, and its inverse.

A refractivity or an n2 measured at a stated pressure and temperature is reduced either
to the *ideal-gas* density at that state or to the *real-gas* density, depending on the
paper. To compare measurements, the package turns each one into a per-particle
quantity, which needs the number density the author actually meant.

Real-gas densities come from CoolProp's Helmholtz-energy equations of state, imported
lazily because CoolProp takes about a second to load. :func:`equation_of_state` returns
a gas's equation as an object with an identity: its forward law
``N = p / (Z k_B T)``, its inverse (the pressure at a number density and temperature),
its speed of sound, its molar mass and its phase domain (OEP 1, section 4.2).

**The phase domain.** optmat is a library of gases. The forward law and its inverse
share one domain: the gas below saturation (``T < T_c``, ``p < p_sigma(T)``) and the
supercritical fluid (``T >= T_c``, ``p <= p_max``, the equation's stated upper
pressure), within the equation's temperatures above its triple point. Wherever the
fluid has a melting line in CoolProp, the domain also ends just below the melting
pressure ``p_melt(T)``, which binds above the critical temperature for several gases
(nitrogen to 284 K, argon to 254 K; TC-811). A state outside it, a liquid, a solid or
a two-phase state among them, raises :class:`~optmat.exceptions.StateError`.

**The CODATA convention.** The density is formed from ``Z`` and CODATA's ``k_B`` rather
than from the equation of state's mass density, because each CoolProp equation carries
the gas constant current when it was published: argon's (Tegeler et al. 1999) and
nitrogen's use R = 8.31451 J/(mol K), which puts their mass densities 5.7e-6 away from
the CODATA ideal gas at vanishing pressure. ``Z`` is dimensionless and free of that
offset.
"""

import atexit
import functools
from dataclasses import dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat.conditions import State
from optmat.constants import K_B, N_A
from optmat.exceptions import StateError
from optmat.gas.gases import Gas, as_gas, info
from optmat.identity import Identified

__all__ = [
    "EVALUATION_VERSION",
    "ROUND_TRIP_TOLERANCE",
    "SATURATION_MARGIN",
    "EquationOfState",
    "air_equation_of_state",
    "compressibility",
    "equation_of_state",
    "ideal_number_density",
    "number_density",
    "pressure",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

ROUND_TRIP_TOLERANCE: Final = 1e-12
"""The bound on ``|N(p(N, T), T) / N - 1|`` that the inverse guarantees."""

SATURATION_MARGIN: Final = 1e-5
"""How far below the vapour pressure the gas's domain ends, relative.

CoolProp's pressure-temperature solve refuses pressures within 1e-6 of the vapour
pressure, so the domain stops just short of it."""

_MELTING_MARGIN: Final = 1e-5
"""How far below the melting pressure the fluid's domain ends, relative.

CoolProp refuses a state with ``T < T_melt(p)``, its melting temperature solved from
the melting line; the domain stops short of the line so that its end is clear of that
solve (TC-811)."""

_R_CODATA: Final = N_A * K_B


def ideal_number_density(state: State) -> float:
    """Return the ideal-gas number density ``p / (k_B T)``.

    Parameters
    ----------
    state : State
        Pressure and temperature.

    Returns
    -------
    float
        Number density in m^-3.

    Examples
    --------
    >>> from optmat.conditions import ONE_ATM_0C
    >>> f"{ideal_number_density(ONE_ATM_0C):.6e}"
    '2.686780e+25'
    """
    return state.pressure / (K_B * state.temperature)


@functools.cache
def _coolprop() -> Any:
    """Return CoolProp's module, imported on first use.

    Returns
    -------
    module
        ``CoolProp.CoolProp``.

    Raises
    ------
    ImportError
        Without the ``gas`` extra.
    """
    try:
        import CoolProp.CoolProp as coolprop  # noqa: N813, PLC0415 - import cost
    except ImportError as exc:
        msg = "real-gas densities need CoolProp: install optmat[gas]"
        raise ImportError(msg) from exc
    return coolprop


@functools.cache
def _abstract_state(fluid: str) -> Any:
    """Return CoolProp's Helmholtz-energy state of a fluid (one per fluid, reused).

    Parameters
    ----------
    fluid : str
        CoolProp's name of the fluid.

    Returns
    -------
    CoolProp.AbstractState
        The ``HEOS`` backend's state.
    """
    return _coolprop().AbstractState("HEOS", fluid)


def _release() -> None:
    """Release the cached CoolProp states before the interpreter's shutdown.

    CoolProp's bindings report the objects still alive when the interpreter exits as
    leaks; the cache holds them, so it is emptied first.
    """
    _compressibility.cache_clear()
    _abstract_state.cache_clear()


atexit.register(_release)


@functools.cache
def _compressibility(fluid: str, pressure: float, temperature: float) -> float:
    """Return the compressibility factor at a state, unchecked (cached).

    Parameters
    ----------
    fluid : str
        CoolProp's name of the fluid.
    pressure : float
        Pa.
    temperature : float
        K.

    Returns
    -------
    float
        ``Z = p / (rho R T)`` in the equation's own units.
    """
    backend = _abstract_state(fluid)
    try:
        backend.update(_coolprop().PT_INPUTS, pressure, temperature)
    except ValueError as exc:
        msg = (
            f"no equation-of-state solution for {fluid} at {pressure} Pa and "
            f"{temperature} K: {exc}"
        )
        raise StateError(msg) from exc
    return float(backend.compressibility_factor())


@dataclass(frozen=True)
class EquationOfState(Identified):
    """The real-gas equation of state of a gas, with its inverse and phase domain.

    Attributes
    ----------
    material : str
        The gas (``"Ar"``) or mixture (``"dry air (Lemmon 2000)"``) of the equation.
    fluid : str
        CoolProp's name of the fluid.
    reference : str
        CoolProp's reference key of the equation (``Tegeler-JPCRD-1999`` for argon).
    backend : str
        The backend and its version, such as ``CoolProp 8.0.0``.
    gas_constant : float
        The equation's own gas constant, J/(mol K).
    critical_temperature : float
        K.
    maximum_pressure : float
        The equation's stated upper pressure, Pa.
    temperature_range : tuple of float
        The equation's temperatures, K.
    status : str
        ``validated``.
    """

    IDENTITY_KIND: ClassVar[str] = "equation_of_state"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.eos",)

    material: str
    fluid: str
    reference: str
    backend: str
    gas_constant: float
    critical_temperature: float
    maximum_pressure: float
    temperature_range: tuple[float, float]

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the equation's name, gas and source records.

        Returns
        -------
        tuple
            ``("eos/<material>", material, ())``: the equation is named by its
            reference key.
        """
        return f"eos/{self.material}", self.material, ()

    @property
    def molar_mass(self) -> float:
        """The molar mass the backend's equation uses, kg/mol.

        It is the mass of a mole in the backend's mass-based quantities, its mass
        density and :meth:`speed_of_sound`, and equals the gas's
        :attr:`~optmat.gas.transport.TransportModel.molar_mass`. The number densities
        do not use it (the CODATA convention). It is resolved from the backend, which
        :attr:`fluid` and :attr:`backend` fix, so it is not a field and is not in the
        identity. For :func:`air_equation_of_state` it is CoolProp's pseudo-pure air,
        0.02896546 kg/mol, not the 28.9586 g/mol Lemmon et al. (2000) print, the mean
        of their composition (TC-832).

        Returns
        -------
        float
            The molar mass, kg/mol.

        Examples
        --------
        >>> equation_of_state("Ar").molar_mass
        0.039948
        >>> air_equation_of_state().molar_mass
        0.02896546
        """
        return float(_abstract_state(self.fluid).molar_mass())

    def saturation_pressure(self, temperature: float) -> float | None:
        """Return the vapour pressure at a temperature below the critical one.

        Parameters
        ----------
        temperature : float
            Temperature, K.

        Returns
        -------
        float or None
            ``p_sigma(T)``, Pa; None at or above the critical temperature.
        """
        if temperature >= self.critical_temperature:
            return None
        backend = _abstract_state(self.fluid)
        backend.update(_coolprop().QT_INPUTS, 1.0, temperature)
        return float(backend.p())

    def _melting_pressure(self, temperature: float) -> float | None:
        """Return CoolProp's melting pressure at a temperature, where its line has one.

        Parameters
        ----------
        temperature : float
            K.

        Returns
        -------
        float or None
            ``p_melt(T)``, Pa; None for a fluid without a melting line in CoolProp
            (xenon, SF6, N2O) and outside the line's temperatures. Above them the
            melting pressure exceeds the line's top, which is at or above the
            equation's pressure for every fluid; below them (deuterium, 18.7-19.7 K)
            the vapour pressure binds (TC-811).
        """
        backend = _abstract_state(self.fluid)
        if not backend.has_melting_line():
            return None
        coolprop = _coolprop()
        low = backend.melting_line(coolprop.iT_min, -1, 0.0)
        high = backend.melting_line(coolprop.iT_max, -1, 0.0)
        if not low <= temperature <= high:
            return None
        return float(backend.melting_line(coolprop.iP, coolprop.iT, temperature))

    def _domain_top(self, temperature: float) -> tuple[float, str]:
        """Return the upper end of the gas's domain at a temperature, and its bound.

        Parameters
        ----------
        temperature : float
            K, within :attr:`temperature_range` and above its lower end.

        Returns
        -------
        tuple
            The pressure, Pa, and the bound that sets it: ``"vapour pressure"``,
            ``"melting pressure"`` or ``"equation's pressure"``.

        Raises
        ------
        StateError
            Outside the equation's temperatures, or at its lower end.
        """
        low, high = self.temperature_range
        if temperature == low:
            # CoolProp refuses the gas at the triple point itself, below the triple
            # pressure (TC-811).
            msg = (
                f"no equation-of-state solution for {self.material} at "
                f"{temperature} K: at its triple point, where CoolProp refuses the gas"
            )
            raise StateError(msg)
        if not low < temperature <= high:
            msg = (
                f"no equation-of-state solution for {self.material} at "
                f"{temperature} K: outside its {low:g}-{high:g} K"
            )
            raise StateError(msg)
        saturation = self.saturation_pressure(temperature)
        if saturation is None:
            top, bound = self.maximum_pressure, "equation's pressure"
        else:
            top, bound = saturation * (1 - SATURATION_MARGIN), "vapour pressure"
        melting = self._melting_pressure(temperature)
        if melting is not None and melting * (1 - _MELTING_MARGIN) < top:
            top, bound = melting * (1 - _MELTING_MARGIN), "melting pressure"
        return top, bound

    def phase_domain(self, temperature: float) -> tuple[float, float]:
        """Return the pressures of the gas's domain at a temperature.

        Parameters
        ----------
        temperature : float
            K, within :attr:`temperature_range`, above its lower end (the triple
            point).

        Returns
        -------
        tuple of float
            ``(0, p_sigma(T) (1 - SATURATION_MARGIN))`` below the critical temperature,
            its upper end excluded, and ``(0, p_max)`` at or above it, Pa; where
            CoolProp's melting pressure, less 1e-5 of it, is lower, the upper end is
            that instead (nitrogen up to 284 K, argon up to 254 K; TC-811).

        Raises
        ------
        StateError
            Outside the equation's temperatures, or at the triple point.

        Examples
        --------
        >>> low, high = equation_of_state("CO2").phase_domain(280.0)
        >>> f"{high / 1e5:.2f} bar"
        '41.61 bar'
        >>> nitrogen = equation_of_state("N2")
        >>> bool(nitrogen.phase_domain(280.0)[1] < nitrogen.maximum_pressure)
        True
        """
        return 0.0, self._domain_top(temperature)[0]

    def _check(self, state: State) -> None:
        """Refuse a state outside the phase domain.

        Parameters
        ----------
        state : State
            Pressure and temperature.

        Raises
        ------
        StateError
            If the state is condensed, solid, above the equation's pressure, or outside
            its temperatures.
        """
        top, bound = self._domain_top(state.temperature)
        if bound == "vapour pressure" and state.pressure >= top:
            msg = (
                f"no equation-of-state solution for {self.material} at {state}: "
                f"condensed, at or near its vapour pressure (the domain ends at "
                f"{top:.6g} Pa)"
            )
            raise StateError(msg)
        if bound == "melting pressure" and state.pressure > top:
            msg = (
                f"no equation-of-state solution for {self.material} at {state}: "
                f"solid, at or above its melting pressure (the domain ends at "
                f"{top:.6g} Pa)"
            )
            raise StateError(msg)
        if state.pressure > top:
            msg = (
                f"no equation-of-state solution for {self.material} at {state}: "
                f"above its equation's pressure, {top:.6g} Pa"
            )
            raise StateError(msg)

    def compressibility(self, state: State) -> float:
        """Return the compressibility factor ``Z = p / (N k_B T)``.

        Parameters
        ----------
        state : State
            Pressure and temperature, in the phase domain.

        Returns
        -------
        float
            Dimensionless; 1 for an ideal gas.
        """
        self._check(state)
        return _compressibility(self.fluid, state.pressure, state.temperature)

    def number_density(self, state: State) -> float:
        """Return the real-gas number density ``p / (Z k_B T)``, the forward law.

        Parameters
        ----------
        state : State
            Pressure and temperature, in the phase domain.

        Returns
        -------
        float
            The number density, m^-3.

        Examples
        --------
        >>> from optmat.conditions import ONE_ATM_0C
        >>> f"{equation_of_state('He').number_density(ONE_ATM_0C):.6e}"
        '2.685349e+25'
        """
        return ideal_number_density(state) / self.compressibility(state)

    def speed_of_sound(self, state: State) -> float:
        """Return the real gas's speed of sound.

        Parameters
        ----------
        state : State
            Pressure and temperature, in the phase domain.

        Returns
        -------
        float
            The speed of sound, m/s.
        """
        self._check(state)
        backend = _abstract_state(self.fluid)
        backend.update(_coolprop().PT_INPUTS, state.pressure, state.temperature)
        return float(backend.speed_sound())

    def _density_at(self, pressure_value: float, temperature: float) -> float:
        """Return the forward law's number density, unchecked.

        Parameters
        ----------
        pressure_value : float
            Pressure, Pa.
        temperature : float
            Temperature, K.

        Returns
        -------
        float
            The number density, m^-3.
        """
        z = _compressibility(self.fluid, pressure_value, temperature)
        return pressure_value / (z * K_B * temperature)

    def _start(self, density: float, temperature: float) -> float:
        """Return the backend's own pressure at a number density, the solve's start.

        CoolProp's (density, temperature) solve, in the equation's molar density and
        its own gas constant, is within the equation's iteration tolerance of the
        root.

        Parameters
        ----------
        density : float
            Number density, m^-3.
        temperature : float
            Temperature, K.

        Returns
        -------
        float
            The pressure, Pa.
        """
        backend = _abstract_state(self.fluid)
        molar = density / N_A * _R_CODATA / self.gas_constant
        backend.update(_coolprop().DmolarT_INPUTS, molar, temperature)
        return float(backend.p())

    def _pressure(self, density: float, temperature: float) -> float:
        """Return the pressure at one number density and temperature.

        Parameters
        ----------
        density : float
            Number density, m^-3, > 0.
        temperature : float
            Temperature, K.

        Returns
        -------
        float
            The pressure, Pa.

        Raises
        ------
        StateError
            For a non-positive density, a density at or above the saturated vapour's
            (condensed), the melting fluid's (solid) or the density at the equation's
            pressure, or a solve that misses the round-trip tolerance.
        """
        if not density > 0:
            msg = f"the pressure of {self.material} needs a density > 0, not {density}"
            raise StateError(msg)
        top, bound = self._domain_top(temperature)
        # The upper end is excluded at the vapour pressure.
        top = top * (1 - 1e-12) if bound == "vapour pressure" else top
        if density >= self._density_at(top, temperature):
            msg = (
                f"no equation-of-state solution for {self.material} at "
                f"{density:.6g} m^-3 and {temperature} K: at or above the density at "
                f"the top of its phase domain, {top:.6g} Pa (its {bound})"
            )
            raise StateError(msg)
        # N(p) rises strictly on the domain: halving the lower end and doubling the
        # upper (up to the domain's top, where N(top) > density) brackets the root.
        start = min(self._start(density, temperature), top)
        low, high = start * (1 - 1e-9), min(start * (1 + 1e-9), top)
        while not (
            self._density_at(low, temperature)
            <= density
            <= self._density_at(high, temperature)
        ):
            low, high = low / 2, min(high * 2, top)
        from scipy import optimize  # noqa: PLC0415 - loaded on the first solve

        root = optimize.brentq(
            lambda p: self._density_at(p, temperature) - density,
            low,
            high,
            xtol=np.finfo(float).tiny,
            rtol=4 * np.finfo(float).eps,
        )
        error = abs(self._density_at(root, temperature) / density - 1)
        if not error <= ROUND_TRIP_TOLERANCE:
            msg = (
                f"the pressure of {self.material} at {density:.6g} m^-3 and "
                f"{temperature} K missed the round trip ({error:.2e})"
            )
            raise StateError(msg)
        return float(root)

    def pressure(
        self, number_density: ArrayLike | float, temperature: ArrayLike | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the pressure at number densities and temperatures, the inverse law.

        It solves ``p = N k_B T Z(p, T)`` on the phase domain, where ``N(p)`` rises
        strictly, so the root is unique, to :data:`ROUND_TRIP_TOLERANCE` in
        ``N(p(N, T), T) / N`` (OEP 1, section 4.2). It keeps the CODATA convention: its
        ``N`` is :meth:`number_density`'s, not the backend's mass density.

        Parameters
        ----------
        number_density : array_like or float
            Number density, m^-3, > 0 (the vacuum limit is refused).
        temperature : array_like or float
            K; broadcast against ``number_density``.

        Returns
        -------
        numpy.ndarray
            Pa, the broadcast shape.

        Raises
        ------
        StateError
            At a non-positive density, a density the gas reaches only condensed, solid
            or above its equation's pressure, or a temperature outside its range or at
            its triple point.

        Examples
        --------
        >>> from optmat.conditions import ONE_ATM_0C
        >>> eos = equation_of_state("Ar")
        >>> n = eos.number_density(ONE_ATM_0C)
        >>> f"{float(eos.pressure(n, 273.15)):.4f}"
        '101325.0000'
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


@functools.cache
def equation_of_state(gas: Gas | str) -> EquationOfState:
    """Return a gas's real-gas equation of state.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    EquationOfState
        CoolProp's Helmholtz-energy equation of the gas, with its identity.

    Raises
    ------
    StateError
        For a species without an equation of state (ozone).

    Examples
    --------
    >>> equation_of_state("Ar").reference
    'Tegeler-JPCRD-1999'
    """
    species = as_gas(gas)
    fluid = info(species).fluid
    if fluid is None:
        msg = f"CoolProp has no equation of state for {species.value}"
        raise StateError(msg)
    return _equation(species.value, fluid)


def _equation(material: str, fluid: str) -> EquationOfState:
    """Return CoolProp's equation of a fluid as an equation-of-state object.

    Parameters
    ----------
    material : str
        What the equation describes.
    fluid : str
        CoolProp's name of the fluid.

    Returns
    -------
    EquationOfState
        The equation, with the backend's version, reference key and limits.
    """
    coolprop = _coolprop()
    backend = _abstract_state(fluid)
    return EquationOfState(
        material=material,
        fluid=fluid,
        reference=str(coolprop.get_fluid_param_string(fluid, "BibTeX-EOS")),
        backend=f"CoolProp {coolprop.get_global_param_string('version')}",
        gas_constant=float(backend.gas_constant()),
        critical_temperature=float(backend.T_critical()),
        maximum_pressure=float(backend.pmax()),
        temperature_range=(float(backend.Tmin()), float(backend.Tmax())),
        status="validated",
        status_reason="CoolProp's reference equation of the fluid",
    )


@functools.cache
def air_equation_of_state() -> EquationOfState:
    """Return CoolProp's pseudo-pure dry air, the density of OEP 1's air models.

    Lemmon et al. (2000) fit dry air as one fluid of the composition N2 0.7812,
    O2 0.2096, Ar 0.0092 (``lemmon_2000``, TC-741), with the unlike-pair interactions
    in it: 0.1 % in density (k = 2) up to 873 K and 70 MPa.

    Returns
    -------
    EquationOfState
        The equation, for the material ``"dry air (Lemmon 2000)"``.

    Examples
    --------
    >>> air_equation_of_state().reference
    'Lemmon-JPCRD-2000'
    """
    return _equation("dry air (Lemmon 2000)", "Air")


def number_density(gas: Gas | str, state: State, *, ideal: bool = False) -> float:
    """Return the number density of a pure gas.

    Parameters
    ----------
    gas : Gas or str
        The species.
    state : State
        Pressure and temperature.
    ideal : bool, optional
        If True, return the ideal-gas density ``p / (k_B T)``; otherwise (the default)
        the real-gas density ``p / (Z k_B T)`` of :func:`equation_of_state`.

    Returns
    -------
    float
        Number density in m^-3.

    Raises
    ------
    StateError
        For a real-gas density of a species without an equation of state (ozone), or at
        a state outside its phase domain: a condensed or solid state is refused, not
        returned.

    Examples
    --------
    >>> from optmat.conditions import ONE_ATM_0C
    >>> n_real = number_density("Xe", ONE_ATM_0C)
    >>> n_ideal = number_density("Xe", ONE_ATM_0C, ideal=True)
    >>> bool(n_real > n_ideal)  # xenon is attractive at 0 degrees C: Z < 1
    True
    """
    if ideal:
        return ideal_number_density(state)
    return equation_of_state(gas).number_density(state)


def pressure(
    gas: Gas | str,
    number_density: ArrayLike | float,
    temperature: ArrayLike | float,
    *,
    ideal: bool = False,
) -> Float[np.ndarray, "*shape"]:
    """Return the pressure of a pure gas at number densities and temperatures.

    Parameters
    ----------
    gas : Gas or str
        The species.
    number_density : array_like or float
        Number density, m^-3, > 0.
    temperature : array_like or float
        Temperature, K.
    ideal : bool, optional
        If True, ``N k_B T``; otherwise (the default) the inverse of the real-gas law,
        :meth:`EquationOfState.pressure`.

    Returns
    -------
    numpy.ndarray
        Pa.

    Examples
    --------
    >>> f"{float(pressure('Ar', 2.5e25, 294.0, ideal=True)):.2f}"
    '101477.70'
    """
    if ideal:
        return np.asarray(
            np.asarray(number_density, dtype=float) * K_B * np.asarray(temperature)
        )
    return equation_of_state(gas).pressure(number_density, temperature)


def compressibility(gas: Gas | str, state: State) -> float:
    """Return the compressibility factor ``Z = p / (N k_B T)`` of the real gas.

    Parameters
    ----------
    gas : Gas or str
        The species.
    state : State
        Pressure and temperature, in the phase domain.

    Returns
    -------
    float
        The dimensionless compressibility factor; 1 for an ideal gas.

    Examples
    --------
    >>> from optmat.conditions import ONE_ATM_0C
    >>> f"{compressibility('He', ONE_ATM_0C):.5f}"
    '1.00053'
    """
    return equation_of_state(gas).compressibility(state)

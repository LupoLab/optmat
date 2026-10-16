# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The dilute gases' transport properties: viscosity, heat capacity, molar mass.

OEP 1, sections 4.3 and 6.5. For each gas a :class:`TransportModel` gives

- the **dilute-gas shear viscosity** ``eta(T)``, Pa s: the zero-density limit, a
  function of temperature alone, as kinetic theory and Wilke's rule assume;
- the **ideal-gas molar isobaric heat capacity** ``c_p0(T)``, J/(mol K): the
  temperature derivative of the ideal-gas part of the gas's equation of state, not a
  value at a finite pressure; for the noble gases ``(5/2) R``;
- the **molar mass** ``M``, kg/mol, the one the equation of state uses.

The source of each gas's viscosity and heat capacity is curated in
``data/gas/transport.toml``, with its temperatures and its stated uncertainty:

- CoolProp's correlation, evaluated at 1e-6 mol/m^3, where it implements a curated
  paper (He, Ar, H2, N2, O2, CH4, CO2, SF6; each reproduces its paper's printed dilute
  form to below 1e-9; CH4's states no uncertainty, ``provisional``);
- Huber's dilute Chapman-Enskog form with his Lennard-Jones parameters and ``F_c`` (Ne,
  Kr, Xe, N2O; ``huber_2018``), the models of the NIST WebBook's tables, which serve as
  the check;
- Assael, Mixafendi and Wakeham's correlation for D2 (``assael_1987``): their
  hydrogen correlation (``assael_1986``, Eqs. 1 and 3) with the D2 mass, times the D2
  paper's quantum factor ``1 - exp(-f T*^(1/3))`` (its Eq. 7), which reproduces its
  Table 2, the check, within the table's rounding (TC-765, TC-839).

``c_p0`` keeps the CODATA convention of :mod:`optmat.gas.eos`: the equation's
dimensionless ``c_p0 / R`` times CODATA's ``R``, so the noble gases give ``(5/2) R``
exactly whatever gas constant their equation was published with.

:func:`mixture_viscosity` is Wilke's rule (``wilke_1950``), a kinetic-theory
approximation: 0.97 % on average for his binary mixtures, 1.9 % for multicomponent ones.
The ideal mixture's heat capacity is ``sum_i x_i c_p0,i``. What combines these for a
flow problem (a sound speed, a mean free path) stays in the caller.
"""

import functools
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike

from optmat._catalogue import records
from optmat.constants import K_B, N_A
from optmat.exceptions import DataError, OptmatError, StateError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.composition import composition
from optmat.gas.eos import (
    EquationOfState,
    _abstract_state,  # pyright: ignore[reportPrivateUsage] - one backend
    _coolprop,  # pyright: ignore[reportPrivateUsage]
    equation_of_state,
)
from optmat.gas.gases import Gas, as_gas
from optmat.identity import FrozenMapping, Identified, parse_status

__all__ = [
    "EVALUATION_VERSION",
    "QUANTITIES",
    "TRANSPORT_GASES",
    "TransportModel",
    "heat_capacity",
    "mixture_viscosity",
    "molar_mass",
    "transport_model",
    "viscosity",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in every transport model's identity;
raised whenever the code changes its output (OEP 1, section 4.8)."""

QUANTITIES: Final = ("viscosity", "heat_capacity")
"""The temperature-dependent quantities of a :class:`TransportModel`."""

_FORMS: Final = frozenset({"coolprop", "chapman_enskog", "assael"})
_DILUTE_DENSITY: Final = 1e-6
"""mol/m^3: where CoolProp's viscosity is taken as the dilute gas's. Its density terms
are below 1e-10 there (CoolProp reports helium's whole viscosity as one term)."""
_NEUFELD: Final = (0.3, 100.0)
"""The reduced temperatures of Neufeld et al.'s collision integral (Huber, Eq. 7)."""
_R_CODATA: Final = N_A * K_B


def _choices() -> Mapping[str, Mapping[str, Any]]:
    """Return ``data/gas/transport.toml``'s entries, by gas.

    Returns
    -------
    Mapping
        Each gas's choices.
    """
    with (GAS_DATA_DIRECTORY / "transport.toml").open("rb") as handle:
        return tomllib.load(handle)["gases"]


TRANSPORT_GASES: Final[tuple[str, ...]] = tuple(_choices())
"""The thirteen gases with a curated transport model."""


def _temperatures(temperature: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
    """Return temperatures as a float array.

    Parameters
    ----------
    temperature : array_like or float
        K.

    Returns
    -------
    numpy.ndarray
        K.
    """
    return np.asarray(temperature, dtype=float)


@dataclass(frozen=True)
class TransportModel(Identified):
    """A gas's dilute viscosity, ideal-gas heat capacity and molar mass.

    Attributes
    ----------
    material : str
        The gas.
    viscosity_form : str
        ``"coolprop"``, ``"chapman_enskog"`` or ``"assael"``.
    viscosity_source : str
        The source record of the viscosity, or ``""`` where none is in hand.
    viscosity_parameters : Mapping
        What the form evaluates: CoolProp's fluid and correlation key; Huber's
        constants; or Assael et al.'s constants, the molar mass (kg/mol) and the
        source of the collision integral (``collision_source``).
    viscosity_range : tuple of float
        The viscosity's temperatures, K.
    viscosity_uncertainty : str
        Its stated uncertainty, as the source states it.
    heat_capacity_source : str
        The source record of the ideal-gas part, or ``""`` where none is in hand.
    heat_capacity_range : tuple of float
        The heat capacity's temperatures, K: the equation of state's, or less where
        CoolProp's ideal part meets its source over less (O2).
    equation_of_state : EquationOfState
        The gas's equation of state, whose ideal part gives ``c_p0`` and whose
        temperatures bound it; its digest enters this model's.
    molar_mass : float
        kg/mol, the equation of state's.
    status : str
        ``validated``, or ``provisional`` where ``status_reason`` says why.
    """

    IDENTITY_KIND: ClassVar[str] = "transport"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.transport",)

    material: str
    viscosity_form: str
    viscosity_source: str
    viscosity_parameters: Mapping[str, Any]
    viscosity_range: tuple[float, float]
    viscosity_uncertainty: str
    heat_capacity_source: str
    heat_capacity_range: tuple[float, float]
    equation_of_state: EquationOfState
    molar_mass: float

    def __post_init__(self) -> None:
        """Check the form and the range.

        Raises
        ------
        DataError
            For an unknown form or an empty range.
        """
        super().__post_init__()
        if self.viscosity_form not in _FORMS:
            msg = f"{self.material}: unknown viscosity form {self.viscosity_form!r}"
            raise DataError(msg)
        low, high = self.viscosity_range
        if not 0 < low < high:
            msg = f"{self.material}: empty viscosity range {self.viscosity_range}"
            raise DataError(msg)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and source records.

        Returns
        -------
        tuple
            ``("transport/<gas>", gas, sources)``.
        """
        found = {
            self.viscosity_source,
            self.heat_capacity_source,
            str(self.viscosity_parameters.get("collision_source", "")),
        } - {""}
        return f"transport/{self.material}", self.material, tuple(sorted(found))

    def temperature_range(self, quantity: str) -> tuple[float, float]:
        """Return the temperatures a quantity is defined over.

        Parameters
        ----------
        quantity : str
            One of :data:`QUANTITIES`.

        Returns
        -------
        tuple of float
            K.

        Raises
        ------
        OptmatError
            For an unknown quantity.

        Examples
        --------
        >>> transport_model("CO2").temperature_range("viscosity")
        (200.0, 2000.0)
        """
        if quantity == "viscosity":
            return self.viscosity_range
        if quantity == "heat_capacity":
            return self.heat_capacity_range
        msg = f"unknown quantity {quantity!r}; one of {QUANTITIES}"
        raise OptmatError(msg)

    def _checked(self, quantity: str, temperature: ArrayLike | float) -> np.ndarray:
        """Return the temperatures, refusing any outside the quantity's range.

        Parameters
        ----------
        quantity : str
            One of :data:`QUANTITIES`.
        temperature : array_like or float
            K.

        Returns
        -------
        numpy.ndarray
            K.

        Raises
        ------
        StateError
            If a temperature is outside the range or not finite.
        """
        values = _temperatures(temperature)
        low, high = self.temperature_range(quantity)
        inside = np.isfinite(values) & (values >= low) & (values <= high)
        if not np.all(inside):
            bad = values[~inside].ravel()[0]
            name = quantity.replace("_", " ")
            msg = (
                f"the {name} of {self.material} is defined for {low:g}-{high:g} K, "
                f"not {bad:g} K"
            )
            raise StateError(msg)
        return values

    def viscosity(self, temperature: ArrayLike | float) -> Float[np.ndarray, "*shape"]:
        """Return the dilute-gas shear viscosity.

        Parameters
        ----------
        temperature : array_like or float
            K, within :meth:`temperature_range` ``("viscosity")``.

        Returns
        -------
        numpy.ndarray
            Pa s, the shape of ``temperature``.

        Raises
        ------
        StateError
            Outside the range.

        Examples
        --------
        >>> f"{float(transport_model('N2').viscosity(294.0)):.5e}"
        '1.75991e-05'
        """
        values = self._checked("viscosity", temperature)
        if self.viscosity_form == "chapman_enskog":
            return _chapman_enskog(self.viscosity_parameters, values)
        if self.viscosity_form == "assael":
            return _assael(self.viscosity_parameters, values)
        state = _abstract_state(self.viscosity_parameters["fluid"])
        out = np.empty(values.shape)
        for index, t in np.ndenumerate(values):
            state.update(_coolprop().DmolarT_INPUTS, _DILUTE_DENSITY, float(t))
            out[index] = state.viscosity()
        return out

    def heat_capacity(
        self, temperature: ArrayLike | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the ideal-gas molar isobaric heat capacity.

        Parameters
        ----------
        temperature : array_like or float
            K, within :meth:`temperature_range` ``("heat_capacity")``.

        Returns
        -------
        numpy.ndarray
            J/(mol K), the shape of ``temperature``: the equation's ``c_p0 / R``
            times CODATA's ``R``.

        Raises
        ------
        StateError
            Outside the range.

        Examples
        --------
        >>> f"{float(transport_model('Ar').heat_capacity(294.0)):.6f}"
        '20.786157'
        """
        values = self._checked("heat_capacity", temperature)
        equation = self.equation_of_state
        state = _abstract_state(equation.fluid)
        out = np.empty(values.shape)
        for index, t in np.ndenumerate(values):
            state.update(_coolprop().DmolarT_INPUTS, _DILUTE_DENSITY, float(t))
            out[index] = state.cp0molar() / equation.gas_constant * _R_CODATA
        return out


def _chapman_enskog(
    parameters: Mapping[str, Any], temperature: np.ndarray
) -> np.ndarray:
    """Return Huber's dilute-gas viscosity, Eqs. (6)-(7) times ``F_c``.

    Parameters
    ----------
    parameters : Mapping
        ``prefactor``, ``molar_mass`` (g/mol), ``epsilon_k`` (K), ``sigma_nm`` and
        ``f_c``.
    temperature : numpy.ndarray
        K.

    Returns
    -------
    numpy.ndarray
        Pa s.
    """
    reduced = temperature / parameters["epsilon_k"]
    omega = (
        1.16145 * reduced**-0.14874
        + 0.52487 * np.exp(-0.77320 * reduced)
        + 2.16178 * np.exp(-2.43787 * reduced)
    )
    micro = (
        parameters["prefactor"]
        * np.sqrt(parameters["molar_mass"] * temperature)
        / (parameters["sigma_nm"] ** 2 * omega)
        * parameters["f_c"]
    )
    return np.asarray(micro * 1e-6)


def _assael(parameters: Mapping[str, Any], temperature: np.ndarray) -> np.ndarray:
    """Return Assael, Mixafendi and Wakeham's dilute-gas viscosity of D2.

    Their Eq. (1), ``eta = (5/16) (m k T / pi)^(1/2) / (sigma^2 Omega*(T*))``, with the
    collision integral of their hydrogen correlation, ``ln Omega* = sum_i a_i (ln
    T*)^i`` (1986, Eq. 3 and Table 2), times the D2 paper's factor ``1 - exp(-f
    T*^(1/3))`` (1987, Eq. 7) at every temperature, as its Table 2 applies it
    (TC-765).

    Parameters
    ----------
    parameters : Mapping
        ``a0``-``a4``, ``sigma_nm``, ``epsilon_k`` (K), ``f`` and ``molar_mass``
        (kg/mol).
    temperature : numpy.ndarray
        K.

    Returns
    -------
    numpy.ndarray
        Pa s.
    """
    reduced = temperature / parameters["epsilon_k"]
    log_reduced = np.log(reduced)
    omega = np.exp(sum(parameters[f"a{i}"] * log_reduced**i for i in range(5)))
    mass = parameters["molar_mass"] / N_A
    eta = (
        5
        / 16
        * np.sqrt(mass * K_B * temperature / math.pi)
        / ((parameters["sigma_nm"] * 1e-9) ** 2 * omega)
    )
    return np.asarray(eta * (1 - np.exp(-parameters["f"] * reduced ** (1 / 3))))


def _formula(source: str, formula: str) -> dict[str, float]:
    """Return a curated formula's coefficients as numbers.

    Parameters
    ----------
    source : str
        The source record.
    formula : str
        The formula's id in it.

    Returns
    -------
    dict
        Coefficient name to value.

    Raises
    ------
    DataError
        If the record has no such formula.
    """
    for entry in records()[source].formulas:
        if entry["id"] == formula:
            coefficients = entry["coefficients"]
            assert isinstance(coefficients, Mapping)
            return {str(k): float(v) for k, v in coefficients.items()}
    msg = f"{source} has no formula {formula!r}"
    raise DataError(msg)


def _build(gas: str, entry: Mapping[str, Any]) -> TransportModel:
    """Return one gas's transport model from its entry in the choices file.

    Parameters
    ----------
    gas : str
        The gas.
    entry : Mapping
        Its ``[gases.<gas>]`` table.

    Returns
    -------
    TransportModel
        The model, with its status.
    """
    equation = equation_of_state(gas)
    fluid = equation.fluid
    form = str(entry["viscosity"])
    source = str(entry["viscosity_source"])
    parameters: dict[str, Any]
    if form == "chapman_enskog":
        parameters = _formula(source, str(entry["formula"]))
        low, high = equation.temperature_range
        epsilon = parameters["epsilon_k"]
        default = (max(low, _NEUFELD[0] * epsilon), min(high, _NEUFELD[1] * epsilon))
    elif form == "assael":
        collision = str(entry["collision_source"])
        parameters = {
            **_formula(collision, str(entry["collision_formula"])),
            "f": _formula(source, str(entry["formula"]))["f"],
            "molar_mass": float(_abstract_state(fluid).molar_mass()),
            "collision_source": collision,
        }
        default = equation.temperature_range
    else:
        key = _coolprop().get_fluid_param_string(fluid, "BibTeX-VISCOSITY")
        parameters = {"fluid": fluid, "reference": str(key)}
        default = equation.temperature_range
    stated = entry.get("viscosity_range_k")
    capacity = entry.get("heat_capacity_range_k", equation.temperature_range)
    status, reason = parse_status(entry, f"transport/{gas}", default="validated")
    return TransportModel(
        material=gas,
        viscosity_form=form,
        viscosity_source=source,
        viscosity_parameters=FrozenMapping(parameters),
        viscosity_range=default if stated is None else (stated[0], stated[1]),
        viscosity_uncertainty=str(entry["viscosity_uncertainty"]),
        heat_capacity_source=str(entry["heat_capacity_source"]),
        heat_capacity_range=(float(capacity[0]), float(capacity[1])),
        equation_of_state=equation,
        molar_mass=float(_abstract_state(fluid).molar_mass()),
        status=status,  # pyright: ignore[reportArgumentType] - parse_status checks it
        status_reason=reason,
    )


@functools.cache
def transport_model(gas: Gas | str) -> TransportModel:
    """Return a gas's transport model.

    Parameters
    ----------
    gas : Gas or str
        One of :data:`TRANSPORT_GASES`.

    Returns
    -------
    TransportModel
        Its viscosity, heat capacity and molar mass, with its identity.

    Raises
    ------
    OptmatError
        For a gas without one (O3, CO).

    Examples
    --------
    >>> transport_model("Xe").viscosity_source
    'huber_2018'
    """
    species = as_gas(gas).value
    choices = _choices()
    if species not in choices:
        msg = (
            f"no transport model for {species}; the gases with one are "
            f"{', '.join(TRANSPORT_GASES)}"
        )
        raise OptmatError(msg)
    return _build(species, choices[species])


def viscosity(
    gas: Gas | str, temperature: ArrayLike | float
) -> Float[np.ndarray, "*shape"]:
    """Return a gas's dilute-gas shear viscosity.

    Parameters
    ----------
    gas : Gas or str
        One of :data:`TRANSPORT_GASES`.
    temperature : array_like or float
        K.

    Returns
    -------
    numpy.ndarray
        Pa s; :meth:`TransportModel.viscosity`.

    Examples
    --------
    >>> f"{float(viscosity('He', 294.0)):.5e}"
    '1.96527e-05'
    """
    return transport_model(gas).viscosity(temperature)


def heat_capacity(
    gas: Gas | str, temperature: ArrayLike | float
) -> Float[np.ndarray, "*shape"]:
    """Return a gas's ideal-gas molar isobaric heat capacity.

    Parameters
    ----------
    gas : Gas or str
        One of :data:`TRANSPORT_GASES`.
    temperature : array_like or float
        K.

    Returns
    -------
    numpy.ndarray
        J/(mol K); :meth:`TransportModel.heat_capacity`.

    Examples
    --------
    >>> f"{float(heat_capacity('N2', 294.0)):.4f}"
    '29.1232'
    """
    return transport_model(gas).heat_capacity(temperature)


def molar_mass(gas: Gas | str) -> float:
    """Return a gas's molar mass, the one its equation of state uses.

    Parameters
    ----------
    gas : Gas or str
        One of :data:`TRANSPORT_GASES`.

    Returns
    -------
    float
        The molar mass, kg/mol.

    Examples
    --------
    >>> molar_mass("Ar")
    0.039948
    """
    return transport_model(gas).molar_mass


def _fractions(mixture: str | Mapping[Any, float]) -> dict[str, float]:
    """Return a mixture's gases and mole fractions.

    Parameters
    ----------
    mixture : str or Mapping
        A named composition (``"dry air"`` is ``"dry air (N2, O2, Ar)"``, as for the
        air model), or gases to mole fractions summing to one.

    Returns
    -------
    dict
        Gas to its positive mole fraction.

    Raises
    ------
    OptmatError
        If a fraction is negative or not finite, or the fractions do not sum to one.
    """
    if isinstance(mixture, str):
        name = "dry air (N2, O2, Ar)" if mixture == "dry air" else mixture
        return {c.material: c.mole_fraction for c in composition(name)}
    merged: dict[str, float] = {}
    for gas, fraction in mixture.items():
        value = float(fraction)
        if not (math.isfinite(value) and value >= 0):
            msg = f"mole fraction of {gas} must be finite and not negative, not {value}"
            raise OptmatError(msg)
        species = as_gas(gas).value
        merged[species] = merged.get(species, 0.0) + value
    total = math.fsum(merged.values())
    # A distribution normalised in floating point sums to one within a few ulps.
    if abs(total - 1) > 1e-12:
        msg = f"mole fractions must sum to one, not {total:.17g}"
        raise OptmatError(msg)
    return {gas: x for gas, x in merged.items() if x > 0}


def mixture_viscosity(
    mixture: str | Mapping[Any, float], temperature: ArrayLike | float
) -> Float[np.ndarray, "*shape"]:
    """Return a gas mixture's dilute viscosity by Wilke's rule.

    ``eta = sum_i x_i eta_i / sum_j x_j Phi_ij`` with
    ``Phi_ij = [1 + (eta_i/eta_j)^(1/2) (M_j/M_i)^(1/4)]^2 / [8 (1 + M_i/M_j)]^(1/2)``
    (Wilke 1950, Eqs. 13-14; average error 0.97 % for his binaries, 1.9 % for
    multicomponent mixtures).

    Parameters
    ----------
    mixture : str or Mapping
        A named composition, or gases to mole fractions summing to one.
    temperature : array_like or float
        K, within every constituent's viscosity range.

    Returns
    -------
    numpy.ndarray
        Pa s, the shape of ``temperature``.

    Examples
    --------
    >>> f"{float(mixture_viscosity({'He': 0.5, 'Ar': 0.5}, 294.0)):.4e}"
    '2.3335e-05'
    """
    fractions = _fractions(mixture)
    gases = list(fractions)
    x = [fractions[g] for g in gases]
    etas = [viscosity(g, temperature) for g in gases]
    masses = [molar_mass(g) for g in gases]
    total = np.zeros(np.shape(etas[0]))
    for i in range(len(gases)):
        denominator = sum(
            x[j]
            * (1 + np.sqrt(etas[i] / etas[j]) * (masses[j] / masses[i]) ** 0.25) ** 2
            / math.sqrt(8 * (1 + masses[i] / masses[j]))
            for j in range(len(gases))
        )
        total = total + x[i] * etas[i] / denominator
    return np.asarray(total)

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Published formulas of solids as models.

A curated source's ``[[formulas]]`` entry for a solid becomes a model when it has a
``form`` this module knows; its coefficients are the printed strings, corrected by the
entry's ``corrections`` with the evidence for each (as for the gases):

- ``n2_sellmeier``: n^2 - 1 = A + sum_i B_i lambda^2 / (lambda^2 - C_i), lambda the
  vacuum wavelength in um, C_i in um^2 (terms ``A``, optional, ``B``, ``C``);
- ``n2_sellmeier_wavelength``: the same with the poles printed as wavelengths,
  C_i = L_i^2 (terms ``B``, ``L``, in um);
- ``n2_sellmeier_temperature``: n^2 - 1 = sum_i S_i(T) lambda^2 / (lambda^2 -
  lambda_i(T)^2), with S_i(T) = sum_j S_ij T^j and lambda_i(T) = sum_j lambda_ij T^j
  (terms ``S1``, ``S2``, ... and ``L1``, ``L2``, ..., lambda in um, T in K), as
  Leviton & Frey fit fused silica; the model is evaluated at the entry's
  ``reference_temperature_k``, and :meth:`FormulaModel.at_temperature` gives it at
  another temperature within ``temperature_range_k``;
- ``udm``: Franta et al.'s universal dispersion model (:mod:`optmat.solid.udm`), with
  the values of its newAD2 parameter file as the terms: ``attributes`` (the file's
  first line, ``he:ex=10:ut:ph=52``) and each parameter in newAD2's names and units.

A formula whose variable is the wavelength in standard air says so with
``wavelength_medium = "standard_air"``; the model converts. A formula whose index is
relative to the room air of the measurement, as minimum-deviation refractometry in air
gives it, says so with ``index_medium = "air"`` and the air's ``air_temperature_k`` and
``air_pressure_pa`` (101 325 Pa when the source states none); the model multiplies by
the index of that air (Edlen's dry air, :func:`optmat.conditions.air_refractivity`) to
give the absolute index (TC-398). A formula whose index is relative to nitrogen, as
Daimon & Masumura measure CaF2, says so with ``index_medium = "nitrogen"`` and
``nitrogen_temperature_k`` and ``nitrogen_pressure_pa``; the model multiplies by the
index of that nitrogen from optmat's recommended N2 model at the ideal-gas density
(TC-545). A Sellmeier formula is transparent (real epsilon) and
meaningless near its poles; its range is the one its source states.
"""

import dataclasses
import functools
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, ClassVar

import numpy as np
from jaxtyping import Complex, Float

from optmat._catalogue import formula_sources, source
from optmat.conditions import State, air_refractivity, standard_air_refractivity
from optmat.constants import E_CHARGE, HBAR, STD_ATMOSPHERE, C
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.materials import materials
from optmat.solid.models._defaults import curated
from optmat.solid.models.base import Model
from optmat.solid.udm import UniversalDispersionModel, parse_attributes
from optmat.solid.udm import from_parameters as udm_from_parameters

__all__ = ["FORMS", "FormulaModel", "UdmModel", "published_model", "published_models"]

FORMS = frozenset(
    {"n2_sellmeier", "n2_sellmeier_wavelength", "n2_sellmeier_temperature", "udm"}
)
"""The formula forms solids use."""

_UM = 1e-6

_AIR_STEP = 1e-4
"""Relative step of the central difference of the reference air's refractivity.

For Edlen's air its truncation error is below 8e-9 of the derivative of the
refractivity from 250 nm to 3.5 um, and its rounding error some 1e-10: both below 1e-11
of a solid's dn/domega, whose air term is below 1e-3 of it (TC-838).
"""


@dataclass(frozen=True)
class FormulaModel(Model):
    """A published Sellmeier formula of a solid.

    Attributes
    ----------
    form : str
        ``n2_sellmeier``, ``n2_sellmeier_wavelength`` or ``n2_sellmeier_temperature``.
    terms : dict of str to tuple of float
        The coefficients, after corrections.
    wavelength_medium : str
        ``vacuum``, or ``standard_air`` for a formula whose variable is the wavelength
        in standard air (converted with Edlen's formula).
    air : State or None
        The state of the gas the formula's index is relative to (room air, or
        nitrogen), or None for an absolute index.
    reference_gas : str
        ``air`` or ``nitrogen``: the gas :attr:`air` is the state of.
    temperature : float or None
        The sample temperature in K of an ``n2_sellmeier_temperature`` formula.
    temperature_range : tuple of float or None
        The temperatures in K its source supports.
    """

    form: str = ""
    terms: Mapping[str, tuple[float, ...]] = field(default_factory=dict, repr=False)
    wavelength_medium: str = "vacuum"
    air: State | None = None
    reference_gas: str = "air"
    temperature: float | None = None
    temperature_range: tuple[float, float] | None = None

    EVALUATION_VERSION: ClassVar[int] = 2
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.conditions",)
    INPUTS: ClassVar[tuple[str, ...]] = ("intervals", "terms")

    @classmethod
    def sellmeier(
        cls,
        material: str,
        coefficients: Mapping[str, Sequence[float]],
        wavelength_range: tuple[float, float],
        *,
        reference_medium: str | State = "vacuum",
        citation: str = "",
        name: str = "",
    ) -> "FormulaModel":
        """Return a Sellmeier formula from its coefficients.

        ``n^2 = 1 + A + sum_i B_i lambda^2 / (lambda^2 - C_i)``, lambda the vacuum
        wavelength in um and C_i in um^2.

        Parameters
        ----------
        material : str
            The material, the identity's (a curated id or a caller's name).
        coefficients : Mapping
            ``B`` and ``C`` (equal lengths), and optionally ``A``.
        wavelength_range : tuple of float
            The vacuum wavelengths the formula holds over, m.
        reference_medium : str or State, optional
            ``"vacuum"`` for an absolute index, or the state of the room air the
            index is relative to.
        citation : str, optional
            Where the coefficients come from.
        name : str, optional
            The model's name; ``user/<material>`` by default.

        Returns
        -------
        FormulaModel
            The model, status ``user``.

        Raises
        ------
        ValueError
            If ``B`` and ``C`` differ in length or a term is missing.

        Examples
        --------
        >>> silica = FormulaModel.sellmeier(
        ...     "SiO2", {"B": (0.6961663, 0.4079426, 0.8974794),
        ...              "C": (0.0684043**2, 0.1162414**2, 9.896161**2)},
        ...     (210e-9, 3.71e-6), citation="Malitson (1965)")
        >>> f"{float(silica.refractive_index(1e-6).real):.5f}"
        '1.45042'
        """
        terms = {
            k: tuple(float(v) for v in values) for k, values in coefficients.items()
        }
        if set(terms) - {"A", "B", "C"} or len(terms.get("B", ())) != len(
            terms.get("C", (0.0,))
        ):
            msg = "a Sellmeier formula takes B and C of equal lengths, and A"
            raise ValueError(msg)
        return cls(
            material=material,
            name=name or f"user/{material}",
            citation=citation,
            wavelength_range=wavelength_range,
            kind="user",
            form="n2_sellmeier",
            terms=terms,
            air=None if reference_medium == "vacuum" else reference_medium,  # pyright: ignore[reportArgumentType] - a State
        )

    @functools.cached_property
    def reference_model(self) -> object:
        """The gas model of a formula relative to nitrogen, bound once.

        A solid formula relative to nitrogen needs the gas subpackage only when it is
        evaluated, or its identity asked for.

        Returns
        -------
        optmat.gas.models.base.Model or None
            The recommended model of N2, or None for a formula relative to air or
            absolute.
        """
        if self.air is None or self.reference_gas != "nitrogen":
            return None
        from optmat.gas.models.registry import default_model  # noqa: PLC0415

        return default_model("N2")

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return the gas model a formula relative to nitrogen reads.

        Returns
        -------
        Mapping
            ``{"reference_model": ...}``, None where there is none.
        """
        return {"reference_model": self.reference_model}

    def at_temperature(self, temperature: float) -> Model:
        """Return the model at another sample temperature.

        Parameters
        ----------
        temperature : float
            The sample temperature, K, within the source's ``temperature_range``.

        Returns
        -------
        FormulaModel
            The same formula evaluated at ``temperature`` (typed as its base, for
            the runtime type checker).

        Raises
        ------
        DataError
            If the formula has no temperature dependence.

        Warns
        -----
        ModelRangeWarning
            If ``temperature`` is outside the range its source supports.
        """
        if self.form != "n2_sellmeier_temperature":
            msg = f"{self.name} has no temperature dependence"
            raise DataError(msg)
        low, high = self.temperature_range or (temperature, temperature)
        if not low <= temperature <= high:
            msg = f"{self.name} evaluated at {temperature} K, outside {low}-{high} K"
            warnings.warn(msg, ModelRangeWarning, stacklevel=2)
        return dataclasses.replace(self, temperature=float(temperature))

    def _air_refractivity(self, omega: np.ndarray) -> np.ndarray:
        """Return the refractivity n - 1 of the reference air of a relative index.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The refractivity n_air - 1 at the vacuum wavelength 2 pi c / omega.
        """
        assert self.air is not None
        wavenumber = omega / (2 * math.pi * C)
        if self.reference_gas == "nitrogen":
            # The ideal-gas density: N2's compressibility at 1 atm and room
            # temperature, Z - 1 = -2e-4, changes n by below 1e-7 (TC-545).
            reference = self.reference_model
            assert reference is not None
            return np.asarray(
                reference.refractivity(  # pyright: ignore[reportAttributeAccessIssue]
                    1 / wavenumber, self.air, ideal=True, warn=False
                )
            )
        return np.asarray(air_refractivity(wavenumber, self.air))

    def _air_index_squared(self, omega: np.ndarray) -> np.ndarray:
        """Return the squared index of the reference air, or 1 for an absolute index.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The squared index n_air^2 at the vacuum wavelength 2 pi c / omega.
        """
        if self.air is None:
            return np.ones_like(omega)
        return np.asarray((1 + self._air_refractivity(omega)) ** 2)

    def _lambda2(self, omega: np.ndarray) -> np.ndarray:
        """Return the formula's variable, lambda^2 in um^2.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The squared wavelength in the formula's medium, um^2.
        """
        vacuum = 2 * math.pi * C / omega
        if self.wavelength_medium == "standard_air":
            vacuum = vacuum / (1 + standard_air_refractivity(1 / vacuum))
        return np.asarray((vacuum / _UM) ** 2)

    def _poles(self) -> tuple[np.ndarray, np.ndarray, float]:
        """Return the strengths B, the poles C (um^2) and the constant A.

        Returns
        -------
        tuple
            B, C and A.
        """
        if self.form == "n2_sellmeier_temperature":
            assert self.temperature is not None
            count = sum(1 for key in self.terms if key.startswith("S"))
            powers = [self.temperature**j for j in range(len(self.terms["S1"]))]
            strengths = np.array(
                [np.dot(self.terms[f"S{i}"], powers) for i in range(1, count + 1)]
            )
            poles = np.array(
                [np.dot(self.terms[f"L{i}"], powers) ** 2 for i in range(1, count + 1)]
            )
            return strengths, poles, 0.0
        strengths = np.asarray(self.terms["B"])
        if self.form == "n2_sellmeier_wavelength":
            poles = np.asarray(self.terms["L"]) ** 2
        else:
            poles = np.asarray(self.terms["C"])
        constant = float(sum(self.terms.get("A", (0.0,))))
        return strengths, poles, constant

    def _permittivity(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return epsilon = n^2 from the formula (real).

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The permittivity.
        """
        lam2 = self._lambda2(omega)
        b, c, a = self._poles()
        eps = 1 + a + sum(bi * lam2 / (lam2 - ci) for bi, ci in zip(b, c, strict=True))
        return np.asarray(eps * self._air_index_squared(omega), dtype=complex)

    def _dpermittivity_domega(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d epsilon / d omega from the formula.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The derivative, s/rad.
        """
        lam2 = self._lambda2(omega)
        b, c, _ = self._poles()
        # d/d(lambda^2) of B lambda^2/(lambda^2 - C) is -B C/(lambda^2 - C)^2.
        slope = sum(-bi * ci / (lam2 - ci) ** 2 for bi, ci in zip(b, c, strict=True))
        # The air's dispersion enters through its refractivity, differenced on its
        # own: r ~ 3e-4 keeps its digits, where 1 + r or (1 + r)^2, rounded at 1e-16
        # of 1, carried 1e-8 of their change over the step into the derivative, a
        # jitter between neighbouring frequencies (TC-838).
        h = _AIR_STEP * omega
        if self.wavelength_medium == "standard_air":
            # lambda_std^2 = lambda_vac^2 / (1 + r_s)^2.
            standard = standard_air_refractivity(omega / (2 * math.pi * C))
            dstandard = (
                standard_air_refractivity((omega + h) / (2 * math.pi * C))
                - standard_air_refractivity((omega - h) / (2 * math.pi * C))
            ) / (2 * h)
            dlam2 = -2 * lam2 * (1 / omega + dstandard / (1 + standard))
        else:
            dlam2 = -2 * lam2 / omega
        derivative = slope * dlam2 * self._air_index_squared(omega)
        if self.air is not None:
            # The relative epsilon times d(n_air^2)/d omega = 2 n_air dr/d omega.
            b, c, a = self._poles()
            relative = (
                1
                + a
                + sum(bi * lam2 / (lam2 - ci) for bi, ci in zip(b, c, strict=True))
            )
            refractivity = self._air_refractivity(omega)
            dr = (
                self._air_refractivity(np.asarray(omega + h))
                - self._air_refractivity(np.asarray(omega - h))
            ) / (2 * h)
            derivative = derivative + relative * 2 * (1 + refractivity) * dr
        return np.asarray(derivative, dtype=complex)


@dataclass(frozen=True)
class UdmModel(Model):
    """Franta et al.'s universal dispersion model of a solid, as a model.

    Attributes
    ----------
    udm : UniversalDispersionModel
        The dielectric function in photon energy.
    """

    udm: UniversalDispersionModel | None = field(default=None, repr=False)

    EVALUATION_VERSION: ClassVar[int] = 1

    def _energy(self, omega: np.ndarray) -> np.ndarray:
        """Return the photon energy in eV.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The photon energy hbar omega / e.
        """
        return np.asarray(HBAR * omega / E_CHARGE)

    def _permittivity(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the UDM's epsilon.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The complex permittivity.
        """
        assert self.udm is not None
        energy = self._energy(np.asarray(omega, dtype=float))
        return np.asarray(self.udm.permittivity(energy.ravel()).reshape(energy.shape))

    def _dpermittivity_domega(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the UDM's d epsilon / d omega (exact).

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The derivative, s/rad.
        """
        assert self.udm is not None
        energy = self._energy(np.asarray(omega, dtype=float))
        slope = self.udm.dpermittivity_denergy(energy.ravel()).reshape(energy.shape)
        return np.asarray(slope * HBAR / E_CHARGE)


def _corrected(
    terms: dict[str, list[str]], corrections: Sequence[Mapping[str, object]], where: str
) -> None:
    """Apply a formula's corrections in place, checking each printed coefficient.

    Parameters
    ----------
    terms : dict of str to list of str
        The printed coefficients.
    corrections : sequence of mapping
        Entries with ``term``, ``index``, ``printed``, ``corrected``, ``reason`` and
        ``evidence``.
    where : str
        The formula, for messages.

    Raises
    ------
    DataError
        If a correction is incomplete or its printed value is not the recorded one.
    """
    for correction in corrections:
        missing = {"term", "index", "printed", "corrected", "reason", "evidence"}
        missing -= set(correction)
        if missing:
            msg = f"{where}: correction lacks {sorted(missing)}"
            raise DataError(msg)
        term, index = str(correction["term"]), int(correction["index"])  # pyright: ignore[reportArgumentType]
        if term not in terms or not 0 <= index < len(terms[term]):
            msg = f"{where}: correction of missing coefficient {term}[{index}]"
            raise DataError(msg)
        if terms[term][index] != correction["printed"]:
            msg = (
                f"{where}: correction's printed {term}[{index}] is not the recorded one"
            )
            raise DataError(msg)
        terms[term][index] = str(correction["corrected"])


def _range(entry: Mapping[str, object], where: str) -> tuple[float, float]:
    """Return a formula's stated range in m.

    Parameters
    ----------
    entry : mapping
        The formula entry, with ``range_nm = [low, high]``.
    where : str
        The formula, for messages.

    Returns
    -------
    tuple of float
        The range, m.

    Raises
    ------
    DataError
        If the range is missing or malformed.
    """
    value = entry.get("range_nm")
    if not isinstance(value, Sequence) or len(value) != 2:
        msg = f"{where}: a formula with a form needs range_nm = [low, high]"
        raise DataError(msg)
    low, high = (float(v) * 1e-9 for v in value)  # pyright: ignore[reportArgumentType]
    return low, high


def _model(source_id: str, citation: str, entry: Mapping[str, object]) -> Model:
    """Build one solid formula model from its record entry.

    Parameters
    ----------
    source_id : str
        The source.
    citation : str
        The source's citation.
    entry : mapping
        The ``[[formulas]]`` entry.

    Returns
    -------
    Model
        The model.

    Raises
    ------
    DataError
        If the entry is malformed.
    """
    where = f"{source_id}/{entry.get('id')}"
    form = str(entry["form"])
    parameters = entry.get("parameters", {})
    axis = str(parameters.get("axis", "")) if isinstance(parameters, Mapping) else ""
    common = {
        "material": str(entry["material"]),
        "name": where,
        "citation": citation,
        "wavelength_range": _range(entry, where),
        "kind": str(entry.get("kind", "published")),
        "axis": axis,
    }
    raw = entry.get("terms")
    if not isinstance(raw, Mapping):
        msg = f"{where}: a formula with a form needs a terms table"
        raise DataError(msg)
    if form == "udm":
        values = {str(k): float(v) for k, v in raw.items() if k != "attributes"}
        attributes = parse_attributes(str(raw.get("attributes", "")), where)
        return UdmModel(**common, udm=udm_from_parameters(attributes, values))
    terms = {
        str(k): [str(x) for x in (v if isinstance(v, list) else [v])]
        for k, v in raw.items()
    }
    _corrected(terms, entry.get("corrections", []), where)  # pyright: ignore[reportArgumentType]
    numbers = {k: tuple(float(x) for x in v) for k, v in terms.items()}
    temperature, temperature_range = (
        _temperature(entry, numbers, where)
        if form == "n2_sellmeier_temperature"
        else (None, None)
    )
    if temperature is None:
        needed = {"B", "L"} if form == "n2_sellmeier_wavelength" else {"B", "C"}
        if not needed <= set(terms) or len(terms["B"]) != len(
            terms[(needed - {"B"}).pop()]
        ):
            msg = f"{where}: a {form} needs the terms {sorted(needed)} of equal length"
            raise DataError(msg)
    medium = str(entry.get("wavelength_medium", "vacuum"))
    if medium not in {"vacuum", "standard_air"}:
        msg = f"{where}: unknown wavelength_medium {medium!r}"
        raise DataError(msg)
    return FormulaModel(
        **common,
        form=form,
        terms=numbers,
        wavelength_medium=medium,
        air=_reference_air(entry, where),
        reference_gas=str(entry.get("index_medium", "vacuum")),
        temperature=temperature,
        temperature_range=temperature_range,
    )


def _reference_air(entry: Mapping[str, object], where: str) -> State | None:
    """Return the state of the gas a formula's index is relative to, or None.

    Parameters
    ----------
    entry : mapping
        The formula entry: ``index_medium`` (``vacuum``, the default, ``air`` or
        ``nitrogen``) and, for a gas, ``<medium>_temperature_k`` and
        ``<medium>_pressure_pa`` (101 325 Pa if absent).
    where : str
        The formula, for messages.

    Returns
    -------
    State or None
        The gas's state; None for an absolute index.

    Raises
    ------
    DataError
        If the medium is unknown or an air medium has no temperature.
    """
    medium = str(entry.get("index_medium", "vacuum"))
    if medium == "vacuum":
        return None
    if medium not in {"air", "nitrogen"} or f"{medium}_temperature_k" not in entry:
        msg = f"{where}: index_medium must be vacuum, or air or nitrogen with "
        msg += "<medium>_temperature_k"
        raise DataError(msg)
    pressure = float(entry.get(f"{medium}_pressure_pa", STD_ATMOSPHERE))  # pyright: ignore[reportArgumentType]
    return State(pressure, float(entry[f"{medium}_temperature_k"]))  # pyright: ignore[reportArgumentType]


def _temperature(
    entry: Mapping[str, object], terms: Mapping[str, tuple[float, ...]], where: str
) -> tuple[float, tuple[float, float]]:
    """Check an ``n2_sellmeier_temperature`` entry and return its temperatures.

    Parameters
    ----------
    entry : mapping
        The formula entry, with ``reference_temperature_k`` and
        ``temperature_range_k = [low, high]``.
    terms : mapping
        Its coefficients: ``S1``, ``S2``, ... and ``L1``, ``L2``, ..., each a
        polynomial in T (lowest power first), all of one length.
    where : str
        The formula, for messages.

    Returns
    -------
    tuple
        The reference temperature and the supported range, K.

    Raises
    ------
    DataError
        If the terms or the temperatures are malformed.
    """
    count = sum(1 for key in terms if key.startswith("S"))
    names = {f"S{i}" for i in range(1, count + 1)} | {
        f"L{i}" for i in range(1, count + 1)
    }
    lengths = {len(terms[name]) for name in names if name in terms}
    if count == 0 or set(terms) != names or len(lengths) != 1:
        msg = f"{where}: n2_sellmeier_temperature needs S1..Sn, L1..Ln of one length"
        raise DataError(msg)
    span = entry.get("temperature_range_k")
    if not isinstance(span, Sequence) or len(span) != 2:
        msg = f"{where}: an n2_sellmeier_temperature needs temperature_range_k"
        raise DataError(msg)
    low, high = (float(v) for v in span)  # pyright: ignore[reportArgumentType]
    reference = float(entry.get("reference_temperature_k", 0.0))  # pyright: ignore[reportArgumentType]
    if not low <= reference <= high:
        msg = f"{where}: reference_temperature_k must lie in temperature_range_k"
        raise DataError(msg)
    return reference, (low, high)


@functools.cache
def published_model(name: str) -> Model:
    """Return one solid formula model, reading only its source's record.

    :func:`published_models` reads every record with formulas, which takes seconds; a
    lookup by name needs one.

    Parameters
    ----------
    name : str
        ``<source>/<formula>``.

    Returns
    -------
    Model
        The model, the same as :func:`published_models` holds under that name.

    Raises
    ------
    KeyError
        If no solid formula has that name.

    Examples
    --------
    >>> published_model("malitson_1965/sellmeier").material
    'SiO2'
    """
    source_id, _, formula_id = name.partition("/")
    try:
        record = source(source_id)
    except KeyError:
        raise KeyError(name) from None
    solids = {key for key, value in materials().items() if value.is_solid}
    for entry in record.formulas:
        if (
            entry.get("id") == formula_id
            and entry.get("form") in FORMS
            and str(entry.get("material")) in solids
        ):
            return curated({name: _model(record.id, record.citation, entry)})[name]
    raise KeyError(name)


@functools.cache
def published_models() -> dict[str, Model]:
    """Return every solid formula model defined in the source records, by name.

    Returns
    -------
    dict of str to Model
        The models, named ``<source>/<formula>``.
    """
    solids = {key for key, value in materials().items() if value.is_solid}
    return curated(
        {
            f"{source.id}/{entry['id']}": _model(source.id, source.citation, entry)
            for source in formula_sources()
            for entry in source.formulas
            if entry.get("form") in FORMS and str(entry.get("material")) in solids
        }
    )

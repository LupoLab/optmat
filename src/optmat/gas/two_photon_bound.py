# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The bound two-photon levels of the noble gases (OEP 6, Phase 2a: curation only).

Below half the ionisation energy, two photons of a field can excite a noble gas's
ground state to an even-parity level. This module states what is curated about those
levels, per gas (:func:`two_photon_bound_levels`):

- the even levels with :math:`J = 0, 1, 2` below the ionisation energy, from the NIST
  level lists held as returned (``nist_asd_2026``; helium's singlets only). For one
  linear polarisation and equal photons, :math:`J = 0` and 2 are reached; :math:`J = 1`
  only with unequal photons or another polarisation;
- each level's curated two-photon strengths, calculated or measured, in one
  convention: :math:`\delta^{TP}` per level (summed over :math:`M_J`), parallel linear
  polarisation, equal photons, in atomic units, and the frequency-integrated cross
  section

  .. math:: \sigma_0^{(2)} = 2\pi^3 \alpha^2 a_0^4\, \omega_\mathrm{au}^2\,
     \delta^{TP},

  with the line profile normalised in the two-photon angular frequency, for coherent
  light (TC-952). A measured :math:`\sigma_0^{(2)}` is corrected for the photon
  statistics its source used;
- each level's curated photoionisation cross sections.

It is a record of facts with no evaluator. A line profile, the excited population and
the (2+1) ionisation that follows need the theory OEP 6 asks for in Q9, and are deferred
to its Phase 2b (TC-945). Helium's and xenon's strengths are checked against independent
values (``assessment = "checked"``). Neon's, argon's and krypton's are ``"unresolved"``,
each with its reason.

Examples
--------
>>> xenon = two_photon_bound_levels("Xe")
>>> level = xenon.level("5p5.(2P*<3/2>).6p", "2[1/2]", 0)
>>> round(level.wavelength * 1e9, 2)  # the exciting photon's vacuum wavelength, nm
249.63
>>> [(s.source, s.kind) for s in xenon.strengths_of(level)]
[('henriksson_2005', 'calculated'), ('kroll_1990', 'measured')]
"""

import functools
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from typing import Any, ClassVar, Final

from optmat.constants import AU_ANGULAR_FREQUENCY, BOHR_RADIUS, E_CHARGE, EPS0, HBAR, C
from optmat.exceptions import DataError, UnavailableModelError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.identity import Identified, source_ids

__all__ = [
    "BOUND_DERIVATION_VERSION",
    "SIGMA0_PER_DELTA",
    "BoundLevel",
    "BoundPhotoionisation",
    "BoundStrength",
    "TwoPhotonBoundLevels",
    "sigma0_from_delta",
    "two_photon_bound_levels",
]

BOUND_DERIVATION_VERSION: Final = 1
"""The version of the stored records' derivation."""
_ALPHA: Final = E_CHARGE**2 / (4 * math.pi * EPS0 * HBAR * C)
"""The fine-structure constant, from the package's CODATA constants."""
SIGMA0_PER_DELTA: Final = 2 * math.pi**3 * _ALPHA**2 * BOHR_RADIUS**4
r""":math:`\sigma_0^{(2)} / (\omega_\mathrm{au}^2 \delta^{TP})`, m^4 (TC-952)."""

_SELECTION: Final = GAS_DATA_DIRECTORY / "two_photon_bound.toml"
_DIRECTORY: Final = GAS_DATA_DIRECTORY / "models" / "two_photon_bound"
_STANDING: Final = frozenset({"status", "status_reason", "modified", "base"})
_GASES: Final = ("He", "Ne", "Ar", "Kr", "Xe")


def sigma0_from_delta(delta: float, omega: float) -> float:
    r"""Return the frequency-integrated cross section of a two-photon strength.

    Parameters
    ----------
    delta : float
        :math:`\delta^{TP}` per level, parallel linear polarisation, equal photons,
        atomic units.
    omega : float
        The photon's angular frequency, rad/s: half the level's over hbar.

    Returns
    -------
    float
        :math:`\sigma_0^{(2)} = 2\pi^3\alpha^2 a_0^4 \omega_\mathrm{au}^2 \delta`, m^4,
        the profile normalised in the two-photon angular frequency, coherent light
        (TC-952).

    Examples
    --------
    >>> f"{sigma0_from_delta(1.0, AU_ANGULAR_FREQUENCY):.5e}"
    '2.58949e-44'
    """
    return SIGMA0_PER_DELTA * (omega / AU_ANGULAR_FREQUENCY) ** 2 * delta


@dataclass(frozen=True)
class BoundLevel:
    """An even level below the ionisation energy, as NIST lists it.

    Attributes
    ----------
    configuration : str
        NIST's configuration, e.g. ``"5p5.(2P*<3/2>).6p"``.
    term : str
        NIST's term, e.g. ``"2[1/2]"``.
    j : int
        The total angular momentum.
    energy : float
        Above the ground level, J.
    printed : str
        NIST's level, cm^-1, as printed (with any brackets).
    reach : str
        ``"linear"`` (J = 0, 2: reached by two equal, linearly polarised photons) or
        ``"other"`` (J = 1: only with unequal photons or another polarisation).
    """

    configuration: str
    term: str
    j: int
    energy: float
    printed: str
    reach: str

    @property
    def wavelength(self) -> float:
        """The vacuum wavelength of each of two equal photons that reach it, m.

        Returns
        -------
        float
            :math:`2 hc / E`.
        """
        return 2 * (2 * math.pi * HBAR * C) / self.energy

    @property
    def key(self) -> str:
        """The level as ``"configuration|term|J"``.

        Returns
        -------
        str
            The key the selection file uses.
        """
        return f"{self.configuration}|{self.term}|{self.j}"


@dataclass(frozen=True)
class BoundStrength:
    r"""A curated two-photon strength of a level, in the canonical convention.

    Attributes
    ----------
    level : str
        The level's key (``BoundLevel.key``), or empty where unassigned.
    unassigned : str
        Why it is not assigned to a level, or empty.
    dataset : str
        The curated dataset.
    row : str
        Its row.
    source : str
        The source record.
    kind : str
        ``"calculated"`` or ``"measured"``.
    delta : float or None
        :math:`\delta^{TP}` per level, a.u.: as printed, or derived from a
        cross section; None for an unassigned cross section.
    sigma0 : float or None
        :math:`\sigma_0^{(2)}`, m^4, at the level's NIST photon frequency; None where
        unassigned.
    convention : str
        How it was obtained: ``"per_level"``, ``"one_component"`` (times 4/3) or
        ``"sigma0"``.
    derived : bool
        True unless the stored quantity is as printed.
    reason : str
        Why this level.
    """

    level: str
    unassigned: str
    dataset: str
    row: str
    source: str
    kind: str
    delta: float | None
    sigma0: float | None
    convention: str
    derived: bool
    reason: str


@dataclass(frozen=True)
class BoundPhotoionisation:
    """A curated photoionisation cross section of an excited level.

    Attributes
    ----------
    level : str
        The level's key, or empty where unassigned.
    unassigned : str
        Why it is not assigned, or empty.
    dataset : str
        The curated dataset.
    source : str
        The source record.
    omega : tuple of float or None
        The ionising photon's angular frequency, rad/s; None where the source states
        none (a threshold value at an unstated energy).
    sigma : tuple of float
        The cross section, m^2.
    use : str
        The dataset's use.
    """

    level: str
    unassigned: str
    dataset: str
    source: str
    omega: tuple[float | None, ...]
    sigma: tuple[float, ...]
    use: str


@dataclass(frozen=True, eq=False)
class TwoPhotonBoundLevels(Identified):
    """The bound two-photon levels of a noble gas and what is curated about them.

    Attributes
    ----------
    gas : Gas
        The species.
    derivation_version : int
        :data:`BOUND_DERIVATION_VERSION` when it was derived.
    ionisation : float
        The ionisation energy, J.
    levels : tuple of BoundLevel
        The even J = 0, 1, 2 levels below it, by energy.
    strengths : tuple of BoundStrength
        The curated two-photon strengths.
    photoionisation : tuple of BoundPhotoionisation
        The curated photoionisation cross sections of the levels.
    assessment : str
        ``"checked"`` or ``"unresolved"``.
    assessment_reason : str
        Why.
    input_digests : Mapping of str to str
        SHA-256 digests of what it was derived from.
    """

    gas: Gas
    derivation_version: int
    ionisation: float
    levels: tuple[BoundLevel, ...]
    strengths: tuple[BoundStrength, ...]
    photoionisation: tuple[BoundPhotoionisation, ...]
    assessment: str
    assessment_reason: str
    input_digests: Mapping[str, str] = field(default_factory=dict, repr=False)

    IDENTITY_KIND: ClassVar[str] = "two_photon_bound_levels"
    EVALUATION_VERSION: ClassVar[int] = 1

    def __post_init__(self) -> None:
        """Take ownership, and check every assigned value names a listed level.

        Raises
        ------
        DataError
            If an assigned strength or cross section names no listed level, or the
            assessment is unknown.
        """
        super().__post_init__()
        keys = {level.key for level in self.levels}
        named = [s.level for s in (*self.strengths, *self.photoionisation) if s.level]
        unknown = sorted(set(named) - keys)
        if unknown or self.assessment not in {"checked", "unresolved"}:
            msg = f"{self.gas.value}: unknown levels {unknown} or assessment"
            raise DataError(msg)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the record's name, gas and source records.

        Returns
        -------
        tuple
            ``("two_photon_bound/<gas>", gas, sources)``.
        """
        keys = [s.dataset for s in (*self.strengths, *self.photoionisation)]
        return (
            f"two_photon_bound/{self.gas.value}",
            self.gas.value,
            source_ids("nist_asd_2026/", *(f"{k}/" for k in keys)),
        )

    @property
    def digest(self) -> str:
        """The record's digest, ``identity.digest``.

        Returns
        -------
        str
            64 hexadecimal digits.
        """
        return self.identity.digest

    def level(self, configuration: str, term: str, j: int) -> BoundLevel:
        """Return a listed level.

        Parameters
        ----------
        configuration : str
            NIST's configuration.
        term : str
            NIST's term.
        j : int
            The total angular momentum.

        Returns
        -------
        BoundLevel
            The level.

        Raises
        ------
        DataError
            If it is not listed.
        """
        key = f"{configuration}|{term}|{j}"
        for level in self.levels:
            if level.key == key:
                return level
        msg = f"{self.gas.value}: no listed level {key}"
        raise DataError(msg)

    def strengths_of(self, level: BoundLevel) -> tuple[BoundStrength, ...]:
        """Return the strengths assigned to a level.

        Parameters
        ----------
        level : BoundLevel
            A listed level.

        Returns
        -------
        tuple of BoundStrength
            Its strengths, calculated and measured.
        """
        return tuple(s for s in self.strengths if s.level == level.key)

    def photoionisation_of(self, level: BoundLevel) -> tuple[BoundPhotoionisation, ...]:
        """Return the photoionisation cross sections assigned to a level.

        Parameters
        ----------
        level : BoundLevel
            A listed level.

        Returns
        -------
        tuple of BoundPhotoionisation
            Its cross sections.
        """
        return tuple(p for p in self.photoionisation if p.level == level.key)


def _plain(value: object) -> Any:
    """Return a value as plain JSON values.

    Parameters
    ----------
    value : object
        A field's value.

    Returns
    -------
    object
        Dicts, lists, strings, numbers, booleans or None.
    """
    if is_dataclass(value) and not isinstance(value, type):
        value = {
            f.name: getattr(value, f.name)
            for f in fields(value)
            if f.name not in _STANDING
        }
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, Gas):
        return value.value
    return value


def _dumps(record: TwoPhotonBoundLevels) -> str:
    """Serialise a record as JSON text, as stored.

    Parameters
    ----------
    record : TwoPhotonBoundLevels
        The record.

    Returns
    -------
    str
        Indented JSON.
    """
    return json.dumps(_plain(record), indent=1, allow_nan=False) + "\n"


def _object(cls: type, raw: Any, where: str) -> Mapping[str, Any]:
    """Return a stored object's fields, checking they are exactly the class's.

    Parameters
    ----------
    cls : type
        The dataclass.
    raw : object
        The JSON object.
    where : str
        For error messages.

    Returns
    -------
    Mapping
        The fields.

    Raises
    ------
    DataError
        If it is not an object with exactly the class's fields.
    """
    names = {f.name for f in fields(cls)} - _STANDING
    if not isinstance(raw, Mapping) or set(raw) != names:
        msg = f"{where}: a {cls.__name__} needs exactly the fields {sorted(names)}"
        raise DataError(msg)
    return raw


def _loads(text: str, where: str) -> TwoPhotonBoundLevels:
    """Read a stored record from its JSON text.

    Parameters
    ----------
    text : str
        The JSON, as :func:`_dumps` writes it.
    where : str
        For error messages.

    Returns
    -------
    TwoPhotonBoundLevels
        The record, status ``validated``.

    Raises
    ------
    DataError
        If the text is not a record of the current derivation version.
    """
    d = _object(TwoPhotonBoundLevels, json.loads(text), where)
    if d["derivation_version"] != BOUND_DERIVATION_VERSION:
        msg = f"{where}: derivation version {d['derivation_version']}; re-derive it"
        raise DataError(msg)
    return TwoPhotonBoundLevels(
        gas=as_gas(d["gas"]),
        derivation_version=int(d["derivation_version"]),
        ionisation=float(d["ionisation"]),
        levels=tuple(BoundLevel(**_object(BoundLevel, x, where)) for x in d["levels"]),
        strengths=tuple(
            BoundStrength(**_object(BoundStrength, x, where)) for x in d["strengths"]
        ),
        photoionisation=tuple(
            BoundPhotoionisation(
                level=str(x["level"]),
                unassigned=str(x["unassigned"]),
                dataset=str(x["dataset"]),
                source=str(x["source"]),
                omega=tuple(None if w is None else float(w) for w in x["omega"]),
                sigma=tuple(float(v) for v in x["sigma"]),
                use=str(x["use"]),
            )
            for x in (
                _object(BoundPhotoionisation, x, where) for x in d["photoionisation"]
            )
        ),
        assessment=str(d["assessment"]),
        assessment_reason=str(d["assessment_reason"]),
        input_digests={str(k): str(v) for k, v in d["input_digests"].items()},
        status="validated",
    )


@functools.cache
def _stored(gas: Gas) -> TwoPhotonBoundLevels:
    """Return a stored record.

    Parameters
    ----------
    gas : Gas
        The species.

    Returns
    -------
    TwoPhotonBoundLevels
        The record of ``data/gas/models/two_photon_bound/<gas>.json``.

    Raises
    ------
    DataError
        If the file is missing, malformed or holds another gas.
    """
    path = _DIRECTORY / f"{gas.value}.json"
    if not path.is_file():
        msg = f"{path.name}: the stored bound-level record is missing; derive it"
        raise DataError(msg)
    record = _loads(path.read_text(encoding="utf-8"), path.name)
    if record.gas != gas:
        msg = f"{path.name} holds the record of {record.gas.value}"
        raise DataError(msg)
    return record


def two_photon_bound_levels(gas: Gas | str) -> TwoPhotonBoundLevels:
    """Return what is curated about a noble gas's bound two-photon levels.

    Parameters
    ----------
    gas : Gas or str
        He, Ne, Ar, Kr or Xe.

    Returns
    -------
    TwoPhotonBoundLevels
        The record.

    Raises
    ------
    UnavailableModelError
        For any other gas.

    Examples
    --------
    >>> two_photon_bound_levels("He").assessment
    'checked'
    """
    species = as_gas(gas)
    if species.value not in _GASES:
        msg = (
            f"{species.value}: bound two-photon levels are curated for the noble gases"
        )
        raise UnavailableModelError(msg)
    return _stored(species)

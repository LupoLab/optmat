# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The identity of a model: a digest of what it evaluates, and its status.

A propagation code that records which model produced its numbers needs a name that
changes whenever the numbers can. Every model, response and bundle of the public
surface therefore has an :attr:`~Identified.identity`, a :class:`ModelIdentity` whose
``digest`` is the SHA-256 of a canonical record of the model (OEP 1, section 4.8):

- the class's qualified name;
- the ``EVALUATION_VERSION`` of the class and of every class it inherits evaluation
  from, and of each shared evaluator module it evaluates through, raised whenever that
  code changes its output;
- every field: the constructor's parameters, including those that fix a state, and
  the values evaluation reads (a fit's terms, a spectrum's grid), not the files that
  produced them;
- extras a class declares: its validity where it is not a field, values it resolves
  once from another component, and the digests of the components it is bound to.

The converse is not promised: an edit that leaves every number unchanged may still
change the digest. The package version is recorded beside the digest, not in it, so a
digest survives a release that did not touch the model.

The record is serialised as JSON with sorted keys and no insignificant whitespace,
UTF-8, floats in Python's shortest round-trip form, non-finite floats as the strings
``"inf"``, ``"-inf"`` and ``"nan"``, arrays as lists, a complex number as ``[re, im]``
and an enumeration by its value.

**Overrides** (OEP 1, section 6.3). ``model.with_parameters(**changes)`` rebuilds a
model with some of its inputs changed: the class declares the inputs an override may
name (``INPUTS``) and refuses the fields derived from them (``DERIVED``). The result
has the status ``user``, records the changes in its identity's ``modified`` and the
curated model's digest in ``base``, and drops the base's stated uncertainty.

A model owns its storage: its arrays are copied at construction and made read-only,
and its mappings are copied into read-only mappings, so the digest stays true for the
object's life. That guards against accidents; it is not a security boundary.

This module is part of the core and imports no medium.
"""

import abc
import dataclasses
import enum
import functools
import hashlib
import importlib
import json
import math
import re
from collections.abc import Iterator, Mapping
from typing import Any, ClassVar, Final, Literal, get_args

import numpy as np

import optmat
from optmat.exceptions import DataError, OptmatError

__all__ = [
    "STATUSES",
    "FrozenMapping",
    "Identified",
    "ModelIdentity",
    "Status",
    "canonical",
    "canonical_json",
    "describe",
    "frozen",
    "parse_status",
    "sha256_digest",
    "source_ids",
]

type Status = Literal["validated", "provisional", "untested", "comparison", "user"]
"""How far a model is to be trusted.

- ``validated``: a default, with a stated uncertainty;
- ``provisional``: a known stopgap;
- ``untested``: no reference to test it against;
- ``comparison``: a literature form kept for comparison;
- ``user``: built from parameters a caller supplied, or changed by an override; a
  curated model's stated uncertainty does not carry over to it.
"""

STATUSES: Final[tuple[str, ...]] = get_args(Status.__value__)
_NEEDS_REASON: Final = frozenset({"provisional", "untested"})


@dataclasses.dataclass(frozen=True)
class ModelIdentity:
    """The identity of a model: what it is, a digest of its numbers, and its status.

    Attributes
    ----------
    name : str
        The model's name, such as ``"fit/Ar"`` or ``"rate/default/Ar@800nm"``.
    kind : str
        The quantity the model gives: ``"gas_refractive"``,
        ``"gas_hyperpolarizability"``, ``"ionisation_rate"``, ``"raman_response"``,
        ``"solid"`` or ``"kerr"``. It is not the model's own ``kind``, which on a
        dispersion model names its provenance (published, fitted, ...).
    material : str
        The material: a gas, a solid, a mixture such as ``"dry air"``, or a custom
        species' name.
    digest : str
        The SHA-256 of the model's canonical record, 64 hexadecimal digits.
    optmat_version : str
        The package version that built the model; informational, not in the digest.
    sources : tuple of str
        The source records the model rests on.
    status : str
        One of ``STATUSES``.
    status_reason : str
        Why the status is what it is, where that needs saying.
    modified : tuple of tuple of str
        The overrides applied, ``(name, value)`` sorted by name; empty for a curated
        model.
    base : str
        The digest of the curated model an override started from, or empty.
    reference_temperature : float or None
        The temperature a gas polarizability refers to, K, where its source declares
        one.
    infrared_temperature : float or None
        The temperature of a molecular fit's fixed infrared terms, K.
    """

    name: str
    kind: str
    material: str
    digest: str
    optmat_version: str
    sources: tuple[str, ...]
    status: str
    status_reason: str
    modified: tuple[tuple[str, str], ...] = ()
    base: str = ""
    reference_temperature: float | None = None
    infrared_temperature: float | None = None


def parse_status(
    entry: Mapping[str, Any], where: str, *, default: Status
) -> tuple[Status, str]:
    """Return the status and reason a curated data entry declares.

    Parameters
    ----------
    entry : Mapping
        The entry, with optional keys ``status`` and ``status_reason``.
    where : str
        The entry's location, for error messages.
    default : str
        The status of an entry that declares none.

    Returns
    -------
    tuple of str
        The status and its reason (empty where none is given).

    Raises
    ------
    DataError
        If the status is unknown, is ``user`` (which only a caller's own parameters
        earn), or is ``provisional`` or ``untested`` without a reason.

    Examples
    --------
    >>> parse_status({"status": "untested", "status_reason": "no data"}, "x",
    ...              default="validated")
    ('untested', 'no data')
    """
    status = entry.get("status", default)
    reason = str(entry.get("status_reason", ""))
    if status not in STATUSES or status == "user":
        msg = f"{where}: status {status!r} is not one of {STATUSES[:-1]}"
        raise DataError(msg)
    if status in _NEEDS_REASON and not reason:
        msg = f"{where}: a {status} status needs a status_reason"
        raise DataError(msg)
    return status, reason


class FrozenMapping[K, V](Mapping[K, V]):
    """A read-only mapping that owns a copy of its items.

    Unlike :class:`types.MappingProxyType` it can be pickled, so a model that holds one
    can be sent to another process.

    Parameters
    ----------
    items : Mapping or iterable of pairs, optional
        The items, copied.

    Examples
    --------
    >>> m = FrozenMapping({"a": 1})
    >>> m["a"], len(m), dict(m)
    (1, 1, {'a': 1})
    """

    __slots__ = ("_items",)

    def __init__(self, items: Mapping[K, V] | Any = ()) -> None:
        self._items: dict[K, V] = dict(items)

    def __getitem__(self, key: K) -> V:
        """Return the value of a key.

        Parameters
        ----------
        key : object
            The key.

        Returns
        -------
        object
            Its value.
        """
        return self._items[key]

    def __iter__(self) -> Iterator[K]:
        """Iterate over the keys.

        Returns
        -------
        iterator
            The keys, in insertion order.
        """
        return iter(self._items)

    def __len__(self) -> int:
        """Return the number of items.

        Returns
        -------
        int
            The number of items.
        """
        return len(self._items)

    def __repr__(self) -> str:
        """Return a representation like a dict's.

        Returns
        -------
        str
            ``FrozenMapping({...})``.
        """
        return f"FrozenMapping({self._items!r})"

    def __reduce__(self) -> tuple[type, tuple[dict[K, V]]]:
        """Pickle as the class and a plain dict of the items.

        Returns
        -------
        tuple
            The class and its constructor's argument.
        """
        return (type(self), (self._items,))


def frozen(value: Any) -> Any:
    """Return a read-only copy of a value a model is to own.

    Parameters
    ----------
    value : object
        A field's value.

    Returns
    -------
    object
        An array copied and marked read-only; a mapping copied into a
        :class:`FrozenMapping` of frozen values; a list or tuple as a tuple of frozen
        values (a named tuple keeps its type); anything else unchanged.

    Examples
    --------
    >>> a = np.zeros(2)
    >>> b = frozen(a)
    >>> b.flags.writeable, b is a
    (False, False)
    >>> frozen({"x": [1, 2]})
    FrozenMapping({'x': (1, 2)})
    """
    if isinstance(value, np.ndarray):
        copy = np.array(value, copy=True)
        copy.flags.writeable = False
        return copy
    if isinstance(value, Mapping):
        return FrozenMapping({k: frozen(v) for k, v in value.items()})
    if isinstance(value, list) or type(value) is tuple:
        return tuple(frozen(v) for v in value)
    return value


def _number(value: float) -> float | str:
    """Return a float, or its name where it is not finite.

    Parameters
    ----------
    value : float
        The number.

    Returns
    -------
    float or str
        The number, or ``"inf"``, ``"-inf"`` or ``"nan"``.
    """
    if math.isfinite(value):
        return value
    return "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")


def _array(value: np.ndarray) -> Any:
    """Return an array as canonical nested lists.

    Parameters
    ----------
    value : numpy.ndarray
        A real array.

    Returns
    -------
    list
        Nested lists of canonical values.
    """
    if value.dtype.kind in "fi" and np.all(np.isfinite(value)):
        return value.astype(float).tolist()
    return canonical(value.tolist())


def _scalar(value: Any) -> tuple[bool, Any]:
    """Return the canonical form of a scalar.

    Parameters
    ----------
    value : object
        Any value.

    Returns
    -------
    tuple
        ``(True, form)`` for a scalar, else ``(False, None)``.
    """
    match value:
        case enum.Enum():
            return True, canonical(value.value)
        case None | bool() | int() | str():
            return True, value
        case float():
            return True, _number(value)
        case complex():
            return True, [_number(value.real), _number(value.imag)]
        case np.generic():
            return True, canonical(value.item())
        case _:
            return False, None


def _record(value: Any) -> dict[str, Any]:
    """Return a dataclass instance as its type, declared version and fields.

    Parameters
    ----------
    value : object
        A dataclass instance.

    Returns
    -------
    dict
        Its canonical fields, ``__type__`` and, where its class declares one,
        ``__version__``.
    """
    record = {
        f.name: canonical(getattr(value, f.name)) for f in dataclasses.fields(value)
    }
    record["__type__"] = _qualname(type(value))
    declared = vars(type(value)).get("EVALUATION_VERSION")
    if declared is not None:
        record["__version__"] = declared
    return record


def _container(value: Any) -> tuple[bool, Any]:
    """Return the canonical form of an array or a collection.

    Parameters
    ----------
    value : object
        Any value.

    Returns
    -------
    tuple
        ``(True, form)`` for an array, mapping, sequence or set, else
        ``(False, None)``.
    """
    match value:
        case np.ndarray() if np.iscomplexobj(value):
            return True, _array(np.stack([value.real, value.imag], axis=-1))
        case np.ndarray():
            return True, _array(value)
        case Mapping():
            return True, {str(canonical(k)): canonical(v) for k, v in value.items()}
        case list() | tuple():
            return True, [canonical(v) for v in value]
        case set() | frozenset():
            return True, sorted((canonical(v) for v in value), key=_dumps)
        case _:
            return False, None


def canonical(value: Any) -> Any:
    """Return a value as the plain JSON types of the canonical record.

    Parameters
    ----------
    value : object
        A field's value.

    Returns
    -------
    object
        None, a bool, an int, a str, a finite float, or lists and dicts of them.

    Raises
    ------
    TypeError
        For a value with no canonical form.

    Notes
    -----
    A model with an identity enters as its type and digest, which is how a model
    records the components it is bound to. Any other dataclass enters as its type, the
    ``EVALUATION_VERSION`` its class declares, and its fields.

    Examples
    --------
    >>> canonical((1.5, float("inf"), 2 + 1j))
    [1.5, 'inf', [2.0, 1.0]]
    """
    is_scalar, form = _scalar(value)
    if is_scalar:
        return form
    is_container, form = _container(value)
    if is_container:
        return form
    if isinstance(value, Identified):
        return {"__type__": _qualname(type(value)), "digest": value.identity.digest}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _record(value)
    msg = f"no canonical form for a {type(value).__qualname__}"
    raise TypeError(msg)


def _dumps(value: Any) -> str:
    """Return canonical JSON text of plain JSON values.

    Parameters
    ----------
    value : object
        Plain JSON values.

    Returns
    -------
    str
        The JSON, keys sorted, with no insignificant whitespace.
    """
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def canonical_json(value: Any) -> bytes:
    """Return the canonical JSON of a value, UTF-8.

    Parameters
    ----------
    value : object
        Anything :func:`canonical` accepts.

    Returns
    -------
    bytes
        The JSON, keys sorted, with no insignificant whitespace.

    Examples
    --------
    >>> canonical_json({"b": 1, "a": [0.1, float("nan")]})
    b'{"a":[0.1,"nan"],"b":1}'
    """
    return _dumps(canonical(value)).encode()


def sha256_digest(value: Any) -> str:
    """Return the SHA-256 of a value's canonical JSON.

    Parameters
    ----------
    value : object
        Anything :func:`canonical` accepts.

    Returns
    -------
    str
        64 hexadecimal digits.

    Examples
    --------
    >>> sha256_digest([])[:12]
    '4f53cda18c2b'
    """
    return hashlib.sha256(canonical_json(value)).hexdigest()


_SOURCE_KEY: Final = re.compile(r"\b([a-z][a-z0-9_]*_\d{4}[a-z]?)/")


def source_ids(*texts: str) -> tuple[str, ...]:
    """Return the source records named by dataset keys in some texts.

    Parameters
    ----------
    *texts : str
        Texts that may contain dataset keys, ``<source>/<dataset>``.

    Returns
    -------
    tuple of str
        The sources, sorted, each once.

    Examples
    --------
    >>> source_ids("nist_asd_2026/ie_Ar_I + nist_asd_2026/level", "tong_2002/cl_Ar")
    ('nist_asd_2026', 'tong_2002')
    """
    return tuple(sorted({m for text in texts for m in _SOURCE_KEY.findall(text)}))


_NOT_RECORDED: Final = frozenset({"status", "status_reason", "modified", "base"})
"""The fields that describe a model's standing, not its numbers: not in the digest."""


def describe(value: Any) -> str:
    """Return a short description of an override's value, for ``modified``.

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    str
        Its canonical JSON if that is short, else its SHA-256 prefixed ``sha256:``.

    Examples
    --------
    >>> describe(2.4)
    '2.4'
    >>> describe(list(range(100)))[:7]
    'sha256:'
    """
    text = canonical_json(value).decode()
    return text if len(text) <= 80 else f"sha256:{sha256_digest(value)}"


def _qualname(cls: type) -> str:
    """Return a class's module and qualified name.

    Parameters
    ----------
    cls : type
        The class.

    Returns
    -------
    str
        ``<module>.<qualname>``.
    """
    return f"{cls.__module__}.{cls.__qualname__}"


@dataclasses.dataclass(frozen=True, eq=False)
class Identified(abc.ABC):
    """A model that owns its storage and has an identity.

    A subclass is a frozen dataclass. It declares ``IDENTITY_KIND``, its own
    ``EVALUATION_VERSION`` (raised whenever its output changes), and in
    ``EVALUATES_THROUGH`` the shared evaluator modules whose module-level
    ``EVALUATION_VERSION`` also governs its output. It implements
    ``_identity_label``, and ``_identity_extras`` where its record needs more
    than its fields. A subclass that defines ``__post_init__`` calls this class's.

    Attributes
    ----------
    status : str
        One of ``STATUSES``; ``"user"`` unless a curated builder says otherwise.
    status_reason : str
        Why, where that needs saying.
    modified : tuple of tuple of str
        The overrides applied, ``(name, value)`` sorted by name.
    base : str
        The digest of the curated model the overrides started from, or empty.
    """

    IDENTITY_KIND: ClassVar[str]
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ()
    INPUTS: ClassVar[tuple[str, ...]] = ()
    """The parameters :meth:`~optmat.identity.Identified.with_parameters` may change."""
    DERIVED: ClassVar[tuple[str, ...]] = ()
    """Fields derived from the inputs, which an override may not name."""

    status: Status = dataclasses.field(default="user", kw_only=True, repr=False)
    status_reason: str = dataclasses.field(default="", kw_only=True, repr=False)
    modified: tuple[tuple[str, str], ...] = dataclasses.field(
        default=(), kw_only=True, repr=False
    )
    base: str = dataclasses.field(default="", kw_only=True, repr=False)

    def __post_init__(self) -> None:
        """Take ownership of every field: copy its arrays and mappings, read-only."""
        for f in dataclasses.fields(self):
            object.__setattr__(self, f.name, frozen(getattr(self, f.name)))

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Restore a pickled model, its storage read-only again.

        Parameters
        ----------
        state : dict
            The pickled attributes.
        """
        names = {f.name for f in dataclasses.fields(self)}
        for key, value in state.items():
            object.__setattr__(self, key, frozen(value) if key in names else value)

    @abc.abstractmethod
    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, material and source records.

        Returns
        -------
        tuple
            ``(name, material, sources)``.
        """

    def _identity_extras(self) -> Mapping[str, Any]:
        """Return what the record needs beyond the fields.

        Returns
        -------
        Mapping
            Empty by default.
        """
        return {}

    def _identity_temperatures(self) -> tuple[float | None, float | None]:
        """Return the reference temperature and the infrared terms' temperature.

        Returns
        -------
        tuple
            Both None by default, K.
        """
        return None, None

    def identity_record(self) -> dict[str, Any]:
        """Return the canonical record the digest is taken of.

        Returns
        -------
        dict
            ``{class, versions, fields, extras}`` in plain JSON types. The status is
            not in it.
        """
        cls = type(self)
        versions = {
            _qualname(klass): vars(klass)["EVALUATION_VERSION"]
            for klass in cls.__mro__
            if "EVALUATION_VERSION" in vars(klass)
        }
        versions |= {
            name: importlib.import_module(name).EVALUATION_VERSION
            for name in cls.EVALUATES_THROUGH
        }
        return {
            "class": _qualname(cls),
            "versions": versions,
            "fields": canonical(
                {
                    f.name: getattr(self, f.name)
                    for f in dataclasses.fields(self)
                    if f.name not in _NOT_RECORDED
                }
            ),
            "extras": canonical(self._identity_extras()),
        }

    def _override[T: "Identified"](self: T, changes: Mapping[str, Any]) -> T:
        """Return the model rebuilt with some of its inputs changed.

        The default replaces the fields, which reruns the class's checks; a class
        whose fields are derived from its inputs rebuilds them instead.

        Parameters
        ----------
        changes : Mapping
            Input names to their new values, each one of :attr:`INPUTS`.

        Returns
        -------
        Identified
            The rebuilt model, its status not yet set.
        """
        return dataclasses.replace(self, **changes)

    def _stated_uncertainty(self) -> Mapping[str, Any]:
        """Return the fields that state the curated uncertainty, reset by an override.

        Returns
        -------
        Mapping
            Field names to the value that states none; empty by default.
        """
        return {}

    def with_parameters[T: "Identified"](self: T, **changes: Any) -> T:
        """Return the model rebuilt with some of its inputs changed (OEP 1, 6.3).

        Parameters
        ----------
        **changes
            New values of inputs the class declares in :attr:`INPUTS`.

        Returns
        -------
        Identified
            The rebuilt model: status ``user``, the changes in ``identity.modified``,
            the curated model's digest in ``identity.base``, and no stated
            uncertainty unless the changes give one.

        Raises
        ------
        OptmatError
            With no change, or a name that is derived or not an input.
        """
        name = self.identity.name
        if not changes:
            msg = f"{name}: with_parameters needs at least one change"
            raise OptmatError(msg)
        derived = sorted(set(changes) & set(self.DERIVED))
        if derived:
            msg = (
                f"{name}: {', '.join(derived)} derived from the inputs; override "
                f"the inputs instead ({', '.join(self.INPUTS)})"
            )
            raise OptmatError(msg)
        unknown = sorted(set(changes) - set(self.INPUTS))
        if unknown:
            inputs = ", ".join(self.INPUTS) or "none"
            msg = f"{name}: no input {', '.join(unknown)}; its inputs are {inputs}"
            raise OptmatError(msg)
        rebuilt = self._override(changes)
        reset = {
            k: v for k, v in self._stated_uncertainty().items() if k not in changes
        }
        recorded = dict(self.modified) | {k: describe(v) for k, v in changes.items()}
        return dataclasses.replace(
            rebuilt,
            **reset,
            status="user",
            status_reason=(
                f"{name} ({self.status}) with {', '.join(sorted(changes))} "
                "overridden; its stated uncertainty does not carry over"
            ),
            modified=tuple(sorted(recorded.items())),
            base=self.base or self.identity.digest,
        )

    @functools.cached_property
    def identity(self) -> ModelIdentity:
        """The model's identity, computed once.

        Returns
        -------
        ModelIdentity
            Its name, kind, material, digest, sources and status.
        """
        name, material, sources = self._identity_label()
        reference, infrared = self._identity_temperatures()
        return ModelIdentity(
            name=name,
            kind=type(self).IDENTITY_KIND,
            material=material,
            digest=sha256_digest(self.identity_record()),
            optmat_version=optmat.__version__,
            sources=tuple(sources),
            status=self.status,
            status_reason=self.status_reason,
            modified=self.modified,
            base=self.base,
            reference_temperature=reference,
            infrared_temperature=infrared,
        )

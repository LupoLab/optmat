# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The public surface of optmat for a propagation code (OEP 1, section 6.1).

    uv run python tools/public_api.py

writes ``tests/data/public_api.json``: for every name in the ``__all__`` of each module
of :data:`SURFACE_MODULES`, its signature, or for a class its bases, constructor,
dataclass fields, public methods and properties, and enum members.
``tests/test_public_api.py`` compares the surface with it. A change to the surface
fails until the file is re-recorded, deliberately, in the same commit, with a
changelog line. Before 1.0 a minor version may break the surface and says so; a
patch version never does.
"""

import annotationlib
import dataclasses
import enum
import functools
import importlib
import inspect
import json
from pathlib import Path, PurePath
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tests" / "data" / "public_api.json"
SURFACE_MODULES = (
    "optmat",
    "optmat.conditions",
    "optmat.exceptions",
    "optmat.identity",
    "optmat.gas.absorption",
    "optmat.gas.air",
    "optmat.gas.comparison",
    "optmat.gas.composition",
    "optmat.gas.eos",
    "optmat.gas.gases",
    "optmat.gas.ionisation",
    "optmat.gas.ionisation_rates",
    "optmat.gas.mixture",
    "optmat.gas.mixture_eos",
    "optmat.gas.mixture_validation",
    "optmat.gas.models.base",
    "optmat.gas.models.compact",
    "optmat.gas.models.fitted",
    "optmat.gas.models.formulas",
    "optmat.gas.models.oscillator",
    "optmat.gas.models.registry",
    "optmat.gas.molecular_ionisation_rates",
    "optmat.gas.nonlinear",
    "optmat.gas.raman_response",
    "optmat.gas.transport",
    "optmat.solid.bundle",
    "optmat.solid.kerr",
    "optmat.solid.models.base",
    "optmat.solid.models.fitted",
    "optmat.solid.models.formulas",
    "optmat.solid.models.registry",
    "optmat.solid.models.tabulated",
    "optmat.solid.raman_response",
    "optmat.solid.udm",
)


def _signature(obj: Any) -> str:
    """Return a callable's signature, its annotations as written in the source."""
    try:
        # Python 3.14's API: this tool runs on 3.14 (the type check targets 3.12).
        return str(
            inspect.signature(obj, annotation_format=annotationlib.Format.STRING)  # pyright: ignore[reportCallIssue]
        )
    except (TypeError, ValueError):  # a builtin without a signature
        return "(...)"


def _annotations(cls: type) -> dict[str, str]:
    """Return the annotation of every field along a class's MRO, as a string."""
    out: dict[str, str] = {}
    for klass in reversed(cls.__mro__):
        try:
            notes = annotationlib.get_annotations(
                klass, format=annotationlib.Format.STRING
            )
        except TypeError:  # a builtin base
            continue
        out.update({k: str(v) for k, v in notes.items()})
    return out


def _default(field: dataclasses.Field[Any]) -> str | None:
    """Return a field's default, ``<factory>``, or None for a required field."""
    if field.default is not dataclasses.MISSING:
        return repr(field.default)
    if field.default_factory is not dataclasses.MISSING:
        return "<factory>"
    return None


def _class(cls: type) -> dict[str, Any]:
    """Return a class's public shape."""
    entry: dict[str, Any] = {
        "kind": "class",
        "bases": [f"{b.__module__}.{b.__qualname__}" for b in cls.__bases__],
    }
    if issubclass(cls, enum.Enum):
        entry["members"] = {m.name: repr(m.value) for m in cls}
        return entry
    entry["init"] = _signature(cls)
    if dataclasses.is_dataclass(cls):
        notes = _annotations(cls)
        entry["fields"] = [
            {
                "name": f.name,
                "type": notes.get(f.name, ""),
                "default": _default(f),
                "init": f.init,
                "kw_only": f.kw_only is True,
            }
            for f in dataclasses.fields(cls)
        ]
    members: dict[str, str] = {}
    for name, value in sorted(vars(cls).items()):
        if name.startswith("_"):
            continue
        if isinstance(value, property | functools.cached_property):
            members[name] = "property"
        elif isinstance(value, staticmethod | classmethod):
            members[name] = f"{type(value).__name__}{_signature(value.__func__)}"
        elif callable(value):
            members[name] = f"method{_signature(value)}"
        elif not dataclasses.is_dataclass(cls) or name not in {
            f.name for f in dataclasses.fields(cls)
        }:
            members[name] = f"attribute: {type(value).__name__}"
    entry["members"] = members
    return entry


def _describe(value: Any) -> dict[str, Any]:
    """Return the public shape of one name of a module's ``__all__``."""
    if inspect.isclass(value):
        return _class(value)
    if callable(value):
        return {"kind": "function", "signature": _signature(value)}
    # A path is a PosixPath or a WindowsPath by platform; the surface is the same.
    kind = "Path" if isinstance(value, PurePath) else type(value).__name__
    return {"kind": "constant", "type": kind}


def snapshot() -> dict[str, dict[str, Any]]:
    """Return the public surface.

    Returns
    -------
    dict
        ``{module: {name: shape}}`` for every module of :data:`SURFACE_MODULES`.
    """
    out: dict[str, dict[str, Any]] = {}
    for name in SURFACE_MODULES:
        module = importlib.import_module(name)
        out[name] = {n: _describe(getattr(module, n)) for n in sorted(module.__all__)}
    return out


def main() -> None:
    """Write the snapshot."""
    OUTPUT.write_text(
        json.dumps(snapshot(), indent=1, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

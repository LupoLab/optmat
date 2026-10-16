# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The materials the curated records name.

Every dataset belongs to one material: a gas, a glass or a crystal. The core knows only
what a material *is* (its id, name, kind, formula and optical symmetry), read from
``data/materials.toml``; what a medium needs beyond that (a gas's electron count and
equation of state, a solid's models) lives in that medium's subpackage (TC-318).

A material's id is the directory of its extracted files,
``data/values/<family>/<id>/``. An anisotropic crystal lists the axes a dataset's
``parameters.axis`` may take (``o`` and ``e`` for a uniaxial crystal).
"""

import functools
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Literal

from optmat.exceptions import OptmatError

__all__ = [
    "KINDS",
    "MATERIALS_FILE",
    "Material",
    "as_material",
    "material",
    "materials",
]

MATERIALS_FILE: Final = Path(__file__).parent / "data" / "materials.toml"
"""The registry of materials."""

KINDS: Final[tuple[str, ...]] = ("gas", "ion", "glass", "crystal", "metal")
"""The kinds of material."""

_AXES: Final[dict[str, int]] = {"isotropic": 0, "uniaxial": 2, "biaxial": 3}


@dataclass(frozen=True)
class Material:
    """What a material is.

    Attributes
    ----------
    id : str
        The id records use, and the directory of the material's files.
    name : str
        A descriptive name.
    kind : {"gas", "ion", "glass", "crystal", "metal"}
        The kind of medium; an ``ion`` is a gas-phase ion, whose data (its
        polarizability, for the ionisation rates) never enter a gas's models. A
        ``metal`` is a solid whose index is complex throughout (TC-543).
    formula : str
        The chemical formula, or an empty string (a multicomponent glass).
    optics : {"isotropic", "uniaxial", "biaxial"}
        The optical symmetry.
    axes : tuple of str
        The axes a dataset's ``parameters.axis`` may name; empty when isotropic.
    aliases : tuple of str
        Other names that :func:`as_material` accepts.
    """

    id: str
    name: str
    kind: Literal["gas", "ion", "glass", "crystal", "metal"]
    formula: str
    optics: Literal["isotropic", "uniaxial", "biaxial"]
    axes: tuple[str, ...]
    aliases: tuple[str, ...]

    def __post_init__(self) -> None:
        """Check the kind, the symmetry and the axes.

        Raises
        ------
        OptmatError
            If the kind or symmetry is unknown, or the axes do not match the symmetry.
        """
        if self.kind not in KINDS:
            msg = f"material {self.id}: unknown kind {self.kind!r}"
            raise OptmatError(msg)
        if self.optics not in _AXES:
            msg = f"material {self.id}: unknown optics {self.optics!r}"
            raise OptmatError(msg)
        if len(self.axes) != _AXES[self.optics]:
            msg = f"material {self.id}: a {self.optics} material has "
            msg += f"{_AXES[self.optics]} axes, not {len(self.axes)}"
            raise OptmatError(msg)

    @property
    def is_solid(self) -> bool:
        """Whether the material is a glass, a crystal or a metal.

        Returns
        -------
        bool
            True for a glass, a crystal or a metal.
        """
        return self.kind in {"glass", "crystal", "metal"}


def _material(material_id: str, entry: Mapping[str, object]) -> Material:
    """Build one material from its registry entry.

    Parameters
    ----------
    material_id : str
        The id (the entry's table name).
    entry : mapping
        The entry.

    Returns
    -------
    Material
        The material.

    Raises
    ------
    OptmatError
        If the entry has unknown keys.
    """
    unknown = set(entry) - {"name", "kind", "formula", "optics", "axes", "aliases"}
    if unknown:
        msg = f"material {material_id}: unknown keys {sorted(unknown)}"
        raise OptmatError(msg)
    axes = entry.get("axes", [])
    aliases = entry.get("aliases", [])
    return Material(
        id=material_id,
        name=str(entry["name"]),
        kind=str(entry["kind"]),  # pyright: ignore[reportArgumentType] - checked
        formula=str(
            entry.get("formula", material_id if entry["kind"] in {"gas", "ion"} else "")
        ),
        optics=str(entry.get("optics", "isotropic")),  # pyright: ignore[reportArgumentType]
        axes=tuple(str(a) for a in axes),  # pyright: ignore[reportGeneralTypeIssues]
        aliases=tuple(str(a) for a in aliases),  # pyright: ignore[reportGeneralTypeIssues]
    )


@functools.cache
def materials() -> dict[str, Material]:
    """Return every registered material, by id.

    Returns
    -------
    dict of str to Material
        The materials, in registry order.

    Examples
    --------
    >>> materials()["SiO2"].kind
    'glass'
    """
    with MATERIALS_FILE.open("rb") as handle:
        registry = tomllib.load(handle)
    return {key: _material(key, value) for key, value in registry.items()}


@functools.cache
def _names() -> dict[str, str]:
    """Map every id and alias to its id.

    Returns
    -------
    dict of str to str
        Names to ids.

    Raises
    ------
    OptmatError
        If a name is claimed twice.
    """
    names: dict[str, str] = {}
    for entry in materials().values():
        for name in (entry.id, *entry.aliases):
            if names.setdefault(name, entry.id) != entry.id:
                msg = (
                    f"material name {name!r} is claimed by {names[name]} and {entry.id}"
                )
                raise OptmatError(msg)
    return names


def as_material(name: str) -> str:
    """Return the id of the material ``name`` names.

    Parameters
    ----------
    name : str
        An id or an alias, e.g. ``"SiO2"`` or ``"fused silica"``.

    Returns
    -------
    str
        The material's id.

    Raises
    ------
    OptmatError
        If ``name`` names no registered material.

    Examples
    --------
    >>> as_material("fused silica")
    'SiO2'
    """
    try:
        return _names()[str(name)]
    except KeyError as exc:
        msg = (
            f"unknown material {name!r}; it must be registered in {MATERIALS_FILE.name}"
        )
        raise OptmatError(msg) from exc


def material(name: str) -> Material:
    """Return the registered material ``name`` names.

    Parameters
    ----------
    name : str
        An id or an alias.

    Returns
    -------
    Material
        The material.

    Examples
    --------
    >>> material("sapphire").axes
    ('o', 'e')
    """
    return materials()[as_material(name)]

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Curated, validated optical properties of materials.

optmat holds measured and computed optical properties with full provenance: each data
point keeps its source, the conditions it was measured under and the technique that
produced it, and every model states its reference state, validity range and citation.

The package has a core and one subpackage per kind of medium:

- the core: the curated records, their vocabulary and their materials
  (:mod:`optmat.datasets`, :mod:`optmat.vocabulary`, :mod:`optmat.materials`), states
  and wavelengths (:mod:`optmat.conditions`), model identities (:mod:`optmat.identity`),
  constants and exceptions;
- :mod:`optmat.gas`: the gases (linear index, absorption, third-order response, Raman
  response and ionisation);
- :mod:`optmat.solid`: glasses, crystals and metals (their dielectric functions, the
  Kerr and Raman response of fused silica, and the data of their nonlinear response,
  absorption and damage).

Every submodule is an attribute of its package, loaded on first access (:pep:`562`),
so ``import optmat`` does not load CoolProp.
"""

import importlib
from importlib.metadata import PackageNotFoundError, version
from types import ModuleType
from typing import Final

from optmat.exceptions import OptmatError

__all__ = ["OptmatError", "__version__"]

_SUBMODULES: Final = (
    "conditions",
    "constants",
    "datasets",
    "exceptions",
    "gas",
    "identity",
    "materials",
    "solid",
    "vocabulary",
)

try:
    __version__: str = version("optmat")
except PackageNotFoundError:  # pragma: no cover - only in an uninstalled source tree
    __version__ = "0.0.0"


def __getattr__(name: str) -> ModuleType:
    """Import a submodule on first access.

    Parameters
    ----------
    name : str
        The attribute asked for.

    Returns
    -------
    types.ModuleType
        The submodule of that name.

    Raises
    ------
    AttributeError
        If ``name`` is not a submodule.
    """
    if name in _SUBMODULES:
        return importlib.import_module(f"{__name__}.{name}")
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)


def __dir__() -> list[str]:
    """List the public names and the submodules.

    Returns
    -------
    list of str
        The sorted public names and submodule names.
    """
    return sorted({*__all__, *_SUBMODULES})

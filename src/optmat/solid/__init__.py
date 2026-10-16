# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Solids: glasses, crystals and metals (TC-318).

The dielectric-function models (:mod:`optmat.solid.models`: the interface, the
published formulas, the metals' tables and the registry of recommended models), Franta
et al.'s universal dispersion model (:mod:`optmat.solid.udm`), and fused silica's Kerr
coefficient and Raman response (:mod:`optmat.solid.kerr`,
:mod:`optmat.solid.raman_response`). The solid models' configurations and stored
fits are in :data:`SOLID_DATA_DIRECTORY`.
"""

import importlib
from pathlib import Path
from types import ModuleType
from typing import Final

__all__ = ["SOLID_DATA_DIRECTORY"]

SOLID_DATA_DIRECTORY: Final = Path(__file__).parents[1] / "data" / "solid"
"""The solid models' configurations and stored fits (the curated records are shared)."""

_SUBMODULES: Final = (
    "bundle",
    "fitting",
    "kerr",
    "models",
    "raman_response",
    "raman_response_validation",
    "udm",
)


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

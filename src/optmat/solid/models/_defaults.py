# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The recommended model of each solid, and the status it gives a model.

Kept apart from the registry so that the model builders can give each model its
status without importing the registry, which imports them.
"""

import dataclasses
import functools
import tomllib
from collections.abc import Mapping
from typing import Final

from optmat.materials import as_material
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models.base import Model

DEFAULTS_FILE: Final = SOLID_DATA_DIRECTORY / "defaults.toml"
"""The recommended model of each solid, with the reason."""


@functools.cache
def defaults() -> dict[str, tuple[dict[str, str], str]]:
    """Return the recommended model names and the reason of each solid.

    Returns
    -------
    dict of str to tuple
        ``({axis: model name}, reason)`` by material id; the axis is ``""`` for an
        isotropic material.
    """
    with DEFAULTS_FILE.open("rb") as handle:
        record = tomllib.load(handle)
    out = {}
    for key, entry in record.items():
        chosen = entry["model"]
        names = (
            {str(k): str(v) for k, v in chosen.items()}
            if isinstance(chosen, dict)
            else {"": str(chosen)}
        )
        out[as_material(key)] = (names, str(entry["reason"]))
    return out


def curated[M: Model](models: Mapping[str, M]) -> dict[str, M]:
    """Return curated models with their statuses.

    Parameters
    ----------
    models : Mapping of str to Model
        The models, by name.

    Returns
    -------
    dict of str to Model
        Each model ``validated`` where ``defaults.toml`` recommends it (for any axis),
        with its reason, and ``comparison`` otherwise.
    """
    out = {}
    for name, model in models.items():
        names, reason = defaults().get(model.material, ({}, ""))
        if name in names.values():
            status, why = "validated", reason
        else:
            status = "comparison"
            why = f"kept for comparison; the default of {model.material} is " + (
                ", ".join(sorted(set(names.values()))) or "none yet"
            )
        out[name] = dataclasses.replace(model, status=status, status_reason=why)
    return out

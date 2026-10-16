# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every model of a solid's dielectric function, and the recommended one.

Models come from the fits (:mod:`optmat.solid.models.fitted`), the published
formulas of the source records (:mod:`optmat.solid.models.formulas`) and the stitched
tables of curated optical constants (:mod:`optmat.solid.models.tabulated`). The
recommended model of each material is chosen, with its reason, in
``data/solid/defaults.toml``.
"""

import functools
from collections.abc import Callable, Mapping
from typing import Final

from optmat.exceptions import OptmatError
from optmat.materials import as_material
from optmat.solid.models._defaults import DEFAULTS_FILE
from optmat.solid.models._defaults import defaults as _defaults
from optmat.solid.models.base import Model
from optmat.solid.models.fitted import fitted_models
from optmat.solid.models.formulas import published_model, published_models
from optmat.solid.models.tabulated import tabulated_models

__all__ = ["DEFAULTS_FILE", "default_model", "default_reason", "model", "models"]


_REGISTRIES: Final[dict[str, Callable[[], Mapping[str, Model]]]] = {
    "fit": fitted_models,
    "table": tabulated_models,
}
"""The registry that holds each prefix of a name; any other name is a published
formula's (``<source>/<formula>``)."""


@functools.cache
def _all_models() -> dict[str, Model]:
    """Return every model by name.

    Returns
    -------
    dict of str to Model
        The fitted, published and tabulated models.
    """
    return {**fitted_models(), **published_models(), **tabulated_models()}


def models(material: str) -> dict[str, Model]:
    """Return every model of a solid, by name.

    Parameters
    ----------
    material : str
        The material's id or an alias.

    Returns
    -------
    dict of str to Model
        The models, named ``<source>/<formula>``.
    """
    wanted = as_material(material)
    return {name: m for name, m in _all_models().items() if m.material == wanted}


def model(name: str) -> Model:
    """Return a model by name.

    Parameters
    ----------
    name : str
        ``<source>/<formula>``.

    Returns
    -------
    Model
        The model.

    Raises
    ------
    OptmatError
        If no model has that name.
    """
    # Only the registry that holds the name is built (as for the gases, TC-810), and
    # a published formula reads only its own source's record.
    registry = _REGISTRIES.get(name.split("/", 1)[0])
    try:
        return registry()[name] if registry else published_model(name)
    except KeyError as exc:
        msg = f"no solid model named {name!r}"
        raise OptmatError(msg) from exc


def default_model(material: str, axis: str | None = None) -> Model:
    """Return the recommended model of a solid.

    Parameters
    ----------
    material : str
        The material's id or an alias.
    axis : str, optional
        For an anisotropic crystal, the axis (``o`` or ``e`` for a uniaxial one);
        required when the recommendation is per axis.

    Returns
    -------
    Model
        The model ``data/solid/defaults.toml`` names for it.

    Raises
    ------
    OptmatError
        If no model is recommended for the material (or the axis) yet, or an axis is
        needed and not given.
    """
    wanted = as_material(material)
    choice = _defaults().get(wanted)
    if choice is None:
        msg = f"no dispersion model is recommended for {wanted} yet"
        raise OptmatError(msg)
    names = choice[0]
    if "" in names:
        return model(names[""])
    if axis is None:
        msg = f"{wanted} is anisotropic: give the axis ({', '.join(sorted(names))})"
        raise OptmatError(msg)
    if axis not in names:
        msg = f"no dispersion model is recommended for {wanted} axis {axis!r} yet"
        raise OptmatError(msg)
    return model(names[axis])


def default_reason(material: str) -> str:
    """Return why the recommended model of a solid was chosen.

    Parameters
    ----------
    material : str
        The material's id or an alias.

    Returns
    -------
    str
        The reason ``data/solid/defaults.toml`` gives.

    Raises
    ------
    OptmatError
        If no model is recommended for the material yet.
    """
    wanted = as_material(material)
    choice = _defaults().get(wanted)
    if choice is None:
        msg = f"no dispersion model is recommended for {wanted} yet"
        raise OptmatError(msg)
    return choice[1]

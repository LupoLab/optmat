# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every dispersion model of a gas, and the recommended one.

Models come from four places: the fits (``models.fitted``), their compact forms
(``models.compact``, OEP 2), the published formulas of the source records
(``models.formulas``) and the refractiveindex.info comparison entries (``comparison``).
The recommended model of each gas is chosen, with its reason, in
``data/gas/defaults.toml``.
"""

import functools
from collections.abc import Callable, Mapping
from typing import Final

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat.conditions import ONE_ATM_0C, State
from optmat.exceptions import OptmatError
from optmat.gas.comparison import database_models
from optmat.gas.gases import Gas, as_gas
from optmat.gas.models._defaults import DEFAULTS_FILE
from optmat.gas.models._defaults import defaults as _defaults
from optmat.gas.models.base import Model
from optmat.gas.models.compact import compact_models
from optmat.gas.models.fitted import fitted_models
from optmat.gas.models.formulas import published_models

__all__ = [
    "DEFAULTS_FILE",
    "default_model",
    "model",
    "models",
    "refractive_index",
    "refractivity",
]


_REGISTRIES: Final[dict[str, Callable[[], Mapping[str, Model]]]] = {
    "fit": fitted_models,
    "compact": compact_models,
    "refractiveindex.info": database_models,
}
"""The registry that holds each prefix of a name; any other name is a published
formula's (``<source>/<formula>``)."""


@functools.cache
def _all_models() -> dict[str, Model]:
    """Return every model by name.

    Returns
    -------
    dict of str to Model
        Fitted, compact, published and comparison models.
    """
    return {
        **fitted_models(),
        **compact_models(),
        **published_models(),
        **database_models(),
    }


def models(gas: Gas | str) -> dict[str, Model]:
    """Return every model of a gas, by name.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    dict of str to Model
        The fitted models (``fit/<gas>``), the compact models of the molecules
        (``compact/<gas>``), the published formulas (``<source>/<formula>``) and the
        refractiveindex.info entries.

    Examples
    --------
    >>> "fit/Ar" in models("Ar")
    True
    """
    wanted = as_gas(gas)
    return {name: m for name, m in _all_models().items() if m.gas is wanted}


def model(name: str) -> Model:
    """Return a model by name.

    Parameters
    ----------
    name : str
        ``fit/<gas>``, ``compact/<gas>``, ``<source>/<formula>`` or
        ``refractiveindex.info/<gas>/<entry>``.

    Returns
    -------
    Model
        The model.

    Raises
    ------
    OptmatError
        If no model has that name.
    """
    # Only the registry that holds the name is built: the published formulas read the
    # whole source catalogue, which takes seconds, and a fit's lookup does not need it.
    registry = _REGISTRIES.get(name.split("/", 1)[0], published_models)
    try:
        return registry()[name]
    except KeyError as exc:
        msg = f"no model named {name!r}"
        raise OptmatError(msg) from exc


def default_model(gas: Gas | str) -> Model:
    """Return the recommended model of a gas.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    Model
        The model ``data/gas/defaults.toml`` names for it.

    Raises
    ------
    OptmatError
        If no model is recommended for the gas (TC-733).

    Examples
    --------
    >>> default_model("He").name
    'fit/He'
    """
    choice = _defaults().get(as_gas(gas))
    if choice is None:
        msg = f"no dispersion model is recommended for {as_gas(gas).value} yet"
        raise OptmatError(msg)
    return model(choice[0])


def refractivity(
    gas: Gas | str,
    wavelength: ArrayLike | float,
    state: State = ONE_ATM_0C,
    *,
    model_name: str | None = None,
) -> Float[np.ndarray, "*shape"]:
    """Return ``n - 1`` of a gas: the one-call entry point.

    Parameters
    ----------
    gas : Gas or str
        The species.
    wavelength : array_like or float
        Vacuum wavelength, m.
    state : State, optional
        Pressure and temperature (real-gas density, Lorentz-Lorenz with the virial
        term); 0 degrees C and 1 atm by default.
    model_name : str, optional
        A model by name; the recommended one by default.

    Returns
    -------
    numpy.ndarray
        ``n - 1``, the shape of ``wavelength``.

    Warns
    -----
    ModelRangeWarning
        Outside the model's range.

    Examples
    --------
    >>> from optmat.conditions import State
    >>> f"{float(refractivity('Ar', 800e-9, State(101325.0, 293.15))) * 1e7:.2f}"
    '2605.28'
    """
    chosen = default_model(gas) if model_name is None else model(model_name)
    return chosen.refractivity(wavelength, state)


def refractive_index(
    gas: Gas | str,
    wavelength: ArrayLike | float,
    state: State = ONE_ATM_0C,
    *,
    model_name: str | None = None,
) -> Complex[np.ndarray, "*shape"]:
    """Return the complex refractive index of a gas: the entry point for absorption.

    Parameters
    ----------
    gas : Gas or str
        The species.
    wavelength : array_like or float
        Vacuum wavelength, m.
    state : State, optional
        Pressure and temperature (real-gas density, Lorentz-Lorenz with the virial
        term); 0 degrees C and 1 atm by default.
    model_name : str, optional
        A model by name; the recommended one by default.

    Returns
    -------
    numpy.ndarray
        ``n`` (complex; ``Im n > 0`` is absorption, zero for a transparent model), the
        shape of ``wavelength``.

    Warns
    -----
    ModelRangeWarning
        Outside the model's range.

    Examples
    --------
    >>> n = refractive_index("Ar", 800e-9)
    >>> f"{float(n.real - 1) * 1e7:.2f} {float(n.imag):.1f}"
    '2796.77 0.0'
    """
    chosen = default_model(gas) if model_name is None else model(model_name)
    return chosen.refractive_index(wavelength, state)

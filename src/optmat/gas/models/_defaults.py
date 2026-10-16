# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The recommended dispersion model of each gas, and the status it gives a model.

Kept apart from the registry so that the model builders can give each model its
status without importing the registry, which imports them.
"""

import functools
import tomllib
from typing import Final

from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.identity import Status

DEFAULTS_FILE: Final = GAS_DATA_DIRECTORY / "defaults.toml"
"""The recommended model of each gas, with the reason."""


@functools.cache
def defaults() -> dict[Gas, tuple[str, str]]:
    """Return the recommended model name and reason of each gas.

    Returns
    -------
    dict of Gas to tuple of str
        ``(model name, reason)``.
    """
    with DEFAULTS_FILE.open("rb") as handle:
        record = tomllib.load(handle)
    return {as_gas(g): (entry["model"], entry["reason"]) for g, entry in record.items()}


def curated_status(gas: Gas, name: str) -> tuple[Status, str]:
    """Return the status of a curated dispersion model.

    Parameters
    ----------
    gas : Gas
        The species.
    name : str
        The model's name.

    Returns
    -------
    tuple of str
        ``("validated", reason)`` for the model ``defaults.toml`` recommends, with its
        reason; ``("comparison", ...)`` for any other curated model, which is kept to
        compare with it.
    """
    choice = defaults().get(gas)
    if choice is not None and choice[0] == name:
        return "validated", choice[1]
    default = choice[0] if choice is not None else "none yet"
    return "comparison", f"kept for comparison; the default of {gas.value} is {default}"

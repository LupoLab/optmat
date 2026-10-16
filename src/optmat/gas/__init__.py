# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The gases: curated data and models of their optical response.

The linear index (the per-particle polarizability, :mod:`optmat.gas.models`, fitted by
:mod:`optmat.gas.fitting` and reduced from the data by :mod:`optmat.gas.reduction` with
the densities of :mod:`optmat.gas.eos` and Lorentz-Lorenz), the absorption and the
extreme ultraviolet (:mod:`optmat.gas.absorption`, :mod:`optmat.gas.euv`), the
infrared bands (:mod:`optmat.gas.infrared`), the second hyperpolarizability
(:mod:`optmat.gas.nonlinear`), the Raman lines, their widths and the response functions
(:mod:`optmat.gas.raman`, :mod:`optmat.gas.linewidths`,
:mod:`optmat.gas.raman_response`), the ionisation parameters and the strong-field
ionisation rates of atoms and molecules (:mod:`optmat.gas.ionisation`,
:mod:`optmat.gas.ionisation_rates`, :mod:`optmat.gas.molecular_ionisation_rates`) and
the composition of dry air (:mod:`optmat.gas.composition`). The gas models'
configurations and stored fits are in :data:`GAS_DATA_DIRECTORY`.
"""

import importlib
from pathlib import Path
from types import ModuleType
from typing import Final

__all__ = ["GAS_DATA_DIRECTORY"]

# Not built from optmat.datasets.DATA_DIRECTORY: the dataset loader imports the gas
# species, and so this package, first.
GAS_DATA_DIRECTORY: Final = Path(__file__).parents[1] / "data" / "gas"
"""The gas models' configurations and stored fits (the curated records are shared)."""

_SUBMODULES: Final = (
    "absorption",
    "air",
    "comparison",
    "composition",
    "crosschecks",
    "eos",
    "euv",
    "fitting",
    "gases",
    "infrared",
    "ionisation",
    "ionisation_rates",
    "ionisation_validation",
    "linewidths",
    "lorentz_lorenz",
    "mixture",
    "mixture_eos",
    "mixture_validation",
    "models",
    "molecular_ionisation_rates",
    "molecular_ionisation_validation",
    "nonlinear",
    "raman",
    "raman_response",
    "raman_response_validation",
    "reduction",
    "transport",
    "two_photon",
    "two_photon_bound",
    "virial",
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

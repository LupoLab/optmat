# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""From a gas dataset's printed values to the per-particle polarizability.

A refractivity measured at a stated pressure and temperature, a polarizability in some
unit, a molar refractivity or a refractivity per density each imply the per-particle
dipole polarizability at their frequency. The conversion needs the number density the
source meant (real-gas unless the source reduced with the ideal gas,
:mod:`optmat.gas.eos`) and, for a refractivity, Lorentz-Lorenz with the gas's
refractivity virial term (:mod:`optmat.gas.lorentz_lorenz`, :mod:`optmat.gas.virial`).
"""

import math
from typing import Final

import numpy as np
from jaxtyping import Float

from optmat.conditions import ONE_ATM_0C, State
from optmat.constants import AU_POLARIZABILITY, EPS0, N_A
from optmat.datasets import Dataset
from optmat.exceptions import DataError
from optmat.gas.eos import number_density
from optmat.gas.lorentz_lorenz import polarizability_from_refractivity, refractivity

__all__ = [
    "POLARIZABILITY_UNITS",
    "density_at_state",
    "polarizability",
    "polarizability_uncertainty",
    "refractivity_at",
]

POLARIZABILITY_UNITS: Final[dict[str, float]] = {
    "au": AU_POLARIZABILITY,
    # A polarizability *volume* alpha / (4 pi eps0), in A^3 or cm^3.
    "angstrom3": 4 * math.pi * EPS0 * 1e-30,
    "cm3": 4 * math.pi * EPS0 * 1e-6,
    "SI": 1.0,
}
"""Factors from each polarizability unit of the vocabulary to C m^2/V."""


def density_at_state(dataset: Dataset) -> float:
    """Return the number density the source meant at the dataset's state.

    Parameters
    ----------
    dataset : Dataset
        A dataset with a state.

    Returns
    -------
    float
        Number density, m^-3: real-gas unless the source reduced with the ideal gas.

    Raises
    ------
    DataError
        If the dataset has no state.
    """
    if dataset.state is None:
        msg = f"{dataset.key}: no state"
        raise DataError(msg)
    return number_density(
        dataset.material, dataset.state, ideal=dataset.density == "ideal"
    )


def _virial(dataset: Dataset) -> float:
    """Return the gas's refractivity virial ratio ``B_R / A_R``.

    Parameters
    ----------
    dataset : Dataset
        The dataset, for its gas.

    Returns
    -------
    float
        The ratio in m^3/mol; 0 where none is curated.
    """
    # Imported here: the virial table is itself read from the curated datasets.
    from optmat.gas.virial import virial_ratio  # noqa: PLC0415

    return virial_ratio(dataset.material)


def polarizability(dataset: Dataset) -> Float[np.ndarray, " n"]:
    """Return the per-particle polarizability each value implies.

    Parameters
    ----------
    dataset : Dataset
        A dataset of one of :data:`~optmat.vocabulary.POLARIZABILITY_QUANTITIES`.

    Returns
    -------
    numpy.ndarray
        Dipole polarizability, C m^2/V.

    Raises
    ------
    DataError
        If the quantity does not determine the polarizability by itself.
    """
    values = dataset.scaled_values()
    match dataset.quantity:
        case "polarizability":
            return values * POLARIZABILITY_UNITS[dataset.unit]
        case "refractivity":
            return polarizability_from_refractivity(
                values, density_at_state(dataset), _virial(dataset)
            )
        case "molar_refractivity":
            # A_R = N_A alpha / (3 eps0), in cm^3/mol.
            return 3 * EPS0 * values * 1e-6 / N_A
        case "refractivity_per_density":
            # (n - 1) / N in cm^3: n - 1 at the state's density, then Lorentz-Lorenz.
            density = density_at_state(dataset)
            return polarizability_from_refractivity(
                values * 1e-6 * density, density, _virial(dataset)
            )
        case _:
            msg = (
                f"{dataset.key}: {dataset.quantity!r} has no polarizability of its own"
            )
            raise DataError(msg)


def polarizability_uncertainty(dataset: Dataset) -> Float[np.ndarray, " n"]:
    """Return the one-sigma uncertainty of :func:`polarizability`.

    Every conversion is linear in the value to within ``n - 1`` (1e-3), so the relative
    uncertainty carries over.

    Parameters
    ----------
    dataset : Dataset
        A dataset of one of :data:`~optmat.vocabulary.POLARIZABILITY_QUANTITIES`.

    Returns
    -------
    numpy.ndarray
        One-sigma uncertainty, C m^2/V; NaN where none is stated.
    """
    values = dataset.scaled_values()
    return polarizability(dataset) * dataset.standard_uncertainty() / np.abs(values)


def refractivity_at(
    dataset: Dataset, state: State = ONE_ATM_0C
) -> Float[np.ndarray, " n"]:
    """Return ``n - 1`` of the real gas at ``state`` implied by each value.

    Parameters
    ----------
    dataset : Dataset
        A dataset of one of :data:`~optmat.vocabulary.POLARIZABILITY_QUANTITIES`.
    state : State, optional
        The state to evaluate at; 0 degrees C and 1 atm by default.

    Returns
    -------
    numpy.ndarray
        ``n - 1`` from the polarizability by Lorentz-Lorenz at the real-gas density,
        with the gas's refractivity virial term.
    """
    return refractivity(
        polarizability(dataset),
        number_density(dataset.material, state),
        _virial(dataset),
    )

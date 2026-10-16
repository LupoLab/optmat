# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The gas species the package curates data for.

The registry contains five noble gases and ten molecular gases.
Each carries the few fixed facts the rest of the package relies on: its electron count
(the Thomas-Reiche-Kuhn sum rule ``S(0)`` of its dipole oscillator-strength
distribution) and the name of its equation of state in CoolProp, where one exists.
"""

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

from optmat.exceptions import OptmatError

__all__ = ["Gas", "GasInfo", "as_gas", "info"]


class Gas(StrEnum):
    """A gas species, named by its chemical formula."""

    He = "He"
    Ne = "Ne"
    Ar = "Ar"
    Kr = "Kr"
    Xe = "Xe"
    H2 = "H2"
    D2 = "D2"
    N2 = "N2"
    O2 = "O2"
    CH4 = "CH4"
    SF6 = "SF6"
    N2O = "N2O"
    CO2 = "CO2"
    O3 = "O3"
    CO = "CO"


@dataclass(frozen=True)
class GasInfo:
    """Fixed facts about one species.

    Attributes
    ----------
    gas : Gas
        The species.
    electrons : int
        Number of electrons per atom or molecule: the Thomas-Reiche-Kuhn sum rule
        ``S(0)`` of the dipole oscillator-strength distribution.
    fluid : str or None
        CoolProp's name of the fluid for the real-gas equation of state, or None where
        CoolProp has none (ozone).
    atomic : bool
        True for the noble gases, which have no vibrational or rotational structure
        and so no infrared resonances.
    """

    gas: Gas
    electrons: int
    fluid: str | None
    atomic: bool


_INFO: Final[dict[Gas, GasInfo]] = {
    info.gas: info
    for info in (
        GasInfo(Gas.He, 2, "Helium", atomic=True),
        GasInfo(Gas.Ne, 10, "Neon", atomic=True),
        GasInfo(Gas.Ar, 18, "Argon", atomic=True),
        GasInfo(Gas.Kr, 36, "Krypton", atomic=True),
        GasInfo(Gas.Xe, 54, "Xenon", atomic=True),
        # Normal hydrogen and deuterium: CoolProp's equilibrium forms, as in solprop.
        GasInfo(Gas.H2, 2, "Hydrogen", atomic=False),
        GasInfo(Gas.D2, 2, "Deuterium", atomic=False),
        GasInfo(Gas.N2, 14, "Nitrogen", atomic=False),
        GasInfo(Gas.O2, 16, "Oxygen", atomic=False),
        GasInfo(Gas.CH4, 10, "Methane", atomic=False),
        GasInfo(Gas.SF6, 70, "SulfurHexafluoride", atomic=False),
        GasInfo(Gas.N2O, 22, "NitrousOxide", atomic=False),
        GasInfo(Gas.CO2, 22, "CarbonDioxide", atomic=False),
        GasInfo(Gas.O3, 24, None, atomic=False),
        GasInfo(Gas.CO, 14, "CarbonMonoxide", atomic=False),
    )
}


def as_gas(gas: Gas | str) -> Gas:
    """Return the :class:`Gas` named by ``gas``.

    Parameters
    ----------
    gas : Gas or str
        A species or its formula, e.g. ``"Ar"``.

    Returns
    -------
    Gas
        The species.

    Raises
    ------
    OptmatError
        If ``gas`` names no curated species.

    Examples
    --------
    >>> as_gas("Ar")
    <Gas.Ar: 'Ar'>
    """
    try:
        return Gas(gas)
    except ValueError as exc:
        known = ", ".join(g.value for g in Gas)
        msg = f"unknown gas {gas!r}; curated species are {known}"
        raise OptmatError(msg) from exc


def info(gas: Gas | str) -> GasInfo:
    """Return the fixed facts about a species.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    GasInfo
        Its electron count, equation-of-state name and whether it is atomic.

    Examples
    --------
    >>> info("Xe").electrons
    54
    """
    return _INFO[as_gas(gas)]

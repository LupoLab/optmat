# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Huot & Bose (1991): dielectric virial coefficients of the atomic gases.

Table I (``tbl_9b08ee215279``, cells file only; image checked): per gas and
temperature, B_eps (cm^6/mol^2) and C_eps (cm^9/mol^3), static. Each goes to
``<coefficient>_<T>_<gas>`` (``B_eps`` order 2, ``C_eps`` order 3), the static limit as
the abscissa. Cells printed "..." (no value) are skipped.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Huot_1991_ExperimentalDeterminationDielectric"
TABLE_I = "tbl_9b08ee215279"
GASES = {"Helium": "He", "Neon": "Ne", "Argon": "Ar", "Krypton": "Kr", "Xenon": "Xe"}


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``- 0.07 ± 0.01``."""
    value, _, sigma = cell.partition("\N{PLUS-MINUS SIGN}")
    return number(value.replace(" ", "")), sigma.strip()


def extract() -> dict[str, list[Point]]:
    """Return Table I's B_eps and C_eps."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        gas = GASES[cells[0]] if cells[0] else gas
        temperature, printed = cells[1], " | ".join(cells)
        for column, coefficient in ((2, "B_eps"), (3, "C_eps")):
            if cells[column] and cells[column] != "...":
                value, sigma = split(cells[column])
                name = f"{coefficient}_{temperature}_{gas}"
                files[f"{gas}/huot_1991__{name}.csv"] = [
                    Point(temperature, "inf", value, sigma, printed)
                ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

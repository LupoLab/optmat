# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Everall, Partanen, Barr and Shaw (1987): the parameters of their SRS table.

The one table (``tbl_4397977c2c69``, p. 396; no CSV, the JSON cells match the image;
caption 'Comparison of calculated (eq. (7)) and measured thresholds at 5 atm') gives
per gas the cross section dsigma/dOmega (10^-29 cm^2 sr^-1, the column factor the
dataset scale), the Raman linewidth Delta nu_R (cm^-1, FWHM), the dephasing time T2 (ps)
and the calculated and measured threshold energies (uJ). The cross sections, linewidths
and dephasing times of the curated gases CH4, H2, SF6, O2 and N2 are written; NH3 is
not a curated species, and the thresholds have no quantity (the record's notes).

Two linewidths are printed with an approximate sign, ``≈4.0`` (NH3) and ``≈1.5`` (N2):
the number is written and the sign kept in the printed column.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Everall_1987_ThresholdMeasurementsStimulated"
TABLE = "tbl_4397977c2c69"
PUMP = "248"  # nm, the KrF pump of every value
PRESSURE = "5"  # atm, the pressure of the table (its caption)
GASES = {"CH_4": "CH4", "H_2": "H2", "SF_6": "SF6", "O_2": "O2", "N_2": "N2"}
APPROXIMATE = "≈"


def extract() -> dict[str, list[Point]]:
    """Return the cross sections, linewidths and T2 of the curated gases."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        if cells[0] not in GASES:
            continue
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        stem = f"{gas}/everall_1987__"
        files[f"{stem}dsigma_{gas}.csv"] = [
            Point(PUMP, PUMP, number(cells[1]), "", printed)
        ]
        width = number(cells[2].removeprefix(APPROXIMATE))
        files[f"{stem}linewidth_{gas}.csv"] = [
            Point(f"{PRESSURE} atm", PRESSURE, width, "", printed)
        ]
        files[f"{stem}T2_{gas}.csv"] = [
            Point(f"{PRESSURE} atm", PRESSURE, number(cells[3]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dunmur & Jessup (1979), Part I: calculated Kerr virial coefficients.

Both tables have no CSV rendering and are read from their JSON cells; both checked
against their images (no differences; decimal points printed as raised dots in Table 2).
Only Ar, Kr, Xe, CH4 and SF6 are taken (CF4 and neopentane are out of scope).

* Table 2 (``tbl_02b495d28eaa``): the dipole-induced-dipole B_K (10^33 C^2 m^8 J^-2) at
  250-330 K per gas (``table2_B_K_<gas>``), and the last row, "expt [19]", the measured
  B_K of Part II with its temperature in parentheses ("3·8(305)"; the image prints the
  temperature on a second line), one point per gas (``table2_expt_<gas>``).
* Table 4 (``tbl_63fdb7ad5cff``): C_K at 300 K (10^38 C^2 m^11 J^-2): the total
  (``table4_C_K_<gas>``, the contributions in the record's notes) and Part II's
  experimental value where printed (``table4_expt_<gas>``; "—" elsewhere).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Dunmur_1979_InfluenceIntermolecularInteractions"
TABLE_2 = "tbl_02b495d28eaa"
TABLE_4 = "tbl_63fdb7ad5cff"
DOT = "\N{MIDDLE DOT}"
GASES = {"Ar": "Ar", "Kr": "Kr", "Xe": "Xe", "CH_4": "CH4", "SF_6": "SF6"}
_EXPT = re.compile(r"^([\d.]+)\((\d+)\)$")


def _table2() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_2)
    header = grid[1]
    files: dict[str, list[Point]] = {}
    for column, label in enumerate(header):
        gas = GASES.get(label)
        if gas is None:
            continue
        for cells in grid[2:]:
            cell = cells[column].replace(DOT, ".")
            printed = " | ".join(cells)
            if cells[0].endswith("expt [19]"):
                match = _EXPT.match(cell)
                if match is None:
                    msg = f"unexpected expt cell {cells[column]!r}"
                    raise ValueError(msg)
                temperature = number(match[2])
                files[f"{gas}/dunmur_1979__table2_expt_{gas}.csv"] = [
                    Point(temperature, temperature, number(match[1]), "", printed)
                ]
                continue
            temperature = number(cells[0])
            files.setdefault(f"{gas}/dunmur_1979__table2_B_K_{gas}.csv", []).append(
                Point(temperature, temperature, number(cell), "", printed)
            )
    return files


def _table4() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_4)
    if not ("TOTAL" in grid[1][5] and grid[1][6].startswith("Expt.")):
        msg = f"unexpected Table 4 header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        files[f"{gas}/dunmur_1979__table4_C_K_{gas}.csv"] = [
            Point("300", "300", number(cells[5]), "", printed)
        ]
        if cells[6] != "\N{EM DASH}":
            files[f"{gas}/dunmur_1979__table4_expt_{gas}.csv"] = [
                Point("300", "300", number(cells[6]), "", printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 2 and Table 4 files."""
    return _table2() | _table4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

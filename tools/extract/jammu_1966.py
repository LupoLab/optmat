# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Jammu, St. John & Welsh (1966): rotational Raman broadening coefficients.

Can. J. Phys. 44, 797. Both tables exist as cells files only; checked against their
images. The coefficient is Delta nu_1/2 / p, the HALF width at half maximum per atm
(the paper's footnote, p. 2), in 10^-3 cm-1/atm.

* Table I (``tbl_f9fde98532d6``): O2 self-broadening for the S-branch lines K = 3-19,
  from the 20- and 40-minute exposures.
* Table II (``tbl_9e05c97e2664``): the comparison. "This investigation" prints the
  extreme values observed over the rotational lines, "a -> b" (note *: "the values ...
  refer to the extreme values observed"); they are split into the largest (lowest J)
  and the smallest (highest J). Mikhailov's N2 cell holds two values ("47.5 and 26.5":
  up to 25 atm and above, p. 17); Pinter's cells list the J = 4, 8, 12 (N2) and J =
  22, 26, 32 (CO2) lines in parentheses, with a note mark. "—" is no value.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Jammu_1966_PRESSUREBROADENINGROTATIONAL"
TABLE_I = "tbl_f9fde98532d6"
TABLE_II = "tbl_9e05c97e2664"
DASH = "\N{EM DASH}"
ARROW = "\N{RIGHTWARDS ARROW}"
PAIRS = {
    "O_{2}-O_{2}": ("O2", "self"),
    "N_{2}-N_{2}": ("N2", "self"),
    "CO_{2}-CO_{2}": ("CO2", "self"),
    "CO-CO": ("CO", "self"),
    "O_{2}-Ar": ("O2", "Ar"),
    "N_{2}-Ar": ("N2", "Ar"),
    "O_{2}-He": ("O2", "He"),
    "N_{2}-He": ("N2", "He"),
}
# The fourth column's author, and the lines his values are for (Table II notes double
# dagger and section sign).
FOURTH = "Pinter"  # codespell:ignore pinter
PINTER_LINES = {"N2": ("4", "8", "12"), "CO2": ("22", "26", "32")}


def _name(material: str, perturber: str) -> str:
    return material if perturber == "self" else f"{material}_{perturber}"


def _table_i() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[1][1:] != ["20-min exposure", "40-min exposure"]:
        msg = f"unexpected Table I header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        row = f"S({number(cells[0])})"
        for column, name in ((1, "20min"), (2, "40min")):
            files.setdefault(f"O2/jammu_1966__table1_{name}.csv", []).append(
                Point(row, row, number(cells[column]), "", " | ".join(cells))
            )
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    heads = [cell.split(" (")[0] for cell in grid[1][1:]]
    expected = ["This investigation*", "Mikhailov", "Lazarev", FOURTH]
    if heads != expected:
        msg = f"unexpected Table II header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        material, perturber = PAIRS[cells[0]]
        name = _name(material, perturber)
        printed = " | ".join(cells)
        high, low = (number(part) for part in cells[1].split(ARROW))
        files[f"{material}/jammu_1966__table2_{name}.csv"] = [
            Point("largest", "largest", high, "", printed),
            Point("smallest", "smallest", low, "", printed),
        ]
        mikhailov = cells[2]
        if mikhailov not in {DASH, ""}:
            values = mikhailov.split(" and ")
            rows = ["single"] if len(values) == 1 else ["up to 25 atm", "above 25 atm"]
            files[f"{material}/jammu_1966__table2_mikhailov_{name}.csv"] = [
                Point(row, row, number(value), "", printed)
                for row, value in zip(rows, values, strict=True)
            ]
        if cells[3] not in {DASH, ""}:
            files[f"{material}/jammu_1966__table2_lazarev_{name}.csv"] = [
                Point(
                    "most intense lines",
                    "most intense lines",
                    number(cells[3]),
                    "",
                    printed,
                )
            ]
        if cells[4] not in {DASH, ""}:
            values = re.sub(r"[()]|\\ddagger|§", "", cells[4]).split(", ")
            files[f"{material}/jammu_1966__table2_pinter_{name}.csv"] = [
                Point(f"S({j})", f"S({j})", number(value), "", printed)
                for j, value in zip(PINTER_LINES[material], values, strict=True)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I and II files."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

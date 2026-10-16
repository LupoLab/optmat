# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Reichert et al. (2015): Table 1, Delta alpha, B and D of N2, O2 (p. 22232).

Table 1 lists, per molecule, the polarizability anisotropy used as the self-reference of
the beam-deflection measurement (1e-40 F m2, quoted from Bridge and Buckingham, the
reference mark "[35]" printed in the cell), the rotational constant B and centrifugal
constant D fitted here (columns marked "a") and the Raman values the text compares
them with (columns "b", Weber [41]). The digits in parentheses are the uncertainty in
the last digits ("1.9896(2)" is 1.9896 +- 0.0002); they become the sigma column. The
CSV matches the table image (checked).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Reichert_2015_BeamDeflectionMeasurement"
TABLE = "tbl_f40cb41739d3"
GASES = {"14N2": "N2", "16O2": "O2"}
HEADER = [
    "Molecule",
    "Δα (10-40 Fm2)",
    "B (cm-1)a",
    "B (cm-1)b",
    "D (10-6 cm-1)a",
    "D (10-6 cm-1)b",
]
# (column, dataset prefix) of the B and D columns.
CONSTANTS = ((2, "B_fit"), (3, "B_raman"), (4, "D_fit"), (5, "D_raman"))
_PARENTHESES = re.compile(r"^([\d.]+)\((\d+)\)$")


def _with_sigma(cell: str) -> tuple[str, str]:
    """Return the value and uncertainty of a "1.9896(2)" cell, in the printed unit."""
    match = _PARENTHESES.match(cell.strip())
    if match is None:
        msg = f"not a 'value(n)' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1))
    decimals = len(value.partition(".")[2])
    sigma = f"{int(match.group(2)) * 10**-decimals:.{decimals}f}"
    return value, sigma


def extract() -> dict[str, list[Point]]:
    """Return the reproduced anisotropy and the B and D of each gas."""
    grid = table(DIRECTORY, TABLE)
    if grid[0] != HEADER:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files = {}
    for cells in grid[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        value = number(cells[1].removesuffix("[35]"))
        path = f"{gas}/reichert_2015__dalpha_bridge_buckingham_{gas}.csv"
        files[path] = [Point("inf", "not stated", value, "", printed)]
        for column, name in CONSTANTS:
            value, sigma = _with_sigma(cells[column])
            label = name.partition("_")[0]
            point = Point(label, label, value, sigma, printed)
            files[f"{gas}/reichert_2015__{name}_{gas}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

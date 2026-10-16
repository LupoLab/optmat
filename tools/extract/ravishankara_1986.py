# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ravishankara et al. (1986): N2O5 photolysis quantum yields of O(3P) and NO3.

Both tables are cells only (``.json``; images checked, pp. 5356 and 5359).

* **Table 1** (``tbl_92212005c108``): per photolysis wavelength, the individual O(3P)
  yields, their average +- 2 sigma, and in the last column "Phi +- Accuracy". The
  wavelength is printed once, on the first row of its group (the 266 nm group spans
  two sets of runs and its final value sits on a row with no wavelength). The last
  column's four values are written.
* **Table 2** (``tbl_eaa82cc4a2ec``): the NO3 yield at 248 nm; the last cell of the
  last column holds three results, of which "(0.96 +- 0.13)^d", the average with random
  and systematic errors combined, is written.

Register: TC-975.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Ravishankara_1986_N2O5PhotolysisQuantum"
TABLE_1 = "tbl_92212005c108"
TABLE_2 = "tbl_eaa82cc4a2ec"
_ACCURACY = re.compile(r"^([\d.]+) \\pm ([\d.]+)$")
_NO3 = re.compile(r"\(([\d.]+) \\pm ([\d.]+)\)\^\{d\}$")


def _o3p() -> list[Point]:
    """Return Table 1's O(3P) yields, "Phi +- Accuracy", by wavelength."""
    points = []
    wavelength = ""
    for cells in json_table(DIRECTORY, TABLE_1)[1:]:
        wavelength = cells[0] or wavelength
        if not cells[-1]:
            continue
        match = _ACCURACY.match(cells[-1])
        if match is None:
            msg = f"unreadable cell {cells[-1]!r}"
            raise ValueError(msg)
        printed = f"{wavelength} nm | " + " | ".join(cells)
        points.append(
            Point(wavelength, wavelength, match.group(1), match.group(2), printed)
        )
    return points


def _no3() -> list[Point]:
    """Return Table 2's NO3 yield at 248 nm, with random and systematic errors."""
    cell = json_table(DIRECTORY, TABLE_2)[-1][-1]
    match = _NO3.search(cell)
    if match is None:
        msg = f"no footnote-d average in {cell!r}"
        raise ValueError(msg)
    return [Point("248", "248", match.group(1), match.group(2), cell)]


def extract() -> dict[str, list[Point]]:
    """Return the O(3P) and NO3 quantum yields."""
    return {
        "N2O5/ravishankara_1986__o3p_yield.csv": _o3p(),
        "N2O5/ravishankara_1986__no3_yield.csv": _no3(),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

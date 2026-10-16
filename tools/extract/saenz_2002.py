# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Saenz (2002): static-field ionisation rates of H2 by vibrational state.

Table I (``tbl_b2f3d07790b3``, JSON cells only: the field cell spans four rows) gives,
for F = 0.05338-0.18 au and v = 0, 1, 2, 4, 9, 14, four rows per field: (a) the
frozen-nuclei ab initio rate, (b) the frozen-nuclei ADK rate with the R-dependent
ionisation potential, (c) the fully adiabatic ab initio rate and (d) the ADK rate with
the R-independent ionisation potential of the state. One file per scheme and v
(``tableI_<scheme>_v<v>``); a blank cell (no bound state, or not computed) is skipped.
Values are printed "m[-k]" (m x 10^-k) and written as "me-k"; the two (d) cells at F =
0.18 au for v = 9 and 14 are printed without a power ("1.203", "1.553"). Checked against
the table image (2026-09-26): the cells agree.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Saenz_2002_BehaviorMolecularHydrogen"
TABLE = "tbl_b2f3d07790b3"
SCHEMES = {"(a)": "a", "(b)": "b", "(c)": "c", "(d)": "d"}
STATES = ("0", "1", "2", "4", "9", "14")
_POWER = re.compile(r"^(\d+\.\d+)\[-(\d+)\]$")


def _value(cell: str) -> str:
    """Return a printed rate "m[-k]" as "me-k"; a plain number as printed."""
    match = _POWER.match(cell)
    if match:
        return f"{match.group(1)}e-{match.group(2)}"
    float(cell)  # a plain number, or an error
    return cell


def extract() -> dict[str, list[Point]]:
    """Return one file per scheme (a)-(d) and vibrational state."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[0] != ["F (a.u.)", "", *(f"v={v}" for v in STATES)]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"H2/saenz_2002__tableI_{scheme}_v{v}.csv": []
        for scheme in SCHEMES.values()
        for v in STATES
    }
    field = ""
    for cells in grid[1:]:
        field = cells[0] or field  # the field cell spans the four rows of an entry
        scheme = SCHEMES[cells[1]]
        printed = f"{field} | " + " | ".join(cells[1:])
        for v, cell in zip(STATES, cells[2:], strict=True):
            if cell:
                point = Point(field, field, _value(cell), "", printed)
                files[f"H2/saenz_2002__tableI_{scheme}_v{v}.csv"].append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Li (1980): recommended n and dn/dT of silicon at 293 K.

Only the 293 K columns of Li's Tables 1 and 2 are extracted (the recommended values
for silicon; the record's notes say why the other columns and the appendix tables of
others' measurements are not). Both tables are JSON only in the library (merged
header cells): ``tbl_67c9bf1134a1`` (Table 1, n, PDF p. 17) and ``tbl_c785eb5d7e50``
(Table 2, dn/dT in 10^-4 K^-1, PDF p. 19); rows are the wavelengths in um, 1.20-14.00
(35 rows), and the second header row the 14 temperatures, 100-750 K.

The pages are scans (one 1-bit image each, 426 ppi), typewriter digits with an
overstrike on the digits Li calls insignificant, and a scratch across rows 1.30-1.60
of both tables. Every 293 K cell was compared with the page image (the scan's own
pixels, enlarged) and with Li's Eq. (22) and its T-derivative (see the record): the
library misreads one cell of each column, corrected here in ``CORRECTIONS`` as
(library cell, cell on the page):

- Table 1, 1.50 um: "3.4796" in the library; the page prints 3.4799 (under the
  scratch, with the last digit overstruck), which Eq. (22) gives (3.47986) and
  which the neighbouring columns' steps confirm.
- Table 2, 1.20 um: "1.963" in the library; the page and the PDF's text layer print
  1.983, which the T-derivative of Eq. (22) gives (1.984).

Table 2's 1.36 um cell (1.892) lies under the scratch; its upper half reads 1.892 and
the derivative of Eq. (22) gives 1.893.
"""

from typing import Final

from _json_table import json_table
from _library import Point, number, write

DIRECTORY: Final = "Li_1980_RefractiveIndexSilicon"
TABLE_1: Final = "tbl_67c9bf1134a1"
TABLE_2: Final = "tbl_c785eb5d7e50"
TEMPERATURES: Final = (
    "100",
    "150",
    "200",
    "250",
    "293",
    "350",
    "400",
    "450",
    "500",
    "550",
    "600",
    "650",
    "700",
    "750",
)
COLUMN: Final = 1 + TEMPERATURES.index("293")

# (library table, row label) -> (library cell, cell on the page), 293 K column.
CORRECTIONS: Final = {
    (TABLE_1, "1.50"): ("3.4796", "3.4799"),
    (TABLE_2, "1.20"): ("1.963", "1.983"),
}


def _column(table_id: str, name: str) -> dict[str, list[Point]]:
    """Return the 293 K column of Table 1 or 2 as one file."""
    rows = json_table(DIRECTORY, table_id)
    if tuple(rows[1][1:15]) != TEMPERATURES or len(rows) != 37:
        msg = f"unexpected layout of {table_id}: {rows[1]!r}, {len(rows)} rows"
        raise ValueError(msg)
    points = []
    for cells in rows[2:]:
        wavelength = number(cells[0])
        cell = cells[COLUMN]
        correction = CORRECTIONS.get((table_id, cells[0]))
        if correction is not None:
            if cell != correction[0]:
                msg = f"{table_id} {cells[0]}: cell {cell!r}, not {correction[0]!r}"
                raise ValueError(msg)
            cell = correction[1]
        printed = " | ".join(c for c in cells[:15])
        points.append(Point(wavelength, wavelength, number(cell), "", printed))
    return {f"Si/li_1980__{name}.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return the 293 K columns of Tables 1 (n) and 2 (dn/dT)."""
    return {**_column(TABLE_1, "n_293K"), **_column(TABLE_2, "dndt_293K")}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

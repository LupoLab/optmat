# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Au, Cooper, Burton, Olney & Brion (1993): CH4 absolute oscillator strengths.

Only methane is a curated species; the other alkanes (C2H6 ... C8H18, Tables 4-6 and
8-14, their rows of Tables 1-2) are described in the record's notes. Images checked
against the cells (2026-09-27).

- Table 7 (``tbl_f5a4c303b9e9``, CSV, unlabelled in the extraction): the CH4
  low-resolution (1 eV FWHM) df/dE in ``10^-2 eV^-1`` at 7.5-220 eV, four column pairs
  read pair by pair (print order).
- Table 1 (``tbl_8bb88b5f26fe``, cells file): the CH4 rows, the vertical ionisation
  potentials the paper quotes, 1t2 (13.6, 14.4, 15.0 eV, from Kimura et al.'s
  handbook, Ref. [32]; the three rows are printed in one cell group) and 2a1 (22.91 eV,
  Potts and Streets, Ref. [33]): one file per reference, rows ``1t2_1``, ``1t2_2``,
  ``1t2_3`` and ``2a1``. The rows below CH4 belong to it until the next molecule's
  name.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Au_1993_ValenceShellPhotoabsorption"
TABLE_1 = "tbl_8bb88b5f26fe"
TABLE_7 = "tbl_f5a4c303b9e9"
PREFIX = "CH4/au_1993__"


def _table_7() -> list[Point]:
    rows = table(DIRECTORY, TABLE_7)[1:]
    return [
        Point(
            number(cells[2 * k]),
            number(cells[2 * k]),
            number(cells[2 * k + 1]),
            "",
            f"{cells[2 * k].strip()} | {cells[2 * k + 1].strip()}",
        )
        for k in range(4)
        for cells in rows
        if cells[2 * k].strip()
    ]


def _table_1() -> dict[str, list[Point]]:
    outer: list[Point] = []  # Ref. [32]
    inner: list[Point] = []  # Ref. [33]
    in_methane = False
    for cells in json_table(DIRECTORY, TABLE_1)[2:]:
        if cells[0]:
            in_methane = cells[0] == "CH_4"
        if not in_methane:
            continue
        printed = " | ".join(cells)
        if cells[2]:  # outer valence, 1t2 (Jahn-Teller components on three rows)
            row = f"1t2_{len(outer) + 1}"
            outer.append(Point(row, "1t2", number(cells[2]), "", printed))
        if cells[4]:  # inner valence, 2a1
            inner.append(Point("2a1", "2a1", number(cells[4]), "", printed))
    return {
        f"{PREFIX}vertical_outer.csv": outer,
        f"{PREFIX}vertical_inner.csv": inner,
    }


def extract() -> dict[str, list[Point]]:
    """Return CH4's df/dE (Table 7) and the quoted vertical IPs (Table 1)."""
    return {f"{PREFIX}continuum.csv": _table_7()} | _table_1()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

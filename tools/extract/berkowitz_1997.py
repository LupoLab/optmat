# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Berkowitz (1997): helium discrete oscillator strengths and sum rules S(p).

- Table 1 (``tbl_989659e9200e``): oscillator strengths of 1 1S -> n 1P, n = 2-20, from
  five sources (columns b Schiff et al. 1971, c Kono and Hattori 1984, d Drake 1996,
  e Fernley et al. 1987, f Chan et al. 1991, measured, "value +- error"). The
  parenthesised column-c entries (n = 10-20) are the author's extrapolation of Kono and
  Hattori's values and go to their own dataset. The rows are indexed by n (the
  excitation energies E/R, in Rydberg units, have no wavelength unit in the format and
  stay in the printed cells). The table prints digit groups separated by thin spaces
  (``0.276 1647``); the spaces are removed. The sum rows are not extracted.
- Table 2 (``tbl_2a3db8fd2b0e``): the excess oscillator strength f_xs (x 1e4) of the
  2snp autoionising series at its photon energies (eV); the weighted columns and the
  sum row are not extracted.
- Table 3 (``tbl_6523ca2ec08b``, cells file only): the row "Total sum" (the author's
  recommended totals, first cell of each S(p)) and the row "Exact sum" (h Miller, i
  Pekeris); S(p) is in Rydberg units, one file per p. The partial sums, whose cells
  merge two values (checked on the image), are not extracted.

The table images agree with the extracted cells.
"""

import json
import re

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Berkowitz_1997_SumRulesOscillator"
TABLE_1 = "tbl_989659e9200e"
TABLE_2 = "tbl_2a3db8fd2b0e"
TABLE_3 = "tbl_6523ca2ec08b"
PREFIX = "He/berkowitz_1997__"
# Table 1 column -> dataset.
TABLE_1_COLUMNS = {
    3: "f_schiff_1971",
    4: "f_kono_hattori_1984",
    5: "f_drake_1996",
    6: "f_fernley_1987",
    7: "f_chan_1991",
}
# Table 3 column (the first cell of each S(p)) -> order p and file tag.
TABLE_3_ORDERS = {
    1: ("-2", "m2"),
    3: ("-1", "m1"),
    5: ("0", "0"),
    7: ("1", "1"),
    9: ("2", "2"),
}
_MARK = re.compile(r"\^[a-i]$")


def _cells_grid(directory: str, table_id: str) -> list[list[str]]:
    """Return a table's rows from its cells file, for tables without a CSV."""
    path = LIBRARY / directory / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    n_rows = 1 + max(c["row"] for c in cells)
    n_columns = max(c["column"] + c["column_span"] for c in cells)
    grid = [[""] * n_columns for _ in range(n_rows)]
    for cell in cells:
        grid[cell["row"]][cell["column"]] = cell["text"]
    return grid


def _compact(cell: str) -> str:
    """Remove the thin-space digit grouping of a printed number."""
    return cell.strip().replace(" ", "").replace("\N{THIN SPACE}", "")


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_1)[1:]:
        n = cells[0].strip()
        if not n.isdigit():
            continue  # the three sum rows
        printed = " | ".join(cells)
        for column, column_dataset in TABLE_1_COLUMNS.items():
            cell = cells[column].strip()
            if not cell:
                continue
            value, _, sigma = cell.partition("\N{PLUS-MINUS SIGN}")
            value, dataset = _compact(value), column_dataset
            if value.startswith("("):  # the author's extrapolation, n = 10-20
                dataset, value = "f_extrapolated", value.strip("()")
            point = Point(n, n, number(value), _compact(sigma), printed)
            files.setdefault(f"{PREFIX}{dataset}.csv", []).append(point)
    return files


def _table_2() -> list[Point]:
    return [
        Point(c[0].strip(), number(c[1]), number(c[2]), "", " | ".join(c))
        for c in table(DIRECTORY, TABLE_2)[1:]
        if c[0].strip().isdigit()
    ]


def _table_3() -> dict[str, list[Point]]:
    rows = {r[0].strip(): r for r in _cells_grid(DIRECTORY, TABLE_3)}
    total, exact = rows["Total sum"], rows["Exact sum"]
    files = {}
    for column, (order, tag) in TABLE_3_ORDERS.items():
        files[f"{PREFIX}sum_rule_{tag}.csv"] = [
            Point(order, order, number(total[column]), "", " | ".join(total))
        ]
        source = "miller" if exact[column].endswith("^h") else "pekeris"
        value = _MARK.sub("", exact[column].strip())
        files[f"{PREFIX}{source}_sum_rule_{tag}.csv"] = [
            Point(order, order, number(value), "", " | ".join(exact))
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1, 2 (f_xs) and 3 (totals and exact values)."""
    files = _table_1()
    files[f"{PREFIX}f_excess.csv"] = _table_2()
    files |= _table_3()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

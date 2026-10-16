# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dalgarno & Kingston (1960): helium sums S(k) and static polarizabilities.

- Table 1 (``tbl_61824749ae46``, cells file only): S(k) for k = 0, -1, -2, -4, -6,
  "exact" and from the paper's interpolation formula, Eq. (18), for H and He. Only
  helium is a curated species (H is atomic hydrogen). S(k) is in Rydberg units (the
  paper measures energies in rydbergs), so each k goes to its own file and the
  record's scale converts it to atomic units.
- Table 2 (``tbl_07f29f0464dd``, cells file only): static polarizabilities of He-Xe
  in units of a0^3, one file per gas.

The tables print the decimal point as a raised dot, read as a point; the images agree
with the cells.
"""

import json

from _library import LIBRARY, Point, number, write

DIRECTORY = "Dalgarno_1960_RefractiveIndicesVerdet"
TABLE_1 = "tbl_61824749ae46"
TABLE_2 = "tbl_07f29f0464dd"
PREFIX = "{gas}/dalgarno_1960__"
TAGS = {"0": "0", "-1": "m1", "-2": "m2", "-4": "m4", "-6": "m6"}
HELIUM_COLUMNS = {3: "s_exact", 4: "s_eq18"}


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


def _decimal(cell: str) -> str:
    """Read the printed raised decimal point as a decimal point."""
    return number(cell.replace("\N{MIDDLE DOT}", "."))


def extract() -> dict[str, list[Point]]:
    """Return helium's S(k) (Table 1) and the five polarizabilities (Table 2)."""
    files = {}
    prefix = PREFIX.format(gas="He")
    for cells in _cells_grid(DIRECTORY, TABLE_1)[2:]:  # two header rows
        k = number(cells[0])
        for column, stem in HELIUM_COLUMNS.items():
            files[f"{prefix}{stem}_{TAGS[k]}.csv"] = [
                Point(k, k, _decimal(cells[column]), "", " | ".join(cells))
            ]
    gases, values = _cells_grid(DIRECTORY, TABLE_2)[:2]
    for cell, value in zip(gases, values, strict=True):
        gas = cell.strip()
        files[f"{PREFIX.format(gas=gas)}alpha_{gas}.csv"] = [
            Point("1", "inf", _decimal(value), "", " | ".join(values))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Jagau (2016), supplement: all-electron static-field ionisation rates.

Tables S3, S4 and S5 of the supplementary material (``supplement_01``, its pp. 6-7)
list, against the field F (au), the Stark shift Delta E and the ionisation rate Gamma
(au) of the ground states of He (``tbl_835beff31b7f``), Ne (``tbl_bdfc58139762``) and
Ar (``tbl_e8cea2aae6b7``), each at three levels of theory: cs-HF, cs-MP2 and the
correlated level (cs-FCI for He, cs-CCSD for Ne and Ar). One file per gas and level
of the rate, ``tableS<n>_<level>_<gas>``. A cell printed as a dash (no resonance
found at that level and field) is skipped. The Stark shifts are kept in the
``printed`` cells only.

The tables have a two-row header with merged cells, so the library holds them as
cells (``.json``) and not as ``.csv``; the script builds the grid from the cells and
checks the header. Every cell of the three tables was checked against the table
image (2026-09-28): they agree.
"""

import json

from _library import LIBRARY, Point, number, write

DIRECTORY = "Jagau_2016_InvestigatingTunnelAbove/supplement_01"
# Table, gas, and the label of the correlated level as the header prints it.
TABLES = (
    ("tableS3", "tbl_835beff31b7f", "He", "FCI"),
    ("tableS4", "tbl_bdfc58139762", "Ne", "CCSD"),
    ("tableS5", "tbl_e8cea2aae6b7", "Ar", "CCSD"),
)
DASH = "\N{EM DASH}"
HEADER_ROWS = 2
RATE_COLUMNS = (4, 5, 6)
WIDTH = 7


def _grid(table_id: str) -> list[list[str]]:
    """Return a table's cells as rows of strings; a merged cell fills its first."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    rows = 1 + max(cell["row"] for cell in cells)
    grid = [[""] * WIDTH for _ in range(rows)]
    for cell in cells:
        grid[cell["row"]][cell["column"]] = cell["text"].strip()
    return grid


def _table(stem: str, table_id: str, gas: str, correlated: str) -> dict[str, list]:
    grid = _grid(table_id)
    levels = ("HF", "MP2", correlated)
    if grid[0] != ["F", "\N{GREEK CAPITAL LETTER DELTA}E", "", "", "Γ", "", ""] or grid[
        1
    ] != ["", *levels, *levels]:
        msg = f"unexpected header of {stem}: {grid[:HEADER_ROWS]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, level in zip(RATE_COLUMNS, levels, strict=True):
        points = []
        for cells in grid[HEADER_ROWS:]:
            if cells[column] == DASH:
                continue  # no value at this level and field
            field = number(cells[0])
            printed = " | ".join(cells)
            points.append(Point(field, field, number(cells[column]), "", printed))
        files[f"{gas}/jagau_2016__{stem}_{level.lower()}_{gas}.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the rates of Tables S3 to S5, one file per gas and level of theory."""
    files: dict[str, list[Point]] = {}
    for stem, table_id, gas, correlated in TABLES:
        files |= _table(stem, table_id, gas, correlated)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

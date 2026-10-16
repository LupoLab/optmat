# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Teal & MacWood (1935), J. Chem. Phys. 3, 760: Raman lines of H2 and D2.

Both tables exist only as JSON cells; both checked against their images (the cells
agree). HD's rows are not taken (the record's notes).

* Table I (``tbl_cd69d07ceea3``): the pure rotational lines (v 0-0) in vacuum
  wavenumbers, per plate, their AVERAGE (the value taken), the calculated value and
  the difference. Rows are the printed J pair, initial-final ("3-1" an anti-Stokes
  line, "1-3" the Stokes line of the same pair of levels); the molecule cell spans its
  rows.
* Table II (``tbl_4c1e371a6dd6``): the lines of the 0-1 band, the AVERAGE (observed)
  and the value calculated with this paper's constants ("T. & M."), which is a check
  and is written also for the unobserved 0-0 lines. The values calculated with the
  older constants ("U. & T.") and the differences stay in ``printed``.

Table III (vibrational constants) is a formula in the record; the Delta G_1 values are
inline points (text, p. 763).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Teal_1935_RamanSpectraIsotopic"
TABLE_1 = "tbl_cd69d07ceea3"
TABLE_2 = "tbl_4c1e371a6dd6"
PREFIX = "{gas}/teal_1935__"
GASES = {"{\\mathrm{H}}_{2}": "H2", "{\\mathrm{D}}_{2}": "D2", "H_2": "H2", "D_2": "D2"}
DASH = "\N{EM DASH}"


def _rows(grid: list[list[str]]) -> list[tuple[str, list[str]]]:
    """Return (molecule, cells) for each row, the spanned molecule cell filled in."""
    rows = []
    molecule = ""
    for cells in grid:
        if cells[0]:
            molecule = GASES.get(cells[0], cells[0])
        rows.append((molecule, cells))
    return rows


def _table1() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_1)
    header = ["", "v", "J", "PLATE I", "PLATE II", "PLATE III", "AVERAGE", "CALC."]
    if grid[0][:8] != header:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for molecule, cells in _rows(grid[1:]):
        if molecule not in {"H2", "D2"}:
            continue  # HD: notes only
        label = cells[2]
        point = Point(label, label, number(cells[6]), "", " | ".join(cells))
        path = PREFIX.format(gas=molecule) + f"table1_{molecule}.csv"
        files.setdefault(path, []).append(point)
    return files


def _table2() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_2)
    if grid[1][7:9] != ["U. & T.", "T. & M."] or grid[0][6] != "AVERAGE":
        msg = f"unexpected Table II header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for molecule, cells in _rows(grid[2:]):
        if molecule not in {"H2", "D2"}:
            continue  # HD: notes only
        label = cells[2]
        printed = " | ".join(cells)
        prefix = PREFIX.format(gas=molecule)
        if cells[6] != DASH:  # the 0-0 lines of H2 and D2 were not observed
            point = Point(label, label, number(cells[6]), "", printed)
            files.setdefault(f"{prefix}table2_{molecule}.csv", []).append(point)
        point = Point(label, label, number(cells[8]), "", printed)
        files.setdefault(f"{prefix}table2_calc_{molecule}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Tables I and II."""
    return _table1() | _table2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

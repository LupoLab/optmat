# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Stoicheff (1957), Can. J. Phys. 35, 730: Raman lines and constants, H2, D2.

All three tables exist only as JSON cells; each was checked against its image. HD's
columns are not taken (the record's notes).

* Table I (``tbl_8dcc7741bb44``, p. 732): the observed wavenumber displacements (cm-1)
  of the pure rotational S_0(J) lines and the Q_1(J) lines of the 1-0 band. **Structural
  fix:** the library merged H2's |Delta nu| column into the designation cell (a LaTeX
  ``cases`` block holding "S_0(0) & 354.381 \\ ..."), and every later column moved one
  place left: cells 2-3 are H2's O-C_1 and O-C_2, 4-5 HD's |Delta nu| and O-C_1, 6-7
  D2's, and cell 8 is empty. The image shows H2's four S_0 and four Q_1 lines, HD's
  five S_0 and four Q_1, and D2's five of each. H2's numbers are read from the ``cases``
  block (and the count checked); the O-C columns stay in ``printed``.
* Table II (``tbl_6bd75ebe07be``, p. 734): B_v, D_v, H_v and nu_{v-0} of H2 (v = 0-3)
  and D2 (v = 0, 1). The last digit of most constants is printed as a subscript (the
  image: 59.339 with a subscript 2); the cells join it to the number, as the record
  does. H2's v = 1 cells print the value from the 1-0 band alone and, in parentheses,
  the value refitted with D_1 from the smooth D_v curve (p. 735): the two go to
  separate files. **Structural fixes** (``MISREAD``): H2's D_2 cell reads "0.0427;"
  where the image prints 0.0427 with a subscript 7 (0.04277), and D2's D_1 cell reads
  "0.01073" where the image prints 0.0107 with a subscript 5 (0.01075, as the text's
  D_0 - D_1 = 0.00059 cm-1 on p. 737 requires). The
  effective I_0 and r_0 rows have no quantity (the record's notes).
* Table III (``tbl_f87804e5f4e4``, p. 736): the experimental and the corrected
  equilibrium rotational constants B_e and B_e* (the dagger on D2's B_e* is note text;
  I_e and r_e rows are notes).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Stoicheff_1957_HIGHRESOLUTIONRAMAN"
TABLE_I = "tbl_8dcc7741bb44"
TABLE_II = "tbl_6bd75ebe07be"
TABLE_III = "tbl_f87804e5f4e4"
DAGGER = "\N{DAGGER}"
# Table II cells whose subscript last digit the library misread (cell, image).
MISREAD = {
    # H2 D_2: the image prints 0.0427 with a subscript 7, read as a semicolon.
    ("H2", "D", "2"): ("0.0427;", "0.04277"),
    # D2 D_1: the image prints 0.0107 with a subscript 5 (the glyph of the 5 in
    # 29.910_5 beside it), read as a 3; the text (p. 737) confirms it: "D_0 - D_1 =
    # 0.00059 cm-1" = 0.01134 - 0.01075.
    ("D2", "D", "1"): ("0.01073", "0.01075"),
}
_CASES_ENTRY = re.compile(r"([SQ])_([01])\((\d)\)\s*&\s*([\d.]+)")


def _printed(cells: list[str]) -> str:
    return " | ".join(cells)


def _table_i() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0][0] != "Designation" or grid[1][2:] != [
        "|\\Delta\\nu|",
        "O-C_1",
        "O-C_2",
        "|\\Delta\\nu|",
        "O-C_1",
        "|\\Delta\\nu|",
        "O-C_1",
    ]:
        msg = f"unexpected Table I header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    # Two blocks of five rows: the pure rotational S_0 lines, then the Q_1 lines.
    for first, branch, band in ((2, "S", "S0"), (7, "Q", "Q1")):
        block = grid[first : first + 5]
        # H2's numbers sit in the cases cell of the block's first row (see above).
        h2 = [
            (f"{kind}_{v}({j})", value)
            for kind, v, j, value in _CASES_ENTRY.findall(block[0][1])
        ]
        if len(h2) != 4 or any(not label.startswith(branch) for label, _ in h2):
            msg = f"unexpected H2 cases cell {block[0][1]!r}"
            raise ValueError(msg)
        for j, cells in enumerate(block):
            label = f"{band}({j})"
            if j < len(h2) and h2[j][0] != f"{branch}_{band[1]}({j})":
                msg = f"H2 row {j}: {h2[j][0]!r} is not {label!r}"
                raise ValueError(msg)
            h2_cell = h2[j][1] if j < len(h2) else ""
            # printed: the row as the image lays it out, H2's number restored.
            printed = _printed([label, h2_cell, *cells[2:8]])
            if h2_cell:
                files.setdefault(f"H2/stoicheff_1957__table1_{band}_H2.csv", []).append(
                    Point(label, label, number(h2_cell), "", printed)
                )
            if cells[6]:  # D2's |Delta nu| (cell 6 after the shift)
                files.setdefault(f"D2/stoicheff_1957__table1_{band}_D2.csv", []).append(
                    Point(label, label, number(cells[6]), "", printed)
                )
    return files


def _split(cell: str) -> tuple[str, str]:
    """Return (value, bracketed value) of a "v(w)" cell, or (v, "")."""
    value, _, bracketed = cell.partition("(")
    return value, bracketed.removesuffix(")")


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    columns = {"H2": [1, 2, 3, 4], "D2": [7, 8]}
    if (
        grid[0][1] != "H_2"
        or grid[0][7] != "D_2"
        or grid[1][1:] != [f"v = {v}" for v in (0, 1, 2, 3, 0, 1, 0, 1)]
    ):
        msg = f"unexpected Table II header {grid[:2]!r}"
        raise ValueError(msg)
    rows = {cells[0]: cells for cells in grid[2:]}
    names = {"B_v": "B", "D_v": "D", "H_v": "H", "\\nu_{v-0}": "nu"}
    files: dict[str, list[Point]] = {}
    for molecule, indices in columns.items():
        for key, name in names.items():
            cells = rows[key]
            for index in indices:
                cell = cells[index]
                if not cell:
                    continue
                v = grid[1][index].removeprefix("v = ")
                if (molecule, name, v) in MISREAD:
                    misread, image = MISREAD[molecule, name, v]
                    if cell != misread:
                        msg = f"{molecule} {name}_{v} cell changed: {cell!r}"
                        raise ValueError(msg)
                    cell = image
                value, bracketed = _split(cell)
                printed = f"{key} | v = {v} | {cells[index]}"
                if name == "nu":
                    stem = f"table2_nu_{molecule}"
                elif v == "0":
                    stem = f"table2_{name}0_{molecule}"
                else:
                    stem = f"table2_{name}v_{molecule}"
                files.setdefault(f"{molecule}/stoicheff_1957__{stem}.csv", []).append(
                    Point(v, v, number(value), "", printed)
                )
                if bracketed:
                    files.setdefault(
                        f"{molecule}/stoicheff_1957__table2_{name}1_refit_{molecule}.csv",
                        [],
                    ).append(Point(v, v, number(bracketed), "", printed))
    return files


def _table_iii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    if grid[0] != ["", "H_{2}", "HD", "D_{2}"]:
        msg = f"unexpected Table III header {grid[0]!r}"
        raise ValueError(msg)
    rows = {cells[0]: cells for cells in grid[1:]}
    files: dict[str, list[Point]] = {}
    for key, stem in (
        ("B_{e} (cm^{-1})", "Be"),
        ("B_{e}^{*} (cm^{-1})", "Be_corrected"),
    ):
        cells = rows[key]
        for molecule, index in (("H2", 1), ("D2", 3)):
            value = cells[index].removesuffix(DAGGER)  # the dagger is D2's note
            files[f"{molecule}/stoicheff_1957__table3_{stem}_{molecule}.csv"] = [
                Point(molecule, molecule, number(value), "", _printed(cells))
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Tables I-III."""
    return _table_i() | _table_ii() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

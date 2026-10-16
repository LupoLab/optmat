# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rice, Taylor, Lee & Almlof (1991): static gamma of the noble gases.

Tables I-V (Ar, Kr, Xe, Ne, He: alpha and gamma by basis set and method) and Table VI
(best estimates and experimental comparisons). Every image checked against the cells:
no differences. Tables I-IV are cells files only (merged "16/26/10 electrons
correlated" and "Relativistic correction" label rows); Tables V and VI have CSV files.

Only the gamma column is taken (the alpha column is linear data, for the linear phase).
Row ids are "<basis> <method>", with " [core]" for the rows under an "... electrons
correlated" label (core correlation included) and " [relativistic]" for the row under
"Relativistic correction" (SCF with the mass-velocity and Darwin correction). The
footnote markers of the basis names (A^a, B^{b}, ...) and the LaTeX spacing are
stripped from the row ids; ``printed`` keeps the raw cells.

Table IV (neon) is reproduced from Taylor et al. (1989) except its MP2 row; it is
written to two files: the MP2 row (this paper's) and the rest (reproduced, excluded in
the record: those values were revised by an erratum).

Table VI: the theory column "gamma^a" (value and, in parentheses, the uncertainty) per
gas, one file per gas; the experimental column "gamma expt." per gas, where a cell can
hold two values with their footnote letters (Ne and Ar), split here.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Rice_1991_DeterminationAccurateDipole"
TABLES = {  # gas -> (table id, CSV available)
    "Ar": ("tbl_6788c2ece063", False),
    "Kr": ("tbl_7b5b5d1dda56", False),
    "Xe": ("tbl_d0c6f1d5dee5", False),
    "Ne": ("tbl_3fc8a4dd7268", False),
    "He": ("tbl_39e1e26d258b", True),
}
TABLE_VI = "tbl_ff5895bcf9bb"
_MARK = re.compile(r"\^\{?[a-d]\}?")


def _prefix(gas: str) -> str:
    return f"{gas}/rice_1991__"


def _basis(cell: str) -> str:
    return _MARK.sub("", cell).replace("\\ ", " ").replace(" ", "")


def _methods(gas: str) -> dict[str, list[Point]]:
    table_id, has_csv = TABLES[gas]
    grid = table(DIRECTORY, table_id) if has_csv else json_table(DIRECTORY, table_id)
    if grid[0] != ["Basis", "Method", "\\alpha", "\\gamma"]:
        msg = f"unexpected header in {table_id}: {grid[0]!r}"
        raise ValueError(msg)
    own: list[Point] = []
    reproduced: list[Point] = []
    suffix = ""
    for cells in grid[1:]:
        if cells[1] == "":  # a label row
            label = cells[0]
            if "electrons correlated" in label:
                suffix = " [core]"
            elif label == "Relativistic correction":
                suffix = " [relativistic]"
            else:
                msg = f"unknown label row {cells!r}"
                raise ValueError(msg)
            continue
        row = f"{_basis(cells[0])} {cells[1]}{suffix}"
        value = cells[3]
        float(value)
        point = Point(row, "inf", value, "", " | ".join(cells))
        # Table IV: every row but MP2 is Taylor et al.'s (footnote a).
        (reproduced if gas == "Ne" and cells[1] != "MP2" else own).append(point)
    rows = [p.row for p in own + reproduced]
    if len(set(rows)) != len(rows):
        msg = f"duplicate row ids in {table_id}: {rows!r}"
        raise ValueError(msg)
    files = {f"{_prefix(gas)}methods_{gas}.csv": own}
    if reproduced:
        files[f"{_prefix(gas)}methods_taylor_{gas}.csv"] = reproduced
    return files


_THEORY = re.compile(r"^([\d.]+)\((?:\\pm )?([\d.]+)\)(\^c)?$")
_EXPT = re.compile(r"([\d.]+)(?:\(\\pm ([\d.]+)\))?\^([dfg])")


def _table_vi() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_VI)
    if grid[0][3] != "\\gamma^{a}theory" or grid[0][4] != "\\gammaexpt.":
        msg = f"unexpected Table VI header: {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        gas = cells[0]
        printed = " | ".join(cells)
        theory = _THEORY.match(cells[3])
        if theory is None:
            msg = f"unparsed theory cell {cells[3]!r}"
            raise ValueError(msg)
        if theory.group(3):  # He: 43.6 (Bishop and Pipin's 43.104 in parentheses)
            point = Point(gas, "inf", theory.group(1), "", printed)
        else:
            point = Point(gas, "inf", theory.group(1), theory.group(2), printed)
        files[f"{_prefix(gas)}best_{gas}.csv"] = [point]
        expt = _EXPT.findall(cells[4])
        if cells[4] and "".join(f"{v}{s}{m}" for v, s, m in expt) == "":
            msg = f"unparsed expt cell {cells[4]!r}"
            raise ValueError(msg)
        for value, sigma, mark in expt:
            name = "expt_dc_kerr" if mark == "f" else "expt"
            path = f"{_prefix(gas)}{name}_{gas}.csv"
            wavelength = "632.8" if mark == "f" else "inf"
            files.setdefault(path, []).append(
                Point(f"{gas}^{mark}", wavelength, value, sigma, printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the gamma columns of Tables I-V and the gamma columns of Table VI."""
    files: dict[str, list[Point]] = {}
    for gas in TABLES:
        files |= _methods(gas)
    return files | _table_vi()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

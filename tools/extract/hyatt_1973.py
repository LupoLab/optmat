# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hyatt, Cherlow, Fenner and Porto (1973): the absolute N2 cross section.

- Table I (``tbl_fa72308ac503``, p. 2; the CSV cells match the image): the theoretical
  differential cross section of the H2 J = 1 -> 3 rotational line at 90 deg, the column
  dsigma/dOmega (cm^2/mol sr), at 488.0 and 514.5 nm. The page prints the 514.5 nm cell
  as ``1.79 x 19^-30``, a misprint of the column factor 10^-30 (the image shows it,
  and Eq. (1), dsigma/dOmega = 21 sigma_0/(80 pi), gives 1.788e-30 from the printed
  sigma_0): the mantissa is written with the dataset's scale 1e-30, the raw cell kept
  in the printed column. The sigma_0 and Phi(1) dsigma/dOmega columns are quoted in
  the record.
- Table II (``tbl_e6410a8fb5ef``, p. 2; no CSV, the JSON cells match the image): the
  absolute N2 Q-branch cross section at 488.0 and 514.5 nm, this work and six others'
  rows, ``...`` for an empty cell; the column factor ``\times 10^{-31}`` is the dataset
  scale. The footnote letters stuck to the row names are dropped from the row labels.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Hyatt_1973_CrossSectionRaman"
TABLE_I = "tbl_fa72308ac503"
TABLE_II = "tbl_e6410a8fb5ef"
WAVELENGTHS = ("488.0", "514.5")  # nm, the argon-ion lines of both tables
EMPTY = "..."
# The 514.5 nm cell of Table I prints the power's base as 19 (see the docstring).
_TABLE_I_CELL = re.compile(r"^(?P<value>[\d.]+)\\times (?:10|19)\^\{-30\}$")
_TABLE_II_CELL = re.compile(
    r"^\(?(?P<value>[\d.]+)(?: \\pm (?P<sigma>[\d.]+)\))? \\times 10\^\{-31\}$"
)
# Table II's rows -> (dataset, row label); the footnotes give the references.
TABLE_II_ROWS = {
    "This work": ("absolute_N2", "this work"),
    "Penney et al.a": ("absolute_N2_others", "Penney et al. (ref. 7)"),
    "Fouche and Changb": ("absolute_N2_others", "Fouche and Chang (ref. 8)"),
    "Fenner et al.c": ("absolute_N2_fenner", "Fenner et al., against benzene"),
    "Fenner et al.d": ("absolute_N2_fenner", "Fenner et al., calibrated spectrometer"),
    "Murphy et al.e,f": ("absolute_N2_others", "Murphy et al. (ref. 6)"),
    "Stansbury et al.e,g": ("absolute_N2_others", "Stansbury et al. (ref. 5)"),
}


def _table_i() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_I)
    points = []
    for cells, wavelength in zip(rows[1:], WAVELENGTHS, strict=True):
        match = _TABLE_I_CELL.match(cells[2])
        if cells[0] != wavelength or match is None:
            msg = f"Table I row {cells!r}"
            raise ValueError(msg)
        printed = " | ".join(cells)
        value = number(match["value"])
        points.append(Point(wavelength, wavelength, value, "", printed))
    return {"H2/hyatt_1973__theory_S1_H2.csv": points}


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    if tuple(grid[1][1:]) != ("488.0 nm", "514.5 nm"):
        msg = f"Table II header {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        dataset, label = TABLE_II_ROWS[cells[0]]
        printed = " | ".join(cells)
        for wavelength, cell in zip(WAVELENGTHS, cells[1:], strict=True):
            if cell == EMPTY:
                continue
            match = _TABLE_II_CELL.match(cell)
            if match is None:
                msg = f"Table II cell {cell!r}"
                raise ValueError(msg)
            files.setdefault(f"N2/hyatt_1973__{dataset}.csv", []).append(
                Point(
                    f"{label}, {wavelength} nm",
                    wavelength,
                    number(match["value"]),
                    match["sigma"] or "",
                    printed,
                )
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I and II as data files."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

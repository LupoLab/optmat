# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tammer & Hüttner (1990): dc Kerr of H2 and D2 (superseded calibration).

Both tables have no CSV rendering and are read from their JSON cells; both checked
against their images (no differences).

* Table 1 (``tbl_1f9c20083bf3``): per temperature (13 for H2, 196.1-411.7 K; one for
  D2, 298.1 K) the fits a, b of K against density, A_K, Delta = A_K - A_K^calc,
  B_K^(4), B_K^(2) and gamma, each with its standard deviation in the last digits in
  parentheses ("8.267(32)"), written here as a sigma in the value's unit. Taken: A_K
  (``table1_A_K_<gas>``), B_K^(4) and B_K^(2) (H2 only; D2 prints "-"), and gamma
  (with the temperature as an extra column). The gas label rows ("H_2", "D_2") and the
  footnote row are skipped. a, b and Delta are in the record's notes.
* Table 2 (``tbl_01dc04374704``): gamma in 10^-62 C^4 m^4 J^-3 for H2 and D2: this work
  (``table2_gamma_<gas>``) and the literature rows (``table2_literature_<gas>``). The
  esu column (the same values converted with 1 esu = 1.2379e-25 C^4 m^4 J^-3) is not
  taken. Footnote markers ("b)", "c)", "d)") are dropped.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Tammer_1990_AnisotropySecondHyperpolarizability"
TABLE_1 = "tbl_1f9c20083bf3"
TABLE_2 = "tbl_01dc04374704"
WAVELENGTH = "632.8"  # nm, the He-Ne line as printed (in air)
GASES = {"H_2": "H2", "D_2": "D2"}
_CELL = re.compile(r"^(-?[\d.]+)\((\d+)\)$")
_MARK = re.compile(r"(?:\^)?[bcd]\)$")
T1_COLUMNS = {
    "A_K": (4, "10^{29}A_K"),
    "B_K4": (6, "10^{34}B_{K}^{(4)}"),
    "B_K2": (7, "10^{34}B_{K}^{(2)}"),
    "gamma": (8, "10^{62}\\gamma"),
}


def _parenthesised(cell: str) -> tuple[str, str]:
    """Split "8.267(32)" into ("8.267", "0.032"): the sigma on the last digits."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value, digits = match[1], match[2]
    decimals = len(value.partition(".")[2])
    sigma = f"{int(digits) / 10**decimals:.{decimals}f}"
    return number(value), sigma


def _table1() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_1)
    for column, start in T1_COLUMNS.values():
        if not grid[0][column].startswith(start):
            msg = f"unexpected Table 1 header {grid[0][column]!r}"
            raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in grid[1:]:
        if cells[0] in GASES:
            gas = GASES[cells[0]]
            continue
        if cells[0].startswith("^a)"):
            continue  # the footnotes, in the record's notes
        temperature = number(cells[0])
        printed = " | ".join(cells)
        prefix = f"{gas}/tammer_1990__table1_"
        for name, (column, _) in T1_COLUMNS.items():
            if cells[column] == "-":
                continue  # D2: no B_K
            value, sigma = _parenthesised(cells[column])
            if name == "gamma":
                point = Point(
                    temperature,
                    WAVELENGTH,
                    value,
                    sigma,
                    printed,
                    {"temperature": temperature},
                )
            else:
                point = Point(temperature, temperature, value, sigma, printed)
            files.setdefault(f"{prefix}{name}_{gas}.csv", []).append(point)
    return files


def _table2() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_2)
    if not grid[0][1].startswith("10^{62}\\gamma"):
        msg = f"unexpected Table 2 header {grid[0][1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in grid[1:]:
        gas = GASES.get(cells[0], gas)
        value, sigma = _parenthesised(_MARK.sub("", cells[1].strip()))
        reference = _MARK.sub("", cells[3].strip()).strip("[] ")
        row = "this work" if reference == "this work" else f"ref {reference}"
        name = f"gamma_{gas}" if row == "this work" else f"literature_{gas}"
        point = Point(row, WAVELENGTH, value, sigma, " | ".join(cells))
        files.setdefault(f"{gas}/tammer_1990__table2_{name}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 1 and Table 2 files."""
    return _table1() | _table2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

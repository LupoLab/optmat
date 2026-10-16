# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tammer, Löblein, Peting & Hüttner (1992): calibrated dc Kerr of H2, He.

All tables checked against their images.

* Table 1 (``tbl_4d048dca076d``, CSV): the H2 dc Kerr data of Tammer & Hüttner (1990,
  their Table 1) re-evaluated with the new field integral (13 temperatures, 196.1-411.7
  K), and a last row measured in this work at 303.7 K. Every cell prints the standard
  deviation in the last digits in parentheses ("8.765(33)"), written here as a sigma in
  the value's unit. Taken: A_K (``table1_A_K_H2``, the 13 re-corrected rows, and
  ``table1_A_K_H2_this_work``, the last row), B_K^(4) and B_K^(2) (all 14 rows) and
  gamma (all 14 rows, with the temperature as an extra column). The temperature is the
  row id and the abscissa; the last row's id is "303.7 (this work)". The header of the
  sixth column prints a capital Delta (10^29 Delta, the residual A_K - A_K^calc of the
  1990 paper) which the extraction reads as "A": the column is not taken (structural;
  the image shows Delta). b', a and b (the fits of K against pressure and density) are
  in the record's notes.
* Table 2 (``tbl_05229ef25c1e``, JSON cells only): gamma(H2) in 10^-62 C^4 m^4 J^-3,
  this work (``table2_gamma_H2``) and the literature rows (``table2_literature_H2``).
  The footnote markers "^{a)}", "^{b)}" are dropped.
* Table 3 (``tbl_2ce09dc63f8c``, JSON cells only): gamma(He) in 10^-64 C^4 m^4 J^-3,
  this work (``table3_gamma_He``) and the literature rows (``table3_literature_He``).
  The footnote marker of the first row, "a)" on the image, is extracted as "^{n)}" and
  dropped.

The unlabelled table of K against pressure (``tbl_2b285a6d3881``) is in the record's
notes only; the text values (He A_K^(2), A_K^(4), K/p and the ratio of the gammas) are
inline points in the record.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Tammer_1992_FieldCalibratedMeasurements"
TABLE_1 = "tbl_4d048dca076d"
TABLE_2 = "tbl_05229ef25c1e"
TABLE_3 = "tbl_2ce09dc63f8c"
WAVELENGTH = "632.8"  # nm, the He-Ne line as printed (in air)
_CELL = re.compile(r"^(-?[\d.]+)\((\d+)\)$")
_MARK = re.compile(r"\^\{(?:a\)|b\)|\\text\{n\)\})\}$")
# Table 1 columns (0-based) and the header each must start with.
T1_COLUMNS = {
    "A_K": (5, "10^{29}A_K"),
    "B_K4": (7, "10^{34}B_{K}^{(4)}"),
    "B_K2": (8, "10^{34}B_{K}^{(2)}"),
    "gamma": (9, "10^{62}\\gamma"),
}
# Table 1's sixth column is 10^29 Delta on the image; the extraction reads "A".
DELTA_HEADER = "10^{29}A"


def _parenthesised(cell: str) -> tuple[str, str]:
    """Split "8.765(33)" into ("8.765", "0.033"): the sigma on the last digits."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value, digits = match[1], match[2]
    decimals = len(value.partition(".")[2])
    sigma = f"{int(digits) / 10**decimals:.{decimals}f}"
    return number(value), sigma


def _table1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_1)
    header = grid[0]
    if not header[6].startswith(DELTA_HEADER):
        msg = f"unexpected Table 1 header {header[6]!r}"
        raise ValueError(msg)
    for column, start in T1_COLUMNS.values():
        if not header[column].startswith(start):
            msg = f"unexpected Table 1 header {header[column]!r}"
            raise ValueError(msg)
    rows = grid[1:]
    if len(rows) != 14 or rows[-1][0] != "303.7":
        msg = "Table 1 should have 13 re-corrected rows and the new 303.7 K row"
        raise ValueError(msg)
    prefix = "H2/tammer_1992__table1_"
    files: dict[str, list[Point]] = {}
    for index, cells in enumerate(rows):
        temperature = number(cells[0])
        new = index == len(rows) - 1
        row = f"{temperature} (this work)" if new else temperature
        printed = " | ".join(cells)
        for name, (column, _) in T1_COLUMNS.items():
            value, sigma = _parenthesised(cells[column])
            if name == "gamma":
                point = Point(
                    row, WAVELENGTH, value, sigma, printed, {"temperature": temperature}
                )
                files.setdefault(f"{prefix}gamma_H2.csv", []).append(point)
                continue
            point = Point(row, temperature, value, sigma, printed)
            if name == "A_K":
                stem = "A_K_H2_this_work" if new else "A_K_H2"
            else:
                stem = f"{name}_H2"
            files.setdefault(f"{prefix}{stem}.csv", []).append(point)
    return files


def _comparison(
    table_id: str, gas: str, table_name: str, header: str
) -> dict[str, list[Point]]:
    """Return a Table 2 or 3 of gamma: this work, and the literature rows."""
    grid = json_table(DIRECTORY, table_id)
    if not grid[0][0].startswith(header):
        msg = f"unexpected header {grid[0][0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        text = _MARK.sub("", cells[0].strip())
        value, sigma = _parenthesised(text) if "(" in text else (number(text), "")
        reference = cells[1].strip("[] ")
        row = "this work" if reference == "this work" else f"ref {reference}"
        point = Point(row, WAVELENGTH, value, sigma, " | ".join(cells))
        suffix = f"gamma_{gas}" if row == "this work" else f"literature_{gas}"
        files.setdefault(f"{gas}/tammer_1992__{table_name}_{suffix}.csv", []).append(
            point
        )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 1, 2 and 3 files."""
    files = _table1()
    files |= _comparison(TABLE_2, "H2", "table2", "\\gamma (H_2)")
    files |= _comparison(TABLE_3, "He", "table3", "\\gamma(\\text{He})")
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

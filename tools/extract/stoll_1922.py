# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Stoll (1922): the dispersion of O2, N2 and CO2 at 4388-9224 A (German).

Every refractivity is ``(n0 - 1) x 10^7`` at 0 C and 760 mm Hg, printed with a decimal
comma (``2755,09``), read here as a decimal point. The gas tables label each line by
its rounded air wavelength (``4388``, ``8264`` or ``8265``); the full wavelengths of the
same He and Ar lines are printed in Table II (air, ``tbl_27924830a3f0``: ``4388,100``,
``8264,521``), and each label is given the Table II wavelength within 1 A of it. The
row id is the label as printed. Table II's 5015 A cell has a damaged digit (face
value 5015,722; the library reads 5015,732); the other air table prints 5015,732
clearly, which is used (the difference is immaterial).

- Table V (``tbl_6ed4c184eeb4``), electrolytic oxygen; Table VII (``tbl_af469b918b2b``),
  nitrogen from the gasometer; Table IX (``tbl_098034e2fc60``), carbon dioxide: the
  mean of the single determinations (``Mittel``) and the value from the author's
  Cauchy formula (``berechnet``). The single determinations, the differences and the
  columns computed from Cuthbertson's and Scheel's formulas are not extracted.
- Tables VI (``tbl_1b1b829fe025``, oxygen from a cylinder, about 4 % N2) and VIII
  (``tbl_91f9a85c75b4``, nitrogen from a cylinder, about 0.5 % O2): means only, printed
  in two column pairs.
- The unlabelled CO2 table (``tbl_9e61e81d2f8d``, cells file only): Koch's four values;
  the 'Stoll' column is the author's formula at Koch's wavelengths (not extracted).

Images checked. Two cells are misread by the library (``3`` read as ``8``); the images
print Table IX's calculated value at 6678 A as ``4473,31`` (the cells ``4473,91``; the
printed difference -0.18 = 4473.13 - 4473.31 confirms it) and Table VIII's mean at
7065 A as ``2951,32`` (the cells ``2951,82``). Both are corrected here.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Stoll_1922_DieDispersionDer"
TABLE_II = "tbl_27924830a3f0"
# Table -> (gas, dataset prefix, mean column, calculated column or None).
MAIN = {
    "tbl_6ed4c184eeb4": ("O2", "electrolytic", 7, 8),
    "tbl_af469b918b2b": ("N2", "gasometer", 4, 5),
    "tbl_098034e2fc60": ("CO2", "cylinder", 5, 6),
}
PAIRS = {
    "tbl_1b1b829fe025": ("O2", "cylinder"),
    "tbl_91f9a85c75b4": ("N2", "cylinder"),
}
KOCH = "tbl_9e61e81d2f8d"
# (table, row label, column): (cell as read by the library, cell as printed).
MISREAD = {
    ("tbl_098034e2fc60", "6678", 6): ("4473,91", "4473,31"),
    ("tbl_91f9a85c75b4", "7065", 1): ("2951,82", "2951,32"),
}


def _decimal(cell: str) -> str:
    """Read a printed decimal comma as a decimal point."""
    return number(cell.replace(",", "."))


def _cell(table_id: str, label: str, column: int, text: str) -> str:
    """Return a cell, corrected where the library misread the printed digit."""
    key = (table_id, label, column)
    if key in MISREAD:
        read, printed = MISREAD[key]
        if text != read:
            msg = f"{key}: expected the misread {read!r}, found {text!r}"
            raise ValueError(msg)
        return printed
    return text


def _wavelengths() -> list[str]:
    """Return Table II's full line wavelengths (A, air), as printed."""
    return [_decimal(cells[0]) for cells in table(DIRECTORY, TABLE_II)[1:]]


def _full(label: str, lines: list[str]) -> str:
    """Return the Table II wavelength of the line a rounded label names."""
    (match,) = [w for w in lines if abs(float(w) - float(label)) < 1.0]
    return match


def extract() -> dict[str, list[Point]]:
    """Return the O2, N2 and CO2 tables and Koch's CO2 values as data files."""
    lines = _wavelengths()
    files: dict[str, list[Point]] = {}
    for table_id, (gas, prefix, mean, calculated) in MAIN.items():
        for cells in table(DIRECTORY, table_id)[1:]:
            label = cells[0].strip()
            wavelength = _full(label, lines)
            printed = " | ".join(cells)
            for column, kind in ((mean, prefix), (calculated, f"fit_{prefix}")):
                value = _decimal(_cell(table_id, label, column, cells[column]))
                files.setdefault(f"{gas}/stoll_1922__{kind}_{gas}.csv", []).append(
                    Point(label, wavelength, value, "", printed)
                )
    for table_id, (gas, prefix) in PAIRS.items():
        points = []
        for cells in table(DIRECTORY, table_id)[1:]:
            for label_column, value_column in ((0, 1), (2, 3)):
                label = cells[label_column].strip()
                if not label:
                    continue
                text = _cell(table_id, label, value_column, cells[value_column])
                points.append(
                    Point(
                        label,
                        _full(label, lines),
                        _decimal(text),
                        "",
                        " | ".join(cells),
                    )
                )
        points.sort(key=lambda p: float(p.x))
        files[f"{gas}/stoll_1922__{prefix}_{gas}.csv"] = points
    files["CO2/stoll_1922__koch_CO2.csv"] = [
        Point(cells[0], cells[0], number(cells[1]), "", " | ".join(cells))
        for cells in json_table(DIRECTORY, KOCH)[2:]
    ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

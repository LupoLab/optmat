# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Matsuoka et al. (1991): the absolute index of silica glass, -165 to 83 C.

Table 1 (``tbl_3e381ac738a4``, cells only) gives the absolute index at ten photon
energies (eV) and seven temperatures (-165.4, -111.3, -43.4, 2.5, 20.5, 45.2 and 83.3
degrees C). One file per temperature; rows are the photon energies as printed.

Structural fix (checked on the table image, pdf p. 2): the scan prints the decimal point
of many cells so faintly that the library reads it as a space ("1 459241" for 1.459241,
"2 270486" for 2.270486, "-165 4" for -165.4); the other cells of the same columns
carry the point, and every value has one digit before it. The script restores the point
where a cell is two digit groups separated by one space, and refuses any other shape.
At -165.4 and -111.3 C the rows 4.504003 and 5.213151 eV print "-" (not measured): no
point is written for them.

Table 2 (``tbl_c69345cf4e41``) holds the parameters of the temperature-dependent
Sellmeier equation, the record's formula; this script checks the record's strings
against the cells, with the same decimal-point fix.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Matsuoka_1991_TemperatureDependenceRefractive"
TABLE_1 = "tbl_3e381ac738a4"
TABLE_2 = "tbl_c69345cf4e41"
PREFIX = "SiO2/matsuoka_1991__"

TEMPERATURES = ("-165.4", "-111.3", "-43.4", "2.5", "20.5", "45.2", "83.3")
NAMES = ("m165.4C", "m111.3C", "m43.4C", "2.5C", "20.5C", "45.2C", "83.3C")
_LOST_POINT = re.compile(r"^(-?\d+) (\d+)$")
_NUMBER = re.compile(r"^-?\d+\.\d+$")

# Table 2 as printed (the image), with the decimal points the scan loses restored.
TABLE_2_PRINTED = (
    ("228.7018", "46.40806", "0.014173"),
    ("4.93 \\times 10^{-5}", "-3.27 \\times 10^{-5}", "-1.704 \\times 10^{-6}"),
    ("1.10 \\times 10^{-7}", "-3.78 \\times 10^{-8}", "-2.14 \\times 10^{-9}"),
    ("18.111630", "10.671082", "0.125"),
    ("9.15 \\times 10^{-6}", "-2.9913 \\times 10^{-4}", "<10^{-5}"),
    ("7.478 \\times 10^{-8}", "-4.8074 \\times 10^{-7}", "<10^{-8}"),
)


def _point(cell: str) -> str:
    """Return a cell with its lost decimal point restored (see the docstring)."""
    match = _LOST_POINT.match(cell)
    text = f"{match.group(1)}.{match.group(2)}" if match else cell
    if not _NUMBER.match(text):
        msg = f"unexpected Table 1 cell {cell!r}"
        raise ValueError(msg)
    return text


def _check_table_2() -> None:
    """Check the record's transcription of Table 2 against the library cells."""
    rows = table(DIRECTORY, TABLE_2)[1:]
    # The lost point is read as a thin space ("4\,93") or a space ("18 111630").
    found = tuple(
        tuple(_LOST_POINT.sub(r"\1.\2", cell.replace("\\,", ".")) for cell in r[1:])
        for r in rows
    )
    if found != TABLE_2_PRINTED:
        msg = f"Table 2 reads {found}, not {TABLE_2_PRINTED}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return the absolute index at each of the seven temperatures."""
    _check_table_2()
    rows = json_table(DIRECTORY, TABLE_1)
    temperatures = tuple(_point(c.removesuffix("°C")) for c in rows[1][1:])
    if temperatures != TEMPERATURES or len(rows) != 12:
        msg = f"unexpected Table 1 layout {rows[1]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        energy = number(_point(cells[0]))
        printed = " | ".join(cells)
        for name, cell in zip(NAMES, cells[1:], strict=True):
            if cell == "-":  # not measured at this temperature
                continue
            point = Point(energy, energy, number(_point(cell)), "", printed)
            files.setdefault(f"{PREFIX}n_{name}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

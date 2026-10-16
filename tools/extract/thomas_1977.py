# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Thomas & Meath (1977): the CH4 dipole oscillator strengths and sums.

The tables print a middle dot as the decimal point (``0·273``) and powers of ten in
parentheses after a space (``7·299 (-3)`` = 7.299e-3); the dot is read as a point and
the number written in e-notation with its digits unchanged (the ``printed`` column
keeps the cell). Images checked against the cells (2026-09-27).

- Table 1 (``tbl_df5a5473a0a4``, CSV): the recommended (constrained) DOSD integrated
  over 17 spectral regions, from the absorption threshold (145.86 nm) to 0.00 nm, as
  ``interval_oscillator_strength`` against wavelength (nm): x the region's long- and
  x_upper its short-wavelength edge, the row the region as printed. Only the
  "Recommended" column is extracted; the reference number in brackets after a value
  and the marks (c), (d) are stripped (the record says what they mean). The
  "unmodified experimental data" column (several references per cell) is not.
- Table 2 (``tbl_010d4e4faeeb``, CSV): S(-1), S(-2) and S(-4) of the recommended DOSD
  and of other sources, one ``cauchy_moment`` file per row, rows by k.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Thomas_1977_DipoleSpectrumSums"
TABLE_1 = "tbl_df5a5473a0a4"
TABLE_2 = "tbl_010d4e4faeeb"
PREFIX = "CH4/thomas_1977__"
TABLE_2_ROWS = {
    "Platzman [80]": "platzman_1962",
    "Santar and Bednář [79]": "santar_1969",
    "Watson and Ramasawamy [69]": "fit_watson_1936",
    "Cuthbertson and Cuthbertson [75]": "fit_cuthbertson_1920",
    "Friberg [70]": "fit_friberg_1927",
    "Frivold et al. [71]": "fit_frivold_1939",
    "Calculated using (4) and the recommended DOSD": "moments",
}
ORDERS = {1: "-1", 2: "-2", 3: "-4"}
_VALUE = re.compile(
    r"^(?P<mantissa>\d+\.\d+)(?:\s*\((?P<exponent>-?\d)\))?"
    r"(?:\s*\[\d+\]|\s*\([cd]\))?$"
)


def _value(cell: str) -> str:
    """Return a cell's number in e-notation, the digits unchanged."""
    text = cell.strip().replace("\N{MIDDLE DOT}", ".").replace("\N{MINUS SIGN}", "-")
    found = _VALUE.match(text)
    if found is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = found["mantissa"], found["exponent"]
    return number(f"{mantissa}e{exponent}" if exponent else mantissa)


def _table_1() -> list[Point]:
    points = []
    for cells in table(DIRECTORY, TABLE_1)[1:]:
        region = cells[0].strip().replace("\N{MIDDLE DOT}", ".").replace(" ", "")
        long, short = region.split("-")
        printed = " | ".join(cells)
        point = Point(
            region, number(long), _value(cells[1]), "", printed, {"x_upper": short}
        )
        points.append(point)
    return points


def _table_2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        stem = TABLE_2_ROWS[cells[0].strip()]
        printed = " | ".join(cells)
        for column, order in ORDERS.items():
            if cells[column].strip():
                point = Point(order, order, _value(cells[column]), "", printed)
                files.setdefault(f"{PREFIX}{stem}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the recommended DOSD (Table 1) and the sums of Table 2."""
    return {f"{PREFIX}dosd.csv": _table_1()} | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

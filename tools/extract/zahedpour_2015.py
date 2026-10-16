# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zahedpour, Wahlstrand & Milchberg (2015): Table 1, n2 of N2, O2, Ar.

Table 1 is on p. 5797. Rows are pump wavelengths 400-2400 nm. Columns (all 1e-20 cm2/W,
"scaled to atmospheric pressure"): N2 with the Delta alpha dispersion adjustment (note
b): this work, [32], n2,rot; N2 without it (note c): this work, n2,rot; O2 (note c):
this work, [32], n2,rot; Ar (note d, referenced to N2): this work, [32]. Measured cells
print ``value ± error``. The library has no CSV; the JSON cells match the table image
(checked).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Zahedpour_2015_MeasurementNonlinearRefractive"
TABLE = "tbl_818914a52f6f"
COLUMNS = {
    1: ("N2", "n2_adjusted"),
    2: ("N2", "shelton"),
    3: ("N2", "n2_rot_adjusted"),
    4: ("N2", "n2"),
    5: ("N2", "n2_rot"),
    6: ("O2", "n2"),
    7: ("O2", "shelton"),
    8: ("O2", "n2_rot"),
    9: ("Ar", "n2"),
    10: ("Ar", "shelton"),
}


def _value_error(cell: str) -> tuple[str, str]:
    """Split ``value ± error`` (or a bare value) into the printed numbers."""
    if "±" in cell:
        value, error = cell.split("±")
        return number(value), number(error)
    return number(cell), ""


def extract() -> dict[str, list]:
    """Return one data file per column, rows by pump wavelength."""
    files: dict[str, list] = {}
    for cells in json_table(DIRECTORY, TABLE)[3:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, (gas, name) in COLUMNS.items():
            value, error = _value_error(cells[column])
            path = f"{gas}/zahedpour_2015__{name}_{gas}.csv"
            files.setdefault(path, []).append(
                Point(wavelength, wavelength, value, error, printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

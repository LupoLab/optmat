# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ityaksov, Linnartz and Ubachs (2008): deep-UV Rayleigh cross sections.

Table 1 (``tbl_0ec99fd2f68e``, cells file only; image checked) gives, at eight vacuum
wavelengths from 197.70 to 270.15 nm, the measured extinction cross section of N2, CH4
and SF6 (10^-25 cm^2, +- one sigma) and the value predicted from refractive-index and
King-factor formulas (+- from those formulas; none for SF6). A dash marks a wavelength
not measured for that gas; it is skipped. The three header rows (the unit, the gases,
Measured/Predicted) are dropped.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Ityaksov_2008_DeepUVRayleigh"
TABLE_1 = "tbl_0ec99fd2f68e"
# Gas -> (measured column, predicted column).
COLUMNS = {"N2": (1, 2), "CH4": (3, 4), "SF6": (5, 6)}
_VALUE = re.compile(r"^\s*([\d.]+)\s*(?:±\s*([\d.]+))?\s*$")


def _value(cell: str) -> tuple[str, str] | None:
    """Return a cell's value and +- (or ""), or None for a dash."""
    if cell.strip() in {"-", "\N{EN DASH}", "\N{EM DASH}"}:
        return None
    match = _VALUE.match(cell)
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2) or ""


def extract() -> dict[str, list[Point]]:
    """Return Table 1's measured and predicted cross sections as data files."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_1)[3:]:
        wavelength = cells[0].strip()
        printed = " | ".join(cells)
        for gas, (measured, predicted) in COLUMNS.items():
            for column, kind in ((measured, "rayleigh"), (predicted, "predicted")):
                value = _value(cells[column])
                if value is None:
                    continue
                files.setdefault(f"{gas}/ityaksov_2008__{kind}_{gas}.csv", []).append(
                    Point(wavelength, wavelength, value[0], value[1], printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

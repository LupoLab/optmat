# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bree, Demircan & Steinmeyer (2012): noble-gas n2 at five laser wavelengths.

Table I (``tbl_934361b8054e``) prints n2 in cm^2/TW as ``m(e)``, "The number in
brackets represents the power of 10"; the value is written as ``m``e``e`` (the form
the WIAS preprint of the paper prints, ``Unverified_df851ddaf7be`` Table 1) and the
dataset's scale turns cm^2/TW into cm^2/W. The table image equals the CSV. Rows are
the laser lines 193, 248, 400, 800 and 2000 nm, labelled with their laser in the first
cell.
"""

import re

from _library import MINUS, Point, table, write

DIRECTORY = "Bree_2012_KramersKronigRelations"
TABLE = "tbl_934361b8054e"
_CELL = re.compile(rf"^\s*([\d.]+)\(([-{MINUS}]\d+)\)\s*$")
_ROW = re.compile(r"^(\d+) nm \(")


def _power_of_ten(cell: str) -> str:
    """Turn a printed ``7.8(-9)`` into ``7.8e-9``."""
    match = _CELL.match(cell)
    assert match, cell
    mantissa, exponent = match.groups()
    return f"{mantissa}e{exponent.replace(MINUS, '-')}"


def extract() -> dict[str, list[Point]]:
    """Return one file per gas with the five wavelengths."""
    header, *rows = table(DIRECTORY, TABLE)
    gases = header[1:]
    assert gases == ["He", "Ne", "Ar", "Kr", "Xe"], header
    files: dict[str, list[Point]] = {gas: [] for gas in gases}
    for cells in rows:
        match = _ROW.match(cells[0])
        assert match, cells
        wavelength = match.group(1)
        for gas, cell in zip(gases, cells[1:], strict=True):
            files[gas].append(
                Point(
                    wavelength, wavelength, _power_of_ten(cell), "", " | ".join(cells)
                )
            )
    return {f"{gas}/bree_2012__kk_{gas}.csv": points for gas, points in files.items()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

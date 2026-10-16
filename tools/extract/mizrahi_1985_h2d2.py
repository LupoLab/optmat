# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mizrahi & Shelton (1985, PRA): chi(3) ratios H2/He, H2/D2 (Tables I, II).

Table I (``tbl_d9c17b1ffdb1``, JSON cells only; image checked) gives, per laser line
(air wavelength, vacuum wavenumber), chi_H2/chi_He at 19 lines and chi_H2/chi_D2 at 9,
each as raw data and after the calculated vibrational-resonance correction, with 'the
uncertainty of the last digits in parentheses' ("17.53(11)" = 17.53 +- 0.11). The first
row (lambda "infinity", nu = 0) holds the extrapolated fit values A; it is not a
measurement and is recorded as the formulas' A in the source record, not in a data
file. One file per column: ``H2_over_He_raw``, ``H2_over_He_corrected``,
``H2_over_D2_raw``, ``H2_over_D2_corrected``.

Table II (``tbl_79b8cb3720a5``, CSV; image checked) prints, per molecule, the
Hunt-Poll-Wolniewicz transition matrix elements and the static vibrational
contribution "\chi^{(3)vib}_{H_2}(0)=184", "...D_2(0)=178" (au): these two go to
``vibrational_static_<gas>`` (the number after "="); the matrix elements stay in the
record's notes.
"""

import re
from decimal import Decimal

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Mizrahi_1985_NonlinearSusceptibilityH2"
TABLE1 = "tbl_d9c17b1ffdb1"
TABLE2 = "tbl_79b8cb3720a5"
COLUMNS = {
    2: "H2_over_He_raw",
    3: "H2_over_He_corrected",
    4: "H2_over_D2_raw",
    5: "H2_over_D2_corrected",
}
_PARENTHESES = re.compile(r"^(\d+\.(\d+))\((\d+)\)$")
_VIBRATIONAL = re.compile(
    r"^\\chi\^\{\(3\)\\text\{vib\}\}_\{([HD])_\{2\}\}\(0\)=(\d+)$"
)


def parenthesised(cell: str) -> tuple[str, str]:
    """Return the value and sigma of a cell printed as "17.53(11)"."""
    match = _PARENTHESES.match(cell.strip())
    if match is None:
        msg = f"not a value(uncertainty) cell: {cell!r}"
        raise ValueError(msg)
    value, decimals, digits = match.groups()
    return value, f"{Decimal(digits).scaleb(-len(decimals)):f}"


def extract() -> dict[str, list[Point]]:
    """Return the four ratio columns and the two static vibrational values."""
    grid = json_table(DIRECTORY, TABLE1)
    if grid[1][2:] != ["Raw data", "Vibr.corrected", "Raw data", "Vibr.corrected"]:
        msg = f"Table I headings changed: {grid[1]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[3:]:  # grid[2] is the extrapolated omega = 0 row
        label, wavenumber = number(cells[0]), number(cells[1])
        for column, name in COLUMNS.items():
            if not cells[column]:
                continue
            value, sigma = parenthesised(cells[column])
            files.setdefault(f"H2/mizrahi_1985_h2d2__{name}.csv", []).append(
                Point(label, wavenumber, value, sigma, " | ".join(cells))
            )
    for cell in table(DIRECTORY, TABLE2)[-1]:
        match = _VIBRATIONAL.match(cell.strip())
        if match is None:
            msg = f"unexpected Table II cell {cell!r}"
            raise ValueError(msg)
        gas = f"{match.group(1)}2"
        files[f"{gas}/mizrahi_1985_h2d2__vibrational_static_{gas}.csv"] = [
            Point("static", "inf", number(match.group(2)), "", cell)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

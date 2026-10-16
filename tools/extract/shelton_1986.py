# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Shelton & Mizrahi (1986): D2/H2 ESHG phase-match density ratios (Table III).

Table III (``tbl_0a217d3795df``, checked against its image) gives, at nine fundamental
wavelengths, the measured ratio of the densities at which D2 and H2 periodically
phase-match electric-field-induced second-harmonic generation, printed with the error
of the last digit in parentheses (``1.0308(7)``), and the ratio Shelton calculated from
Rychlewski's ab initio polarizabilities. The parenthesised error is written to the
``sigma`` column as the number it stands for (``(7)`` on ``1.0308`` is 0.0007).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Shelton_1986_RefractiveIndexDispersion"
TABLE = "tbl_0a217d3795df"

_PARENTHESISED = re.compile(r"^\s*(\d+)\.(\d+)\((\d+)\)\s*$")


def _value_sigma(cell: str) -> tuple[str, str]:
    """Split ``"1.0308(7)"`` into ``"1.0308"`` and ``"0.0007"``."""
    match = _PARENTHESISED.match(cell)
    if match is None:
        msg = f"not a value with a parenthesised error: {cell!r}"
        raise ValueError(msg)
    whole, fraction, error = match.groups()
    decimals = len(fraction)
    sigma = f"{int(error) * 10.0**-decimals:.{decimals}f}"
    return f"{whole}.{fraction}", sigma


def extract() -> dict[str, list[Point]]:
    """Return the measured and the calculated ratio columns."""
    measured = []
    calculated = []
    for cells in table(DIRECTORY, TABLE)[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        value, sigma = _value_sigma(cells[1])
        measured.append(Point(wavelength, wavelength, value, sigma, printed))
        calculated.append(Point(wavelength, wavelength, number(cells[2]), "", printed))
    return {
        "D2/shelton_1986__ratio_D2_H2.csv": measured,
        "D2/shelton_1986__rychlewski_ratio_D2_H2.csv": calculated,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

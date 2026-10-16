# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hohm & Kerl (1991): the polarizability of CH4 at 632.99 nm, 128-885 K.

Table 2 (``tbl_ca1d21fe5b0e``, CSV; image checked, no differences) gives per
temperature T/K the mean polarizability 10^30 alpha_0(lambda, T)/m^3 at 632.99 nm and
-B(T)/(cm^3 mol^-1), each with its standard deviation on the last digits in
parentheses ("2.6151(44)" is 2.6151 +- 0.0044, p. 36). The header's alpha extracts as
"x_0"; the image prints alpha_0 (structural: the header is not used). Only the
polarizability column is extracted, one file with the temperature as the row id (as
printed, without its own uncertainty); B(T), the thermodynamic second virial
coefficient, has no quantity and is transcribed in the record's notes.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Hohm_1991_InterferometricDeterminationMean"
TABLE_2 = "tbl_ca1d21fe5b0e"
_UNCERTAIN = re.compile(r"^([\d.]+)\((\d+)\)$")


def _split(cell: str) -> tuple[str, str]:
    """Return "2.6151(44)" as ("2.6151", "0.0044")."""
    match = _UNCERTAIN.match(cell.strip())
    if match is None:
        msg = f"no uncertainty in {cell!r}"
        raise ValueError(msg)
    value, digits = match[1], match[2]
    decimals = len(value.partition(".")[2])
    return number(value), f"{int(digits) / 10**decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return Table 2's polarizabilities."""
    grid = table(DIRECTORY, TABLE_2)
    if grid[0][0] != "T/K" or "(\\lambda, T)/m^3" not in grid[0][1]:
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[1:]:
        temperature, _ = _split(cells[0])
        value, sigma = _split(cells[1])
        points.append(Point(temperature, "632.99", value, sigma, " | ".join(cells)))
    return {"CH4/hohm_1991__alpha_T_CH4.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

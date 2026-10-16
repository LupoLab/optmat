# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Brodersen (1979): Table 2.3, the constants of 14N2 from Raman spectra.

Chapter 2 (pp. 7-70) of *Raman Spectroscopy of Gases and Liquids*, A. Weber, ed.,
Topics in Current Physics 11 (Springer, 1979); the library's copy is the whole book
(``Weber_1979_RamanSpectroscopyGases``, formerly ``Unverified_6d54463e97c7``; PDF
page n is book page n - 12 here).

Table 2.3 (``tbl_ec367a2c12bf``, p. 41; the CSV matches the image): 'Constants for
isotopic species of nitrogen, determined from Raman spectra [2.62]' (Bendtsen 1974).
Only the 14N2 row is taken (B_0, D_0, nu_0, alpha); 14N15N, 15N2 and r_e (no quantity)
are in the record's notes. 'Standard errors are given in parentheses in units of the
last decimal. nu_0 values are accurate to +-0.01 cm-1' (the table note): the
parenthesised digits are converted to the printed unit by the value's decimals
('1.989574(12)' is 0.000012), and nu_0's +-0.01 is the dataset's uncertainty.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Weber_1979_RamanSpectroscopyGases"
TABLE_2_3 = "tbl_ec367a2c12bf"
PREFIX = "N2/brodersen_1979__"
_CELL = re.compile(r"^(?P<value>[\d.]+)\s*\((?P<digits>\d+)\)(?:\\times 10\^\{-6\})?$")


def _with_sigma(cell: str) -> tuple[str, str]:
    """Return the value and the standard error of a '1.989574(12)' cell."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"not a 'value(n)' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group("value"))
    decimals = len(value.partition(".")[2])
    return value, f"{int(match.group('digits')) * 10**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the 14N2 files of Table 2.3."""
    grid = table(DIRECTORY, TABLE_2_3)
    header = [
        "",
        "B_0[cm^{-1}]",
        "D_0[cm^{-1}]",
        "v_0[cm^{-1}]",
        "\\alpha[cm^{-1}]",
        "r_e[Å]",
    ]
    if grid[0] != header or grid[1][0] != "^{14}N_2":
        msg = f"unexpected Table 2.3 layout {grid[:2]!r}"
        raise ValueError(msg)
    cells = grid[1]
    if not cells[2].endswith("\\times 10^{-6}"):
        msg = f"D_0 is not x 10^-6: {cells[2]!r}"
        raise ValueError(msg)
    printed = " | ".join(cells)
    files = {}
    for name, column in (("B0", 1), ("D0", 2), ("alpha", 4)):
        value, sigma = _with_sigma(cells[column])
        files[f"{PREFIX}table2_3_{name}_N2.csv"] = [
            Point("14N2", "14N2", value, sigma, printed)
        ]
    files[f"{PREFIX}table2_3_nu0_N2.csv"] = [
        Point("14N2", "14N2", number(cells[3]), "", printed)
    ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zeiss, Meath, MacDonald & Dawson (1977): dipole sums S_k of H2, N2, O2, N2O.

Table 3 (``tbl_81ed6f008ff0``, JSON only) lists S_k for k = 2, 1, 0, -1/2, -1, -3/2,
-2, -5/2, -3, -4, -6, ..., -14 (then L_k and I_k), atomic units, with powers of ten in
parentheses (``1.230(4)`` is 1.230e4, written in e-notation; the printed cell is kept).
Rows are named by k in decimals (-1/2 is ``-0.5``). Only the curated species are taken
(Li, N, O, NH3, H2O and NO are skipped); L_k and I_k have no quantity.

Structural fix, from the table image: the library's table has lost its last column,
N2O (it ends with NO; the image has ten species). The N2O values are read from the
image into ``N2O_FROM_IMAGE``; the script asserts that the rows they belong to are the
library's k rows. They agree with Margoliash & Meath's (1978) N2O pseudospectrum, which
was fitted to these sums, to within the printed four figures (source notes). The library
also misreads Li's S_-4 ("3.504(3)"; the image prints 3.504(4)); Li is not extracted.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Zeiss_1977_DipoleOscillatorStrength"
TABLE_3 = "tbl_81ed6f008ff0"
COLUMNS = {"H2": 5, "N2": 6, "O2": 7}
# Table 3's N2O column (S_k rows, k = 2 ... -14), read from the table image.
N2O_FROM_IMAGE = (
    ("2", "2.291(4)"),
    ("1", "2.313(2)"),
    ("0", "2.200(1)"),
    ("-1/2", "1.561(1)"),
    ("-1", "1.481(1)"),
    ("-3/2", "1.630(1)"),
    ("-2", "1.970(1)"),
    ("-5/2", "2.547(1)"),
    ("-3", "3.469(1)"),
    ("-4", "7.211(1)"),
    ("-6", "4.107(2)"),
    ("-8", "2.847(3)"),
    ("-10", "2.211(4)"),
    ("-12", "1.887(5)"),
    ("-14", "1.793(6)"),
)
_CELL = re.compile(r"^(?P<mantissa>-?\d+(?:\.\d*)?)(?:\((?P<exponent>-?\d+)\))?$")


def _value(cell: str) -> str:
    """Return a printed ``mantissa(power of ten)`` cell's number in e-notation."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = match["mantissa"], match["exponent"]
    return f"{mantissa}e{exponent}" if exponent else mantissa


def _order(cell: str) -> str:
    """Return k as a decimal string: "-1/2" -> "-0.5"."""
    numerator, _, denominator = cell.strip().partition("/")
    return str(int(numerator) / int(denominator)) if denominator else numerator


def extract() -> dict[str, list[Point]]:
    """Return S_k for H2, N2, O2 (library cells) and N2O (read from the image)."""
    grid = json_table(DIRECTORY, TABLE_3)
    rows = grid[2:17]  # the S_k block, k = 2 ... -14
    assert [c[1] for c in rows] == [k for k, _ in N2O_FROM_IMAGE]
    files = {}
    for gas, column in COLUMNS.items():
        files[f"{gas}/zeiss_1977__moments_{gas}.csv"] = [
            Point(_order(c[1]), _order(c[1]), _value(c[column]), "", " | ".join(c))
            for c in rows
        ]
    files["N2O/zeiss_1977__moments_N2O.csv"] = [
        Point(_order(k), _order(k), _value(v), "", f"{k} | N2O (table image): {v}")
        for k, v in N2O_FROM_IMAGE
    ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

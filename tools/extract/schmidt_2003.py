# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Schmidt & Moldover (2003): A_eps,273 of Ar, N2, O2, CH4 and CO2 (Table V).

Table V (``tbl_d4ad17b1acdf``, CSV; checked against the image) lists, per gas and
equation of state, the Eq. (9) coefficients A_eps,273, b, c, A_tau and q with the
uncertainty of their last digits in parentheses "for the more satisfactory fits", and
the standard deviation 10^6 x sigma of the fit. Structural problems of the CSV (the
image has eight columns, Gas, Eq. of state ref., A_eps,273, b, c, A_tau, q, 10^6 x
sigma):

- the b and c columns are merged into one cell ("1.250 (34) -44.7 (68)"), and so are
  q and sigma ("114.21 (61) 1.3");
- the headers are garbled ("b c cm3·mol-1 cm6·mol-2", "Ay" for A_tau, "106 s" for
  10^6 x sigma).

Only the A_eps,273 column is extracted here (the full fits are formulas in the record);
it is the third cell of each row and is not touched by the merges. The in-scope rows
are Ar, the two N2, two O2, two CH4 and two CO2 fits (propane and ethane are not
curated species). Each fit is one file, with the fitted parameters as the row id; the
parenthesised digits become the sigma column ("4.14203 (19)" -> 0.00019).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Schmidt_2003_DielectricPermittivityEight"
TABLE_V = "tbl_d4ad17b1acdf"
_ENTRY = re.compile(r"^(?P<value>-?\d+\.\d*)(?: \((?P<digits>\d+)\))?$")

# (row index in the CSV, gas, dataset suffix, fitted parameters): the image's rows.
ROWS = (
    (1, "Ar", "Ar", "A_eps/b"),
    (2, "N2", "two_N2", "A_eps/b"),
    (3, "N2", "tau_N2", "A_eps/b/A_tau"),
    (4, "O2", "two_O2", "A_eps/b"),
    (5, "O2", "tau_O2", "A_eps/b/A_tau"),
    (6, "CH4", "two_CH4", "A_eps/b"),
    (7, "CH4", "four_CH4", "A_eps/b/c/A_tau"),
    (10, "CO2", "three_CO2", "A_eps/b/A_tau"),
    (11, "CO2", "four_CO2", "A_eps/b/c/q"),
)


def last_digits(value: str, digits: str) -> str:
    """Return the uncertainty "(19)" of the last digits of ``value`` as a number."""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) / 10**decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file per in-scope Table V fit: its A_eps,273."""
    rows = table(DIRECTORY, TABLE_V)
    files: dict[str, list[Point]] = {}
    for index, gas, suffix, fitted in ROWS:
        cells = rows[index]
        if cells[0] != gas:
            msg = f"row {index} is {cells[0]!r}, not {gas}"
            raise ValueError(msg)
        match = _ENTRY.match(cells[2].strip())
        if match is None:
            msg = f"unreadable A_eps,273 cell {cells[2]!r}"
            raise ValueError(msg)
        value = number(match["value"])
        sigma = last_digits(value, match["digits"]) if match["digits"] else ""
        path = f"{gas}/schmidt_2003__A_eps_{suffix}.csv"
        files[path] = [Point(fitted, "inf", value, sigma, " | ".join(cells))]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

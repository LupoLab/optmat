# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Yu, Miller, Drouin and Mueller (2012): the O2 X-state constants.

- Table III (``tbl_844a6de53312``, p. 8; no CSV, the JSON cells match the image): the
  16O16O X 3Sigma_g- constants (MHz) of v = 0 and 1 from two band-by-band fits, Fit #1
  and Fit #2, one data file per fit and constant, one row per v. A constant whose SPFIT
  operator ends in 9999 (4 lambda_D0, 4 lambda_H0, gamma_D1, gamma_H1) is common to
  v = 0 and 1, one row. ``\times 10^{-k}`` factors are the datasets' scales. The value's
  uncertainty in parentheses is in units of its last printed digit and is written as
  the absolute sigma. T_1 (the band origin) has no quantity in MHz and the b 1Sigma_g+
  rows are not taken (the record's notes).
- Table I (``tbl_f249a8f9f01f``, p. 7; no CSV, the JSON cells match the image): the
  X-state Dunham constants of 16O16O, the coefficients of the record's formula
  ``x_state_dunham``. ``extract`` checks that the formula's coefficient strings are the
  printed cells (thin spaces removed, ``\times 10^{-k}`` written ``e-k``, the footnote
  mark ``a`` of U_03 dropped), so the formula cannot drift from the table.
"""

import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, write

DIRECTORY = "Yu_2012_HighResolutionSpectral"
TABLE_I = "tbl_f249a8f9f01f"
TABLE_III = "tbl_844a6de53312"

# Table III's X-state rows: printed label -> (constant, row, scale exponent).
TABLE_III_ROWS = {
    "B_0": ("B", "v=0", 0),
    "B_1": ("B", "v=1", 0),
    "-D_0": ("minus_D", "v=0", 0),
    "-D_1": ("minus_D", "v=1", 0),
    "H_0": ("H", "v=0", -6),
    "H_1": ("H", "v=1", -6),
    "4\\lambda_0": ("4lambda", "v=0", 0),
    "4\\lambda_1": ("4lambda", "v=1", 0),
    "4\\lambda_{D0}": ("4lambda_D", "v=0 and 1", 0),
    "4\\lambda_{H0}": ("4lambda_H", "v=0 and 1", -6),
    "\\gamma_0": ("gamma", "v=0", 0),
    "\\gamma_1": ("gamma", "v=1", 0),
    "\\gamma_{D1}": ("gamma_D", "v=0 and 1", -3),
    "\\gamma_{H1}": ("gamma_H", "v=0 and 1", -9),
}
_CELL = re.compile(
    r"^(?P<value>-?[\d. ]+?)(?:\((?P<sigma>\d+)\))?"
    r"(?: *(?:\\times|\N{MULTIPLICATION SIGN}) *10\^\{(?P<power>-\d+)\})?"
    r"(?P<mark>a)?$"
)

# Table I's X-state rows: printed label -> the formula's coefficient key.
TABLE_I_ROWS = {
    "U_{1,0} \\times \\mu^{-0.5}/2": "Y10",
    "U_{1,0}\\Delta^{O}_{1,0}m_e \\times \\mu^{-1.5}/2": "Y10_BOB",
    "U_{2,0} \\times \\mu^{-1.0}/4": "Y20",
    "U_{3,0} \\times \\mu^{-1.5}/8": "Y30",
    "U_{4,0} \\times \\mu^{-2.0}/16": "Y40",
    "U_{5,0} \\times \\mu^{-2.5}/32": "Y50",
    "U_{6,0} \\times \\mu^{-3.0}/64": "Y60",
    "U_{0,1} \\times \\mu^{-1.0}": "Y01",
    "U_{0,1}\\Delta^{O}_{0,1}m_e \\times \\mu^{-2.0}": "Y01_BOB",
    "U_{1,1} \\times \\mu^{-1.5}/2": "Y11",
    "U_{2,1} \\times \\mu^{-2.0}/4": "Y21",
    "U_{3,1} \\times \\mu^{-2.5}/8": "Y31",
    "U_{4,1} \\times \\mu^{-3.0}/16": "Y41",
    "U_{5,1} \\times \\mu^{-3.5}/32": "Y51",
    "U_{0,2} \\times \\mu^{-2.0}": "Y02",
    "U_{1,2} \\times \\mu^{-2.5}/2": "Y12",
    "U_{2,2} \\times \\mu^{-3.0}/4": "Y22",
    "U_{3,2} \\times \\mu^{-3.5}/8": "Y32",
    "U_{0,3} \\times \\mu^{-3.0}": "Y03",
    "U_{1,3} \\times \\mu^{-3.5}/2": "Y13",
    "4 \\times U_{0,0}^{\\lambda} \\times \\mu^{-0.0}": "L00",
    "4 \\times U_{0,0}^{\\lambda}\\Delta^{O}_{0,0}^{\\lambda,O}m_e "
    "\\times \\mu^{-1.0}": "L00_BOB",
    "4 \\times U_{1,0}^{\\lambda} \\times \\mu^{-0.5}/2": "L10",
    "4 \\times U_{2,0}^{\\lambda} \\times \\mu^{-1.0}/4": "L20",
    "4 \\times U_{3,0}^{\\lambda} \\times \\mu^{-1.5}/8": "L30",
    "4 \\times U_{4,0}^{\\lambda} \\times \\mu^{-2.0}/16": "L40",
    "4 \\times U_{0,1}^{\\lambda} \\times \\mu^{-1.0}": "L01",
    "4 \\times U_{1,1}^{\\lambda} \\times \\mu^{-1.5}/2": "L11",
    "4 \\times U_{2,1}^{\\lambda} \\times \\mu^{-2.0}/4": "L21",
    "4 \\times U_{0,2}^{\\lambda} \\times \\mu^{-2.0}": "L02",
    "4 \\times U_{0,3}^{\\lambda} \\times \\mu^{-3.0}": "L03",
    "U_{0,0}^{\\gamma} \\times \\mu^{-1.0}": "G00",
    "U_{0,0}^{\\gamma}\\Delta^{O}_{0,0}^{\\gamma,O}m_e \\times \\mu^{-2.0}": "G00_BOB",
    "U_{1,0}^{\\gamma} \\times \\mu^{-1.5}/2": "G10",
    "U_{0,1}^{\\gamma} \\times \\mu^{-2.0}": "G01",
    "U_{0,2}^{\\gamma} \\times \\mu^{-3.0}": "G02",
}
X_ROWS = range(3, 39)  # Table I's X-state rows (16O16O); the 16O17O rows follow.


def _parse(cell: str) -> re.Match[str]:
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return match


def _sigma(value: str, digits: str | None) -> str:
    """Return the parenthesised last-digits uncertainty as an absolute number."""
    if digits is None:
        return ""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _table_iii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    if tuple(grid[0][2:]) != ("SPFIT operator", "Fit #1", "Fit #2"):
        msg = f"Table III header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:16]:  # the X 3Sigma_g- block (its label spans 15 rows)
        if cells[1] == "T_1":
            continue
        constant, row, power = TABLE_III_ROWS[cells[1]]
        printed = " | ".join(cells[1:])
        for fit, cell in (("fit1", cells[3]), ("fit2", cells[4])):
            match = _parse(cell)
            if int(match["power"] or 0) != power or match["mark"]:
                msg = f"Table III cell {cell!r}: expected 10^{power}"
                raise ValueError(msg)
            value = match["value"].strip()
            files.setdefault(f"O2/yu_2012__{fit}_{constant}_O2.csv", []).append(
                Point(row, row, value, _sigma(value, match["sigma"]), printed)
            )
    return files


def _coefficient(cell: str) -> tuple[str, str]:
    """Return a Table I cell as the formula's (value, sigma) strings."""
    match = _parse(cell)
    value = match["value"].replace(" ", "")
    sigma = _sigma(value, match["sigma"])
    if match["power"]:
        value = f"{value}e{match['power']}"
        sigma = f"{sigma}e{match['power']}" if sigma else ""
    return value, sigma


def _check_table_i() -> None:
    """Raise unless the record's formula coefficients are Table I's X-state cells."""
    with (DATA / "sources" / "yu_2012.toml").open("rb") as handle:
        record = tomllib.load(handle)
    (formula,) = [f for f in record["formulas"] if f["id"] == "x_state_dunham"]
    recorded = formula["coefficients"]
    grid = json_table(DIRECTORY, TABLE_I)
    expected: dict[str, str] = {}
    for index in X_ROWS:
        cells = grid[index]
        key = TABLE_I_ROWS[cells[1]]
        value, sigma = _coefficient(cells[2])
        expected[key] = value
        if sigma:
            expected[f"{key}_sigma"] = sigma
    if recorded != expected:
        diff = {
            k: (recorded.get(k), v) for k, v in expected.items() if recorded.get(k) != v
        }
        extra = set(recorded) - set(expected)
        msg = f"x_state_dunham coefficients differ from Table I: {diff}, extra {extra}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Table III as data files, after checking the formula against Table I."""
    _check_table_i()
    return _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

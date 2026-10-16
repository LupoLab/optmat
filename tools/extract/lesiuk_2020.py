# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lesiuk, Przybytek & Jeziorski (2020): neon polarizability, Cauchy moments.

- Table IX (``tbl_68954d2edbcd``, cells file only; image checked): the error budget of
  alpha_0, alpha_2, alpha_4. Only the row "Total" is extracted: the coefficients of
  alpha(omega) = alpha_0 + alpha_2 omega^2 + alpha_4 omega^4, i.e. S(-2), S(-4), S(-6)
  in atomic units, rows indexed by the moment order k. The contribution rows are
  components of the total and are not extracted.
- Table XI (``tbl_c9bc797bccfd``, cells file only; image checked): literature values.
  One file per row; "This work" repeats Table IX and is not extracted. Orcutt and Cole
  have two rows under one (row-spanning) label, the second marked ^a ("after applying
  the correction for the compressibility"); the label is carried down and the marked
  value goes to its own file. Static-only rows become polarizabilities at infinite
  wavelength, rows with dispersion coefficients S(k); the '-' cells are skipped.

Uncertainties printed in the last digits, ``2.66080(36)``, go to the sigma column in the
value's unit.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Lesiuk_2020_TheoreticalDeterminationPolarizability"
TABLE_IX = "tbl_68954d2edbcd"
TABLE_XI = "tbl_c9bc797bccfd"
PREFIX = "Ne/lesiuk_2020__"
ORDERS = {1: "-2", 2: "-4", 3: "-6"}
# Table XI label (start) -> dataset id; the label of Gaiser's two rows ends in the
# reference number.
LITERATURE = {
    "Chan et al.": "chan_1992",
    "Kumar and Meath": "kumar_2010",
    "Orcutt and Cole": "orcutt_1967",
    "Gaiser and Fallmuth [91]": "gaiser_2010",
    "Gaiser and Fallmuth [7]": "gaiser_2018",
    "Pawlowski": "pawlowski_2005",
    "Klopper": "klopper_2004",
    "Larsen": "larsen_1999",
}
_WITH_ERROR = re.compile(r"^([-\d.]+)\((\d+)\)$")
_MARK_A = "^a"


def _split(cell: str) -> tuple[str, str]:
    """Split ``2.66080(36)`` into the value and its uncertainty in the value's unit."""
    text = cell.strip()
    match = _WITH_ERROR.match(text)
    if match is None:
        return number(text), ""
    value, digits = match.groups()
    decimals = len(value.split(".")[1]) if "." in value else 0
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _moments(cells: list[str]) -> list[Point]:
    printed = " | ".join(cells)
    points = []
    for column, order in ORDERS.items():
        if cells[column].strip() in {"-", ""}:
            continue
        value, sigma = _split(cells[column])
        points.append(Point(order, order, value, sigma, printed))
    return points


def extract() -> dict[str, list[Point]]:
    """Return Table IX's totals and Table XI's literature values."""
    files: dict[str, list[Point]] = {}
    (total,) = (r for r in json_table(DIRECTORY, TABLE_IX) if r[0].strip() == "Total")
    files[f"{PREFIX}cauchy_Ne.csv"] = _moments(total)
    label = ""
    for cells in json_table(DIRECTORY, TABLE_XI):
        label = cells[0].strip() or label  # Orcutt and Cole's label spans two rows
        name = next((v for k, v in LITERATURE.items() if label.startswith(k)), None)
        if name is None or not cells[1].strip():
            continue  # the header, the group rows and "This work"
        printed = " | ".join(cells)
        if cells[1].strip().endswith(_MARK_A):
            value, sigma = _split(cells[1].strip().removesuffix(_MARK_A))
            points = [Point("alpha_0", "inf", value, sigma, printed)]
            name += "_compressibility"
        elif all(c.strip() == "-" for c in cells[2:]):
            value, sigma = _split(cells[1])
            points = [Point("alpha_0", "inf", value, sigma, printed)]
        else:
            points = _moments(cells)
        files[f"{PREFIX}{name}_Ne.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

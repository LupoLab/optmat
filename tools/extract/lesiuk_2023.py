# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lesiuk & Jeziorski (2023): argon static polarizability and Cauchy moments.

- Table XI (``tbl_cb5cfdf665dc``, cells file only; image checked): the error budget of
  alpha_0, alpha_2 and alpha_4 (columns n = 0, 2, 4). Only the row "Total" is
  extracted: the coefficients of alpha(omega) = alpha_0 + alpha_2 omega^2 + alpha_4
  omega^4 + ..., i.e. the dipole sums S(-2), S(-4), S(-6) in atomic units; the rows are
  indexed by the moment order k = -(n + 2). The contribution rows (Hartree-Fock,
  coupled-cluster increments, core correlation, relativistic, QED, finite nuclear size
  and mass) are components of the total and are not extracted.
- Table XII (``tbl_787fa31e3c55``, cells file only; image checked): literature values of
  alpha_0, alpha_2, alpha_4. One file per row; "This work" repeats Table XI and is not
  extracted. The static-only rows (dielectric measurements, one finite-field
  calculation) become polarizabilities at infinite wavelength, the rows with
  dispersion coefficients become S(k) like Table XI.

Every cell prints the uncertainty in the last digits in parentheses, ``11.0775(19)``;
it is written to the sigma column in the value's unit (0.0019).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Lesiuk_2023_FirstPrinciplesCalculation"
TABLE_XI = "tbl_cb5cfdf665dc"
TABLE_XII = "tbl_787fa31e3c55"
PREFIX = "Ar/lesiuk_2023__"
# Column n of alpha_n -> moment order k of S(k).
ORDERS = {1: "-2", 2: "-4", 3: "-6"}
# Table XII row label (start) -> dataset id; the static-only rows are polarizabilities.
LITERATURE = {
    "Kumar and Thakkar": "kumar_2010",
    "Orcutt and Cole": "orcutt_1967",
    "Buckley": "buckley_2000",
    "Gaiser and Fellmuth": "gaiser_2018",
    "Pawłowski": "pawlowski_2005",
    "Lupinetti": "lupinetti_2005",
}
_WITH_ERROR = re.compile(r"^([-\d.]+)\((\d+)\)$")
_DASH = "\N{EM DASH}"


def _split(cell: str) -> tuple[str, str]:
    """Split ``11.0775(19)`` into the value and its uncertainty in the value's unit."""
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
        if cells[column].strip() in {_DASH, ""}:
            continue
        value, sigma = _split(cells[column])
        points.append(Point(order, order, value, sigma, printed))
    return points


def extract() -> dict[str, list[Point]]:
    """Return Table XI's totals and Table XII's literature values."""
    files: dict[str, list[Point]] = {}
    (total,) = (r for r in json_table(DIRECTORY, TABLE_XI) if r[0].strip() == "Total")
    files[f"{PREFIX}cauchy_Ar.csv"] = _moments(total)
    for cells in json_table(DIRECTORY, TABLE_XII):
        name = next((v for k, v in LITERATURE.items() if cells[0].startswith(k)), None)
        if name is None:
            continue  # the header, the two group rows and "This work"
        if all(c.strip() == _DASH for c in cells[2:]):
            value, sigma = _split(cells[1])
            points = [Point("alpha_0", "inf", value, sigma, " | ".join(cells))]
        else:
            points = _moments(cells)
        files[f"{PREFIX}{name}_Ar.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

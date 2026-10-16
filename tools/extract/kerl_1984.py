# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kerl and Häusler (1984): molar refractivities of Ar and CH4 at 546.226 nm.

Only Ar and CH4 are curated species; C2H6, C3H8 and C(CH3)4 are not extracted.

Table 1 (``tbl_6ccbf5b10dc1``, p. 994, uncaptioned in the library and without a CSV
rendering: read from its JSON cells, checked against the image, whose caption reads
"[n(lambda_0) - 1] = a_1 (p/RT) + a_2 (p/RT)^2 for lambda_0 = 546.226 nm") gives per
isotherm T, p_max, K, |s|, a_1, a_2, A_R and -B. The A_R column (cm^3/mol) is taken,
one file per temperature, with its standard deviation of the last digits in
parentheses. Values in brackets are "uncertain ... by no means [used] in all further
calculations" (p. 994): the brackets are dropped here and the dataset is excluded in
the record. A star before T marks the improved sample cell; it is dropped.

Table 2 (``tbl_44fdf7d2e6c3``, p. 994) gives the arithmetic mean A_R-bar of Table 1's
A_R values (the bracketed ones left out), with the standard deviation of the values.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Kerl_1984_MeanPolarizabilitiesSecond"
TABLE_1 = "tbl_6ccbf5b10dc1"
TABLE_2 = "tbl_44fdf7d2e6c3"
WAVELENGTH = "546.226"  # nm, vacuum (Table 1 caption and the text, pp. 993 and 995)
GASES = {"Ar": "Ar", "\\mathrm{CH_4}": "CH4", "CH_4": "CH4"}
_CELL = re.compile(r"^\[?(?P<value>[\d.]+)\((?P<sigma>\d+)\)\]?$")
_TEMPERATURE = re.compile(r"^\*?\s*(?P<value>[\d.]+)\(\d+\)$")


def _value(cell: str) -> tuple[str, str]:
    """Return a cell's value and its last-digits standard deviation as a number."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value = number(match["value"])
    decimals = len(value.partition(".")[2])
    return value, f"{int(match['sigma']) * 10.0**-decimals:.{decimals}f}"


def _table_1() -> dict[str, list[Point]]:
    """Return one A_R file per gas and temperature."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_1)[1:]:
        if cells[0]:
            gas = GASES.get(cells[0], "")
        if not gas:
            continue  # C2H6, C3H8, C(CH3)4 are not curated
        match = _TEMPERATURE.match(cells[1].strip())
        if match is None:
            msg = f"unexpected temperature {cells[1]!r}"
            raise ValueError(msg)
        temperature = match["value"].replace(".", "_")
        value, sigma = _value(cells[7])
        point = Point(WAVELENGTH, WAVELENGTH, value, sigma, " | ".join(cells))
        files[f"{gas}/kerl_1984__A_R_{temperature}K_{gas}.csv"] = [point]
    return files


def _table_2() -> dict[str, list[Point]]:
    """Return the mean A_R of Ar and CH4."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        gas = GASES.get(cells[0], "")
        if not gas:
            continue
        value, sigma = _value(cells[5])
        point = Point(WAVELENGTH, WAVELENGTH, value, sigma, " | ".join(cells))
        files[f"{gas}/kerl_1984__A_R_mean_{gas}.csv"] = [point]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 1 and Table 2 molar refractivities of Ar and CH4."""
    return _table_1() | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

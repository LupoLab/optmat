# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gaiser & Fellmuth (2018): DCGT molar polarizabilities of He, Ne and Ar.

- Table III (``tbl_a15621c7c166``, cells file only; image checked): A_eps in cm^3/mol,
  this work and the literature values of Fig. 2. Digits are printed in groups of three
  (``0.517 254 4(10)``): the spaces are removed; the parenthesised last-digit
  uncertainty (a standard uncertainty, caption) goes to the sigma column in the value's
  unit. The "Experiment" cells hold several values, each followed by its reference
  number (``0.517 40(70) [16]0.517 47(20) [17]...``); they are split on the reference
  numbers, one file per reference and gas.
- Table II (``tbl_ec04560e1d2e``, cells file only; image checked): the per-capacitor
  results A_eps/R x 10^8 in m^3 K/J (row "A_eps / R x 10^8", the first of the two rows
  of that name; the second is the weighted mean, which Table III restates as A_eps).
  One file per gas, rows by capacitor. Their conversion to A_eps (scale in the record)
  multiplies by the molar gas constant.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Gaiser_2018_PolarizabilityHeliumNeon"
TABLE_II = "tbl_ec04560e1d2e"
TABLE_III = "tbl_a15621c7c166"
GASES = {1: "He", 2: "Ne", 3: "Ar"}
# Table III reference number -> file stem.
REFERENCES = {
    "16": "orcutt_1967",
    "17": "buckley_2000",
    "3": "schmidt_2007",
    "10": "gaiser_2010",
    "1,2": "piszczatowski_2015",
    "18,19": "klopper_2004",
    "15": "lupinetti_2005",
}
_ENTRY = re.compile(
    r"(?P<number>[\d. ]+)\((?P<digits>\d+)\)\s*(?:\[(?P<ref>[\d,]+)\])?"
)


def _split(number: str, digits: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    value = number.replace(" ", "")
    decimals = len(value.split(".")[1])
    return value, f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_III)[1:]:
        label = cells[0].strip()
        printed = " | ".join(cells)
        for column, gas in GASES.items():
            entries = list(_ENTRY.finditer(cells[column]))
            if not entries or "".join(m[0] for m in entries) != cells[column]:
                msg = f"unreadable cell {cells[column]!r}"
                raise ValueError(msg)
            for entry in entries:
                value, sigma = _split(entry["number"], entry["digits"])
                if label == "This work":
                    name = f"molar_polarizability_{gas}"
                else:
                    name = f"{REFERENCES[entry['ref']]}_{gas}"
                path = f"{gas}/gaiser_2018__{name}.csv"
                files[path] = [Point("static", "inf", value, sigma, printed)]
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    gas_of_column: dict[int, str] = {}
    gas = ""
    for column, cell in enumerate(grid[0][1:], start=1):
        gas = cell.strip() or gas  # the gas header spans its capacitors
        gas_of_column[column] = gas
    per_capacitor = next(r for r in grid if r[0].startswith("A_\\epsilon / R"))
    files: dict[str, list[Point]] = {}
    for column in range(1, len(per_capacitor)):
        gas, capacitor = gas_of_column[column], grid[1][column].strip()
        printed = f"{gas} | {capacitor} | {per_capacitor[0]} | {per_capacitor[column]}"
        path = f"{gas}/gaiser_2018__per_capacitor_{gas}.csv"
        files.setdefault(path, []).append(
            Point(capacitor, "inf", per_capacitor[column].strip(), "", printed)
        )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table III's molar polarizabilities and Table II's per-capacitor values."""
    return _table_iii() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Nibbering et al. (1997): n2 of gases, Tables 1-3.

- Table 1 (p. 651, CSV): literature n2 (1e-23 m2/W) at 1.0 bar; only the curated gases,
  all at 248 nm with 10-ps pulses (Shaw et al., Ref. 9; Xe: Lehmberg et al., Ref. 8).
  Structural fixes: the pdf's symbol font maps the plus-minus sign to "6" and the minus
  sign to "2" (the table image shows "7.5 +- 2.5", "-0.02 +- 0.01", "-963 +- 112"; the
  pdf text layer has the same "6" and "2"); the library cells inherit them. The script
  splits on " 6 " and restores the two minus signs, checking the raw cells.
- Table 2 (p. 653, CSV): n2,total and n2,electronic (1e-23 m2/W) at 1 atm, 293 K,
  800 nm, "Taken from Ref. 1" (Hellwarth et al. 1990). Air skipped.
- Table 3 (p. 658, JSON; no CSV): this paper's n2,instant (1e-23 m2/W), n2,rot-Raman
  (not extracted, m2/J) and beta = Delta alpha derived and from the literature (1e-40
  C2 m2/J) at 1.0 bar, 300 K, 120 fs, 800 nm. Cells print ``value \pm error`` or "-".
  Air skipped.

The CSV and JSON cells match the table images (checked).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Nibbering_1997_DeterminationInertialContribution"
LITERATURE = "tbl_5f0a7a95fe44"
HELLWARTH = "tbl_aa252cc4e383"
THIS_WORK = "tbl_e8b12922f707"
CURATED = {"N2", "O2", "Ar", "CH4", "Ne", "He", "Xe", "H2", "SF6"}
TEX = {r"\text{SF}_{6}": "SF6", r"\text{N}_{2}": "N2", r"\text{O}_{2}": "O2"}
# The lost minus signs (the raw cell, what the image shows).
MINUS = {"He": ("20.02 6 0.01", "-0.02"), "Xe": ("2963 6 112", "-963")}


def _literature() -> dict[str, list[Point]]:
    """Table 1: the curated gases' literature values at 248 nm."""
    files = {}
    for cells in table(DIRECTORY, LITERATURE)[2:]:
        gas = cells[0]
        if gas not in CURATED:
            continue  # air
        value, error = cells[4].split(" 6 ")
        if gas in MINUS:
            raw, value = MINUS[gas]
            if cells[4] != raw:
                raise ValueError(cells[4])
        files[f"{gas}/nibbering_1997__literature_{gas}.csv"] = [
            Point(cells[3], cells[3], number(value), number(error), " | ".join(cells))
        ]
    return files


def _hellwarth() -> dict[str, list[Point]]:
    """Table 2: Hellwarth et al.'s total and electronic n2 converted to m2/W."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, HELLWARTH)[1:]:
        gas = cells[0]
        if gas not in CURATED:
            continue  # air
        for column, name in ((3, "total"), (4, "electronic")):
            if cells[column]:
                path = f"{gas}/nibbering_1997__hellwarth_{name}_{gas}.csv"
                # Footnote c: the electronic values are "Estimated for the
                # long-wavelength limit", so the table's 800 nm does not apply to them.
                row = "LWL" if name == "electronic" else cells[2]
                wavelength = "inf" if name == "electronic" else cells[2]
                files[path] = [
                    Point(row, wavelength, number(cells[column]), "", " | ".join(cells))
                ]
    return files


def _this_work() -> dict[str, list[Point]]:
    """Table 3: n2,instant and the anisotropy beta, derived and literature."""
    files = {}
    for cells in json_table(DIRECTORY, THIS_WORK)[2:]:
        gas = TEX.get(cells[0], cells[0])
        if gas not in CURATED:
            continue  # air
        printed = " | ".join(cells)
        value, error = cells[1].split(r"\pm")
        files[f"{gas}/nibbering_1997__n2_instant_{gas}.csv"] = [
            Point("800", "800", number(value), number(error), printed)
        ]
        if cells[3] != "-":
            files[f"{gas}/nibbering_1997__beta_derived_{gas}.csv"] = [
                Point("800", "800", number(cells[3]), "", printed)
            ]
            files[f"{gas}/nibbering_1997__beta_lin_{gas}.csv"] = [
                Point("inf", "not stated", number(cells[4]), "", printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return every data file of the source."""
    return {**_literature(), **_hellwarth(), **_this_work()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

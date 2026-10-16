# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Watson and Ramaswamy (1936): refractivities at 25 C and molar polarizations.

The paper prints the decimal point as a raised dot (``402·94``), read here as a point.
Every refractivity of the authors' own is ``(mu - 1)`` at 25 C and 760 mm reduced to the
ideal gas; wavelengths are vacuum wavelengths (p. 3). Images of every table checked.

- Table II (``tbl_56f6683d6c06`` with its continuation ``tbl_5065fe5949d9``; the
  continuation's CSV takes its first data row, C2H2, as the header): ``(mu - 1) x 10^6``
  at five lines for CH4, SF6 and CO2 (the other 20 gases are not curated). Under the
  red, Cd green and Cd blue values the image prints on a second line the difference
  calculated - observed x 10^8; the cells run the two together (``402·94-7``,
  ``411·130`` for 411.13 with difference 0). The value is the number with two
  decimals; the difference is not extracted (the formula is recorded instead).
- Table III (``tbl_97d70f8c24f3``): at 5462.25 A, the authors' own values (column 1,
  CO, N2O and N2; O2 has none) and Cuthbertson's (column 2, recalculated to 25 C and
  the ideal gas, marked with the reference's footnote sign, stripped here). NO, NH3
  are not curated.
- Table V (``tbl_23af967b0794``), ``(mu - 1) x 10^6`` at 0 C and 5462 A by several
  authors: CH4 (Tausz and Hornung, Friberg, the authors) and SF6 (Klemm and Henkel, the
  authors). The author is named per row; the gas on the first row of its block.
- Table VI (``tbl_06bdd0f09e6a``), the authors' values reduced to 0 C (x 1.09157):
  CO2, CO, N2O, N2 (right-hand column pair).
- Table VII (``tbl_d67d2c70949c``), ``P_A + P_E`` (the molar polarization from the
  dielectric coefficient, cm^3) of the curated gases; CH4 is printed twice with the
  same values (both blocks), taken once. ``A`` is argon.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Watson_1936_RefractiveIndexDispersion"
TABLE_II = ("tbl_56f6683d6c06", "tbl_5065fe5949d9")
TABLE_III = "tbl_97d70f8c24f3"
TABLE_V = "tbl_23af967b0794"
TABLE_VI = "tbl_06bdd0f09e6a"
TABLE_VII = "tbl_d67d2c70949c"
WAVELENGTHS = ("6440.24", "5462.25", "5087.23", "4801.25", "4359.56")
# Printed labels (with or without braces) -> curated gas.
GASES = {
    "CH_4": "CH4",
    "SF_6": "SF6",
    "CO_2": "CO2",
    "CO": "CO",
    "N_2O": "N2O",
    "N_2": "N2",
    "O_2": "O2",
    "H_2": "H2",
    "He": "He",
    "Ne": "Ne",
    "A": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
}
# Table V: author abbreviation -> dataset suffix.
AUTHORS = {
    "T. and H.": "tausz_hornung",
    "F.": "friberg",
    "W. and R.": "zero_celsius",
    "K. and H.": "klemm_henkel",
}
_VALUE_AND_DIFFERENCE = re.compile(r"^(\d+·\d\d)([+-]?\d+)?$")


def _gas(label: str) -> str | None:
    """Return the curated gas a printed label names, or None."""
    return GASES.get(label.replace("{", "").replace("}", "").strip())


def _decimal(cell: str) -> str:
    """Read the printed raised decimal point as a decimal point."""
    return number(cell.replace("·", ".").rstrip("*†‡"))


def _table_ii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_II[0])[2:] + table(DIRECTORY, TABLE_II[1])
    for cells in rows:
        gas = _gas(cells[0])
        if gas is None:
            continue
        points = []
        for wavelength, cell in zip(WAVELENGTHS, cells[1:6], strict=True):
            match = _VALUE_AND_DIFFERENCE.match(cell.strip())
            if match is None:
                msg = f"unexpected Table II cell {cell!r}"
                raise ValueError(msg)
            value = _decimal(match.group(1))
            points.append(Point(wavelength, wavelength, value, "", " | ".join(cells)))
        files[f"{gas}/watson_1936__dispersion_{gas}.csv"] = points
    return files


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_III)[1:]:
        gas = _gas(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        prefix = f"{gas}/watson_1936__"
        if cells[1]:
            files[f"{prefix}green_{gas}.csv"] = [
                Point("5462.25", "5462.25", _decimal(cells[1]), "", printed)
            ]
        files[f"{prefix}cuthbertson_{gas}.csv"] = [
            Point("5462.25", "5462.25", _decimal(cells[2]), "", printed)
        ]
    return files


def _table_v() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    # The two column blocks (gas, value, author) each name a gas on the first row of
    # its group only; CH4 is in the left-hand block, SF6 in the right-hand one.
    blocks = ((0, 1, 3), (4, 5, 7))
    current: dict[int, str | None] = {0: None, 4: None}
    for cells in json_table(DIRECTORY, TABLE_V)[1:]:
        for label, value, author in blocks:
            if cells[label]:
                current[label] = _gas(cells[label])
            gas = current[label]
            if gas not in {"CH4", "SF6"} or not cells[value]:
                continue
            suffix = AUTHORS[cells[author].strip()]
            files[f"{gas}/watson_1936__{suffix}_{gas}.csv"] = [
                Point("5462", "5462", _decimal(cells[value]), "", " | ".join(cells))
            ]
    return files


def _table_vi() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_VI)[1:]:
        gas = _gas(cells[2])
        if gas is None:
            continue
        files[f"{gas}/watson_1936__zero_celsius_{gas}.csv"] = [
            Point("5462", "5462", _decimal(cells[3]), "", " | ".join(cells))
        ]
    return files


def _table_vii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_VII)[2:]:
        for label, total in ((0, 1), (5, 6)):
            gas = _gas(cells[label])
            if gas is None:
                continue
            path = f"{gas}/watson_1936__dielectric_{gas}.csv"
            if path in files:  # CH4 is printed twice with the same values
                continue
            files[path] = [
                Point("inf", "inf", _decimal(cells[total]), "", " | ".join(cells))
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the curated gases' values of Tables II, III, V, VI and VII."""
    return _table_ii() | _table_iii() | _table_v() | _table_vi() | _table_vii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

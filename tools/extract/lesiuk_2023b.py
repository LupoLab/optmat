# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lesiuk & Jeziorski (arXiv 2023): the magnetic susceptibility of Ne and Ar.

Table IX (``tbl_070dd39bd6ec``, cells file only; image checked) is the budget of the
static diamagnetic susceptibility chi_0 of neon and argon, every entry multiplied by
10^5 (caption); only its last row, the total chi_0, becomes a data file (the
contributions are components of it). The parenthesised last-digit uncertainty goes to
the sigma column.

Table X (``tbl_e0613af46035``, CSV; image checked) compares chi_0 with the literature in
atomic units, each cell with its own power of ten (``-8.574(9) · 10-5``): the mantissa
is the value, and the power is the dataset's scale in the source record, which the
script checks. Only the experimental rows (Havens, Mann, Barter et al.) and the
authors' own earlier values (Lesiuk and Jeziorski [25, 27]) are extracted; the older
theoretical rows, printed without uncertainties, and the "this work" row (Table IX's
totals again) are not.

Barter et al.'s argon cell is garbled in the extraction: the image shows two lines,
``-2.16(2)^a · 10^-4`` (note a: the original error estimate) over
``-2.16(15)^b · 10^-5`` (note b: the revised estimate of Rourke 2021), and the library
interleaves their characters. The script checks the garbled string and takes the first
line, the value as Barter et al. printed it; the second is in the record's notes.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Lesiuk_2023_DiamagneticSusceptibilityNeon"
TABLE_IX = "tbl_070dd39bd6ec"
TABLE_X = "tbl_e0613af46035"
COLUMNS = {"Ne": 1, "Ar": 2}
# Table X row (reference as printed) -> dataset name, and the power of ten each gas's
# cell prints (the dataset's scale in the record).
REFERENCES = {
    "Havens [58]": ("havens_1933", {"Ne": -5, "Ar": -4}),
    "Mann [59]": ("mann_1936", {"Ne": -5, "Ar": -4}),
    "Barter el al. [60]": ("barter_1960", {"Ne": -4, "Ar": -4}),
    "Lesiuk and Jeziorski [25, 27]": ("lesiuk_jeziorski", {"Ne": -5, "Ar": -4}),
}
# The two interleaved lines of Barter et al.'s argon cell (see the module docstring).
BARTER_AR_GARBLED = "--22..1166((125)a)b··1100--45"
BARTER_AR_FIRST_LINE = "-2.16(2) · 10-4"
_VALUE = re.compile(r"^(?P<number>-?\d+\.\d+)(?:\((?P<digits>\d+)\))?$")
_TIMES_POWER = re.compile(r"^(?P<mantissa>\S+) · 10(?P<power>-\d+)$")


def _split(cell: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty (or "") in the value's unit."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["number"], match["digits"]
    if digits is None:
        return value, ""
    decimals = len(value.split(".")[1])
    return value, f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _totals() -> dict[str, list[Point]]:
    """Return Table IX's total chi_0 of each gas (x 10^5)."""
    rows = json_table(DIRECTORY, TABLE_IX)
    assert rows[0] == ["contribution", "Ne", "Ar"], rows[0]
    assert rows[-2][0] == "total" and rows[-1][0] == "\\chi_0", rows[-2:]
    files = {}
    for gas, column in COLUMNS.items():
        value, sigma = _split(rows[-1][column])
        files[f"{gas}/lesiuk_2023b__chi0_{gas}.csv"] = [
            Point("chi_0", "", value, sigma, " | ".join(rows[-1]))
        ]
    return files


def _literature() -> dict[str, list[Point]]:
    """Return Table X's experimental rows and the authors' earlier values."""
    rows = table(DIRECTORY, TABLE_X)
    files = {}
    for cells in rows:
        reference = cells[0].strip()
        if reference not in REFERENCES:
            continue
        name, powers = REFERENCES[reference]
        for gas, column in COLUMNS.items():
            cell = cells[column].strip()
            if reference == "Barter el al. [60]" and gas == "Ar":
                assert cell == BARTER_AR_GARBLED, cell
                cell = BARTER_AR_FIRST_LINE
            match = _TIMES_POWER.match(cell)
            if match is None or int(match["power"]) != powers[gas]:
                msg = f"{reference} {gas}: unexpected cell {cell!r}"
                raise ValueError(msg)
            value, sigma = _split(match["mantissa"])
            files[f"{gas}/lesiuk_2023b__{name}_{gas}.csv"] = [
                Point("chi_0", "", value, sigma, " | ".join(cells))
            ]
    if len(files) != 2 * len(REFERENCES):
        msg = f"found {len(files)} of {2 * len(REFERENCES)} Table X values"
        raise ValueError(msg)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the totals of Table IX and the extracted rows of Table X."""
    return _totals() | _literature()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

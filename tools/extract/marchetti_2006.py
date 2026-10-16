# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Marchetti & Simili (2006): CO2, N2, O2, He refractivity at 10.57 um.

Table 1 (``tbl_fcf07d146f15``, OCR re-extraction, all numbers agreeing with the text
layer; image checked) gives the pressure coefficient ``alpha_P`` (n = 1 + alpha_P P)
in ``(atm^-1) x 10^-4`` at 23 degrees C (captioned "@ 296 °C": 296 K), with the
uncertainty in the last digits in parentheses on the next row. Each gas's value
becomes one point, with the uncertainty converted to the value's unit (``(12)`` on
2.725 is 0.012).

Table 2 (``tbl_a64b03e950d7``) gives ``alpha_rho`` = (n - 1)/rho in cm^3/mol: this work
(Table 1 times R T), Simmons (1978) at 10.59 and 0.633 um, Achtermann et al. (1991) at
0.633 um, and theory (Ciddor, Edlen, Old et al.). The row labels are carried by the
row above when a note letter spans two rows (Simmons, Edlen, Old): the script keeps
the last note letter. The image prints Simmons's N2 value at 0.633 um as "7.219" with
"73)" wrapped onto the next line, i.e. 7.219(73), which the library cell reads as
``7.21973)``; the script splits it.

Air is not a curated species and is skipped.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Marchetti_2006_AccurateMeasurementRefractive"
TABLE_1 = "tbl_fcf07d146f15"
TABLE_2 = "tbl_a64b03e950d7"
GASES = {
    "CO_{2}": "CO2",
    "CO2": "CO2",
    "N_{2}": "N2",
    "N2": "N2",
    "O_{2}": "O2",
    "O2": "O2",
    "He": "He",
}
# Table 2 note letter -> dataset stem.
NOTES = {
    "a": "alpha_rho",
    "b": "simmons_1978",
    "c": "achtermann_1991",
    "e": "edlen_theory",
    "f": "old_theory",
}
# The cell the image shows as 7.219 over a wrapped "73)".
WRAPPED = {"7.21973)": ("7.219", "73")}
_WITH_ERROR = re.compile(r"^([\d.]+)\((\d+)\)$")


def _sigma(value: str, digits: str) -> str:
    """Give a last-digits uncertainty in the value's unit: 2.725 (12) -> 0.012."""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _table_1() -> dict[str, list[Point]]:
    header, values, errors = table(DIRECTORY, TABLE_1)
    files = {}
    for column, name in enumerate(header):
        gas = GASES.get(name.strip())
        if gas is None:
            continue  # the label column and AIR
        value = number(values[column])
        digits = errors[column].strip().strip("()")
        printed = f"{values[column]} | {errors[column]}"
        files[f"{gas}/marchetti_2006__alpha_p_{gas}.csv"] = [
            Point("10.57", "10.57", value, _sigma(value, digits), printed)
        ]
    return files


def _table_2() -> dict[str, list[Point]]:
    header, *rows = table(DIRECTORY, TABLE_2)
    files: dict[str, list[Point]] = {}
    note = ""
    for cells in rows:
        note = cells[-1].strip() or note  # a note letter spans its second row
        if note not in NOTES:
            continue  # d: Ciddor's theory for air only
        found = re.search(r"@([\d.]+)", cells[0])
        if found is None:
            msg = f"no wavelength in row label {cells[0]!r}"
            raise ValueError(msg)
        wavelength = found.group(1)
        for column, name in enumerate(header[1:-1], start=1):
            gas, cell = GASES.get(name.strip()), cells[column].strip()
            if gas is None or not cell:
                continue
            if cell in WRAPPED:
                value, digits = WRAPPED[cell]
            elif match := _WITH_ERROR.match(cell):
                value, digits = match.groups()
            else:
                value, digits = cell, ""
            sigma = _sigma(value, digits) if digits else ""
            path = f"{gas}/marchetti_2006__{NOTES[note]}_{gas}.csv"
            files.setdefault(path, []).append(
                Point(wavelength, wavelength, number(value), sigma, " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1 and 2 as data files (air skipped)."""
    return _table_1() | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Marchetti & Simili (2007): CO2 refractivity at 27 laser lines, 9.2-10.7 um.

The single table (``tbl_33adeac1f5b0``; image checked) is printed as three pairs of
columns, "CO2 line (frequency)" and "dn", the frequency a vacuum wavenumber in cm^-1
and dn the pressure coefficient (n = 1 + dn P) in 1e-4 per atm at 21 degrees C, with
"the dn errors ... in parenthesis in the unit of the last digit". Each line becomes a
row with the line label as its id (the 10-um-band lines print as ``1R..``/``1P..``);
the error is converted to the value's unit (``4.242 (18)`` -> 0.018).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Marchetti_2007_AccurateMeasurementCO2"
TABLE = "tbl_33adeac1f5b0"
_LINE = re.compile(r"^(\w+)\s*\(([\d.]+)\)$")
_VALUE = re.compile(r"^([\d.]+)\s*\((\d+)\)$")


def _sigma(value: str, digits: str) -> str:
    """Give a last-digits uncertainty in the value's unit: 4.242 (18) -> 0.018."""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the 27 lines, read down each pair of columns in turn."""
    rows = table(DIRECTORY, TABLE)[1:]
    points = []
    for pair in range(3):
        for cells in rows:
            line_cell, value_cell = cells[2 * pair], cells[2 * pair + 1]
            line, value = (
                _LINE.match(line_cell.strip()),
                _VALUE.match(value_cell.strip()),
            )
            if line is None or value is None:
                msg = f"cannot read {line_cell!r} | {value_cell!r}"
                raise ValueError(msg)
            label, wavenumber = line.groups()
            dn, digits = value.groups()
            points.append(
                Point(
                    label,
                    number(wavenumber),
                    number(dn),
                    _sigma(dn, digits),
                    f"{line_cell} | {value_cell}",
                )
            )
    return {"CO2/marchetti_2007__experiment_CO2.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

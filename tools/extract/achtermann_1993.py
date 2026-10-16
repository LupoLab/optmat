# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Achtermann, Hong, Magnus, Aziz & Slaman (1993): atoms' refractivity virials.

Table II (``tbl_d3914692b976``, cells file only; image checked): per gas and reference,
the temperature, A_R (cm^3/mol), B_R (cm^6/mol^2) and C_R (cm^9/mol^3), and other
authors' A_eps and B_eps. Only this work's measured rows are extracted: A_R to
``A_R_this_work_<gas>`` (rows are the temperatures), B_R and C_R to
``<coefficient>_this_work_<T>_<gas>``, all at the vacuum wavelength 632.99 nm (the
caption). The literature and calculated rows are other sources' values, kept in the
record's notes.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Achtermann_1993_ExperimentalDeterminationRefractivity"
TABLE_II = "tbl_d3914692b976"
WAVELENGTH = "632.99"


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``1.75±0.05``."""
    value, _, sigma = cell.partition("\N{PLUS-MINUS SIGN}")
    return number(value.strip()), sigma.strip()


def extract() -> dict[str, list[Point]]:
    """Return this work's A_R, B_R and C_R of Table II."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_II)[1:]:
        gas = cells[0] or gas
        if cells[7] != "This work":
            continue
        temperature, printed = cells[1], " | ".join(cells)
        if cells[2]:
            value, sigma = split(cells[2])
            files.setdefault(
                f"{gas}/achtermann_1993__A_R_this_work_{gas}.csv", []
            ).append(Point(temperature, WAVELENGTH, value, sigma, printed))
        for column, coefficient in ((3, "B_R"), (4, "C_R")):
            if cells[column]:
                value, sigma = split(cells[column])
                name = f"{coefficient}_this_work_{temperature}_{gas}"
                files[f"{gas}/achtermann_1993__{name}.csv"] = [
                    Point(temperature, WAVELENGTH, value, sigma, printed)
                ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hohm (1994): mean polarizabilities of H2, O2, N2O, CO2 at 325-633 nm.

Table 2 (uncaptioned in the library, no CSV rendering) is read from its JSON cells and
was checked against the image. For each gas a first row holds the input data (mean
polarizabilities alpha_1..alpha_5 in atomic units, then depolarization ratios) and a
second row, with an empty gas cell, the values of the fitted two-term formula. H2
appears twice (fit sets "a" and "b") with the same input row, taken once. A letter
after the first input cell (``5.53900b)``) is the table note giving its reference and is
stripped; "-" is no value. The depolarization ratios (anisotropies) are not extracted
here. NH3, C2H6 and cyclo-C3H6 are not curated species.

The wavelengths of the five polarizability columns are those of the table note a)
(632.990, 611.971, 594.096, 543.516 and 325.130 nm, with 0.0719812 ... 0.1401389 au).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hohm_1994_DispersionPolarizabilityAnisotropy"
TABLE_2 = "tbl_14b1cb679065"
WAVELENGTHS = ("632.990", "611.971", "594.096", "543.516", "325.130")  # note a), nm
GASES = {
    "H_2^{(a)}": "H2",
    "H_2^{(b)}": "H2",
    "O_2": "O2",
    "N_2O": "N2O",
    "CO_2": "CO2",
}
FIT_SET = {"H_2^{(a)}": "_a", "H_2^{(b)}": "_b"}
_NOTE = re.compile(r"[a-i]\)$")


def _points(cells: list[str]) -> list[Point]:
    """Return the five polarizability cells of a row, skipping "-"."""
    printed = " | ".join(cells)
    points = []
    for wavelength, cell in zip(WAVELENGTHS, cells[1:6], strict=True):
        text = _NOTE.sub("", cell.strip())
        if text != "-":
            points.append(Point(wavelength, wavelength, number(text), "", printed))
    return points


def extract() -> dict[str, list[Point]]:
    """Return the input and fitted polarizabilities per gas."""
    rows = json_table(DIRECTORY, TABLE_2)[1:]
    files: dict[str, list[Point]] = {}
    for first, second in zip(rows[::2], rows[1::2], strict=True):
        if first[0] not in GASES:
            continue
        gas = GASES[first[0]]
        suffix = FIT_SET.get(first[0], "")
        if suffix != "_b":  # set b repeats set a's input row
            files[f"{gas}/hohm_1994__input_{gas}.csv"] = _points(first)
        files[f"{gas}/hohm_1994__fit_{gas}{suffix}.csv"] = _points(second)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

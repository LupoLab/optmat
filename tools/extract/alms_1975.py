# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Alms, Burnham and Flygare (1975): mean polarizabilities of N2, N2O and CO2.

Table I (``tbl_24841f218204``, cells file only; image checked) gives, per molecule and
laser line, the depolarization ratio, the mean polarizability ``alpha`` (A^3), the
anisotropy and the two tensor components. Only the mean polarizability of the curated
gases (N2, N2O, CO2) is extracted here; the anisotropies are left for later, and the
other nine molecules (OCS, CS2, CH3Cl, CHCl3, HCCH, propyne, 2-butyne, C6H6, CH3CN) are
not curated species. The mean polarizabilities are not the authors' measurements but
Landolt-Bornstein's refractive indices converted with n - 1 = 2 pi N alpha (p. 3).

What the image shows and the cells lose: the last digit of each nitrogen polarizability
is printed as a subscript (``1.78_7``), a guard digit; the cells read it as an ordinary
digit (``1.787``), which is kept. Nitrogen's third wavelength is printed ``5148`` in
the image too (every other molecule has ``5145``): kept as printed, corrected by an
erratum in the source record. The molecule name is printed on the first row of its
block only (a row span in the cells).

The ``inf`` row of each molecule (the zero-frequency extrapolation) is written with
the wavelength ``inf``; the source record excludes it from the dataset.

Table V (``tbl_e60dab9bc511``, image checked) prints Stuart's (1952) static electronic
polarizabilities with this work's in parentheses, ``1.76(1.74)``: Stuart's mean
``alpha`` (the number before the parenthesis) is extracted for N2, N2O and CO2 as a
reproduced value at wavelength ``inf``; the parenthesised value repeats Table I.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Alms_1975_MeasurementDispersionPolarizability"
TABLE_I = "tbl_24841f218204"
TABLE_V = "tbl_e60dab9bc511"
# Printed molecule label -> curated gas.
GASES = {"N_2": "N2", "N_2O": "N2O", "CO_2": "CO2"}
WAVELENGTH, ALPHA = 1, 3


def _wavelength(cell: str) -> str:
    """Return the printed wavelength, ``inf`` for the zero-frequency row."""
    text = cell.strip()
    return "inf" if text == "\N{INFINITY}" else number(text)


def extract() -> dict[str, list[Point]]:
    """Return the mean polarizabilities of Tables I and V as data files."""
    files: dict[str, list[Point]] = {}
    molecule = ""
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        if cells[0]:
            molecule = cells[0]
        if molecule not in GASES:
            continue
        gas = GASES[molecule]
        wavelength = _wavelength(cells[WAVELENGTH])
        point = Point(
            wavelength, wavelength, number(cells[ALPHA]), "", " | ".join(cells)
        )
        files.setdefault(f"{gas}/alms_1975__alpha_{gas}.csv", []).append(point)
    for cells in table(DIRECTORY, TABLE_V)[1:]:
        if cells[0] not in GASES:
            continue
        gas = GASES[cells[0]]
        stuart = number(cells[1].split("(")[0])
        files[f"{gas}/alms_1975__stuart_alpha_{gas}.csv"] = [
            Point("inf", "inf", stuart, "", " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

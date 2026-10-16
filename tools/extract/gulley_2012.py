# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gulley et al. (2012): the fused-silica parameters of their simulations.

Table I (``tbl_763380c3ceb4``), "Terms used to simulate the propagation of the
ultrashort pulse", checked against the table image: the cells agree with the print. Each
value cell carries its unit ("32 fs", "9 eV", "2.48 \times 10^{-16} \text{cm}^2
\text{W}^{-1}"); the script splits number and unit and checks the unit. The GVD row is
printed "361 fs^2 W^-1" in the image itself: a misprinted unit (see the record), kept in
the ``printed`` column. The collision time row reads "Multiple values used" (no value).

Table II (``tbl_f6dc38765622``) is scrambled in the library: one row holding the
characters of each column interleaved. STRUCTURAL FIX: it is transcribed below from the
table image, and the script checks that each scrambled cell is exactly the multiset of
the characters of the transcribed column (spaces dropped). Table II holds FROG-retrieval
diagnostics (G_error and the offset compensation), not material data: no dataset is
written from it; the record quotes the transcription.
"""

import re
from collections import Counter

from _library import Point, number, table, write

DIRECTORY = "Gulley_2012_InteractionUltrashortLaser"
TABLE1 = "tbl_763380c3ceb4"
TABLE2 = "tbl_f6dc38765622"

# Table I symbol -> (file stem, the unit printed after the number).
ROWS = {
    "\\tau_1": ("raman_tau_1", "fs"),
    "\\tau_2": ("raman_tau_2", "fs"),
    "f_r": ("raman_fraction", ""),
    "n_0": ("refractive_index", ""),
    "k_2": ("gvd", "fs^2W^{-1}"),
    "n_2": ("n2", "\\text{cm}^2 \\text{W}^{-1}"),
    "U": ("band_gap", "eV"),
}
N2_CELL = re.compile(r"^([\d.]+) \\times 10\^\{(-?\d+)\} (.*)$")

# Table II as printed (the image): pulse energy (uJ), trace dimensions ("256 x 128",
# with the multiplication sign of the print), G_error, OC.
TRACE = "256 \N{MULTIPLICATION SIGN} 128"
TABLE2_IMAGE = [
    ("1", TRACE, "0.00319", "0.0115"),
    ("10", TRACE, "0.00417", "0.0135"),
    ("17", TRACE, "0.0036", "0.0105"),
    ("25", TRACE, "0.0041", "0.008"),
    ("35", TRACE, "0.00623", "0.01"),
    ("46", TRACE, "0.0051", "0.011"),
    ("50", TRACE, "0.00536", "0.0095"),
    ("55", TRACE, "0.0062", "0.01"),
]


def _check_table2() -> None:
    """Check the transcription of Table II against the scrambled library cells."""
    scrambled = table(DIRECTORY, TABLE2)[2]
    for column, cell in enumerate(scrambled):
        transcribed = "".join(row[column] for row in TABLE2_IMAGE)
        if Counter(cell.replace(" ", "")) != Counter(transcribed.replace(" ", "")):
            msg = f"Table II column {column}: {cell!r} is not {transcribed!r} mixed"
            raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Table I's material values as data files (and check Table II)."""
    _check_table2()
    rows = {cells[0]: cells for cells in table(DIRECTORY, TABLE1)[1:]}
    wavelength_cell = rows["\\lambda"][2]
    if wavelength_cell != "800 nm":
        msg = f"unexpected wavelength {wavelength_cell!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for symbol, (stem, unit) in ROWS.items():
        cells = rows[symbol]
        cell = cells[2]
        printed = " | ".join(cells)
        if symbol == "n_2":
            match = N2_CELL.match(cell)
            if match is None or match.group(2) != "-16" or match.group(3) != unit:
                msg = f"unexpected n2 cell {cell!r}"
                raise ValueError(msg)
            value = number(match.group(1))  # the dataset's scale is 1e-16
        else:
            text, _, printed_unit = cell.partition(" ")
            if printed_unit != unit:
                msg = f"{symbol}: unit {printed_unit!r}, expected {unit!r}"
                raise ValueError(msg)
            value = number(text)
        wavelength = symbol in {"n_0", "k_2", "n_2"}
        x = "800" if wavelength else ""
        files[f"SiO2/gulley_2012__{stem}.csv"] = [
            Point("800" if wavelength else "1", x, value, "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

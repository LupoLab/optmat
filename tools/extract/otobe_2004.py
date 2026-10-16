# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Otobe, Yabana & Iwata (2004): Xe 5p Gamow-state widths in a static field.

Table II (``tbl_a857b4484f15``, p. 053404-6) prints, against the "laser intensity"
(W/cm^2) of the static field, the complex Kohn-Sham orbital energies epsilon_i =
epsilon_i^R + i Gamma_i (eV) of the Xe valence p orbitals p_x, p_y and p_z (the field
along z), as ``-12.27-i1.07[-5]``: ``[k]`` is a power of ten (the caption), so
Gamma_i = -1.07e-5 eV. The headings say "2p"; the orbitals are Xe's 5p.

STRUCTURAL FIX: the library's cells are scrambled. The first data row sits in the
header row ("Laser intensity sW/cm2d431013", "2px seVd -12.27-i1.07f-5g", ...: the
text layer's parentheses and brackets read as s/d and f/g, "4 x 10^13" as "431013"),
and the other five rows are merged into one row, each column's five values
concatenated in one cell. The table is transcribed below from the image, and the script
checks that the library's cells are exactly the transcription in that scrambled form
(a mismatch stops the extraction). The image and the PDF's text layer agree.

One file per orbital, ``tableII_<orbital>_Xe``: the row id and ``x`` are the intensity
as a number (``4e13``); the value is Gamma_i as printed, with its sign (eV); the
dataset's scale turns it into the rate w_i = -2 Gamma_i / hbar in atomic units;
``printed`` keeps the intensity and the complex energy as printed.
"""

import re

from _library import Point, table, write

DIRECTORY = "Otobe_2004_FirstPrinciplesCalculations"
TABLE_II = "tbl_a857b4484f15"
ORBITALS = ("px", "py", "pz")
# Table II as printed (the image, p. 053404-6): intensity (W/cm^2), then the complex
# orbital energies of 2p_x, 2p_y, 2p_z (eV), "[k]" a power of ten.
TABLE_II_IMAGE = (
    ("4", "-12.27-i1.07[-5]", "-12.27-i1.07[-5]", "-12.30-i1.99[-4]"),
    ("5", "-12.33-i2.22[-5]", "-12.33-i2.22[-5]", "-12.38-i7.04[-4]"),
    ("6", "-12.39-i5.13[-5]", "-12.39-i5.13[-5]", "-12.47-i1.85[-3]"),
    ("7", "-12.47-i1.32[-4]", "-12.47-i1.32[-4]", "-12.56-i3.68[-3]"),
    ("8", "-12.56-i2.55[-4]", "-12.56-i2.55[-4]", "-12.69-i8.65[-3]"),
    ("9", "-12.67-i5.93[-4]", "-12.67-i5.93[-4]", "-12.79-i1.11[-2]"),
)
EXPONENT = "13"  # every intensity is printed as n x 10^13 W/cm^2
_ENERGY = re.compile(r"^(-\d+\.\d+)-i(\d+\.\d+)\[(-\d+)\]$")


def _scrambled(energy: str) -> str:
    """Return a printed energy as the library's text layer renders it."""
    return energy.replace("[", "f").replace("]", "g")


def _check_library() -> None:
    """Check the library's scrambled cells against the transcription."""
    grid = table(DIRECTORY, TABLE_II)
    first, rest = TABLE_II_IMAGE[0], TABLE_II_IMAGE[1:]
    expected = [
        [
            f"Laser intensity sW/cm2d{first[0]}310{EXPONENT}",
            *(
                f"2{name} seVd {_scrambled(cell)}"
                for name, cell in zip(ORBITALS, first[1:], strict=True)
            ),
            "IP (eV)12.13",
        ],
        [
            "".join(f"{row[0]}310{EXPONENT}" for row in rest),
            *("".join(_scrambled(row[k]) for row in rest) for k in (1, 2, 3)),
            "",
        ],
    ]
    if grid != expected:
        msg = f"Table II's library cells are not the transcription: {grid!r}"
        raise ValueError(msg)


def _gamma(energy: str) -> str:
    """Return the imaginary part Gamma_i of a printed energy, with its sign, as text."""
    match = _ENERGY.match(energy)
    if match is None:
        msg = f"not a complex energy: {energy!r}"
        raise ValueError(msg)
    return f"-{match.group(2)}e{match.group(3)}"


def extract() -> dict[str, list[Point]]:
    """Return Gamma_i of each p orbital against the intensity."""
    _check_library()
    files: dict[str, list[Point]] = {}
    for k, orbital in enumerate(ORBITALS, start=1):
        files[f"Xe/otobe_2004__tableII_{orbital}_Xe.csv"] = [
            Point(
                f"{row[0]}e{EXPONENT}",
                f"{row[0]}e{EXPONENT}",
                _gamma(row[k]),
                "",
                f"{row[0]} x 10^{EXPONENT} | {row[k]}",
            )
            for row in TABLE_II_IMAGE
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

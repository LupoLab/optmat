# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dodge (1984): MgF2, ordinary and extraordinary indices at 19 degrees C.

Table II (``tbl_722b5ecebf2e``) gives n_O, n_E and the birefringence n_E - n_O at 36
regular wavelengths, 0.20-7.00 um, computed from the Table I Sellmeier formulas (text,
p. 1981): the measured indices themselves are not printed. One file per ray (rows are
the wavelengths in um as printed); the birefringence has no quantity (it is the
difference of the printed columns to the last digit, checked). The cells agree with the
image.

Table I (``tbl_8bf3bb2b0eba``, cells only) holds the constants of the two formulas; this
script checks the record's strings against the cells.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Dodge_1984_RefractivePropertiesMagnesium"
TABLE_1 = "tbl_8bf3bb2b0eba"
TABLE_2 = "tbl_722b5ecebf2e"
PREFIX = "MgF2/dodge_1984__"

# Table I as printed: (label, O ray, E ray).
TABLE_1_PRINTED = (
    ("A_{1}", "0.48755108", "0.41344023"),
    ("A_{2}", "0.39875031", "0.50497499"),
    ("A_{3}", "2.3120353", "2.4904862"),
    ("\\lambda_{1}", "0.04338408", "0.03684262"),
    ("\\lambda_{2}", "0.09461442", "0.09076162"),
    ("\\lambda_{3}", "23.793604", "23.771995"),
)


def _check_table_1() -> None:
    """Check the record's transcription of Table I against the library cells."""
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[1][1:] != ["O ray", "E ray"]:
        msg = f"unexpected Table I header {rows[1]!r}"
        raise ValueError(msg)
    found = tuple(tuple(r) for r in rows[2:])
    if found != TABLE_1_PRINTED:
        msg = f"Table I reads {found}, not {TABLE_1_PRINTED}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Table II's ordinary and extraordinary indices."""
    _check_table_1()
    rows = table(DIRECTORY, TABLE_2)
    if rows[0][1:3] != ["n_O", "n_E"] or len(rows) != 37:
        msg = f"unexpected Table II layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, ray in ((1, "o"), (2, "e")):
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}n_{ray}_19C.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

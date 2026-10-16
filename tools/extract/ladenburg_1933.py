# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ladenburg & Van Voorhis (1933): Ladenburg and Wolfsohn's O2 refractivity.

Table I (``tbl_0010b8557526``, CSV; image checked, the cells agree) gives at nine wave
numbers nu (cm^-1, printed with a thousands comma, "52,084") and the corresponding
wavelengths (A, 1/nu rounded; the first printed "1818A") the contribution G(nu) of the
Schumann-Runge continuum to n - 1 computed from the measured absorption, the O2
dispersion (n - 1) measured by Ladenburg and Wolfsohn (1932), H = (n - 1) - G and the
one-term fit I, each "x 10^7" (the column factor sits in the header, not in the
cells). The 55,000 cm^-1 row has G only.

Only the (n - 1) column is a vocabulary quantity: it is extracted (eight rows), with
the wave number as row id and abscissa (the comma removed). G, H and I are in the
record's notes.
"""

from _library import Point, number, table, write

DIRECTORY = "Ladenburg_1933_ContinuousAbsorptionOxygen"
TABLE_I = "tbl_0010b8557526"
PREFIX = "O2/ladenburg_1933__"

HEADER = [
    "\\nu",
    "\\lambda",
    "G(\\nu) \\times 10^{7}",
    "(n(\\nu)-1) \\times 10^{7}",
    "H(\\nu) \\times 10^{7}",
    "I(\\nu) \\times 10^{7}",
]


def extract() -> dict[str, list[Point]]:
    """Return Ladenburg and Wolfsohn's (n - 1) x 10^7 against wave number."""
    rows = table(DIRECTORY, TABLE_I)
    if rows[0] != HEADER:
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in rows[1:]:
        nu_cell, wavelength, g, refractivity, h, i = (c.strip() for c in cells)
        if not refractivity:
            continue
        nu = number(nu_cell.replace(",", ""))
        points.append(
            Point(
                nu,
                nu,
                number(refractivity),
                "",
                " | ".join((nu_cell, wavelength, g, refractivity, h, i)),
            )
        )
    return {f"{PREFIX}lw_refractivity.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

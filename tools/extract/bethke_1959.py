# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bethke (1959): O2 Schumann-Runge band oscillator strengths.

Table I (``tbl_172bb90e7da2``, cells file only; image checked, the cells agree) lists
for the bands (v', 0), v' = 2-17, the integrated absorption coefficient K_0 (cm^-2
atm^-1), the band oscillator strength f(v', v''), the Franck-Condon factor q, the
electronic oscillator strength f_el and Ditchburn and Heddle's f(v', v''), then the
"Remainder" of the band system to the dissociation limit and the "Total". Every f cell
is "<mantissa> (-5)", the power of ten in parentheses (note c): the mantissa is kept,
the exponent checked to be -5 and taken into the dataset scale. Footnote marks are
stripped: "c" on the first f and on the Ditchburn-Heddle remainder (the image prints
"e" there, note e: their extrapolated f(18,0)+f(19,0)+f(20,0)), "d" on "Remainder".
The image prints the (9, 0) value as 1.07 with a subscript 1 (a guard digit); the cell
reads "1.071", kept.

Files: the present band strengths (``band_f``), the remainder (``remainder``) and the
total (``total``), and Ditchburn and Heddle's band strengths (``ditchburn_heddle``);
the rows are the bands "(v',0)", the abscissa the index v'. K_0, q and f_el are not
extracted (see the record).
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Bethke_1959_OscillatorStrengthsFar"
TABLE_I = "tbl_172bb90e7da2"
PREFIX = "O2/bethke_1959__"

_F = re.compile(r"^([\d.]+)\s*\(-5\)[ce]?$")


def _mantissa(cell: str) -> str:
    match = _F.match(cell.strip())
    if match is None:
        msg = f"unexpected f cell {cell!r}"
        raise ValueError(msg)
    return match.group(1)


def extract() -> dict[str, list[Point]]:
    """Return the band, remainder, total and Ditchburn-Heddle oscillator strengths."""
    band, remainder, total, heddle = [], [], [], []
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        label, k0, f, q, f_el, dh, _ = (c.strip() for c in cells)
        printed = " | ".join((label, k0, f, q, f_el, dh))
        if label.startswith("Remainder"):
            remainder.append(Point("remainder", "18", _mantissa(f), "", printed))
        elif label == "Total":
            total.append(Point("total", "2", _mantissa(f), "", printed))
        else:
            upper, lower = (s.strip() for s in label.split(","))
            if lower != "0":
                msg = f"unexpected band {label!r}"
                raise ValueError(msg)
            row = f"({upper},0)"
            band.append(Point(row, upper, _mantissa(f), "", printed))
            heddle.append(Point(row, upper, _mantissa(dh), "", printed))
    return {
        f"{PREFIX}band_f.csv": band,
        f"{PREFIX}remainder.csv": remainder,
        f"{PREFIX}total.csv": total,
        f"{PREFIX}ditchburn_heddle.csv": heddle,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

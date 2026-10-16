# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Weissler & Lee (1952): O2 absorption coefficients, 303-1306 A.

Table I ("Absorption coefficients of oxygen") runs over two pages:
``tbl_56bc106d3886`` (page 2, 38 rows) and its unlabelled continuation
``tbl_3c7b2b7dceeb`` (page 3, 38 rows), both cells files only (images checked: every
cell agrees). Columns: the wave number (cm^-1), the source line (ion; a leading
``2x``/``3x`` is the grating order), its wavelength (A), k (cm^-1) printed as ``0.07 x
10^2`` and the location in a known band (B.H. head, B.M. middle, B.T. tail, with the
band head's wavelength and the reference).
Rows are numbered in print order (``1`` ... ``76``): the paper prints the line OI
990.797 twice (rows 21 and 24, identical cells; image checked), so the wavelength is
not a unique row id. The abscissa is the source line's wavelength (the wave-number
column is 1e8/lambda, kept in ``printed``); the value is the printed mantissa, the
``x 10^2`` column factor going to the dataset's scale. Some wave numbers carry italic
final digits in the image (82891.2, 83364.8, 90090.0, 101030, 109170, 119746, 188900),
apparently marking digits not significant for the rounder wavelengths; italics are lost
in the cells and not needed.

Table II (``tbl_3deee1da1d4e``, CSV; image checked) compares k from first- and
second-order lines on plate LE 25: the wavelength cell holds the ion and the wavelength
("OII 644.188"), the values carry the unit ("640 \text{ cm}^{-1}"), which is stripped.
The percent-difference column is derived and not extracted.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Weissler_1952_AbsorptionCoefficientsOxygen"
TABLE_I = ("tbl_56bc106d3886", "tbl_3c7b2b7dceeb")
TABLE_II = "tbl_3deee1da1d4e"
PREFIX = "O2/weissler_1952__"

# "0.07 \times 10^{2}" (page 2 cells) or "6.7 <multiplication sign>102" (page 3
# cells): the mantissa.
_K = re.compile(r"^([\d.]+)\s*(?:\\times\s*10\^\{2\}|\N{MULTIPLICATION SIGN}\s*102)$")
_UNIT = re.compile(r"^(\d+)\s*\\text\{\s*cm\}\^\{-1\}$")


def _table_i() -> list[Point]:
    points = []
    for table_id in TABLE_I:
        for cells in json_table(DIRECTORY, table_id)[1:]:
            wavenumber, ion, wavelength, k, band, reference = (c.strip() for c in cells)
            match = _K.match(k)
            if match is None:
                msg = f"unexpected k cell {k!r}"
                raise ValueError(msg)
            float(wavelength)
            points.append(
                Point(
                    str(len(points) + 1),
                    wavelength,
                    match.group(1),
                    "",
                    " | ".join((wavenumber, ion, wavelength, k, band, reference)),
                )
            )
    return points


def _table_ii() -> tuple[list[Point], list[Point]]:
    first, second = [], []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        line, one, two, percent = (c.strip() for c in cells)
        wavelength = line.split()[-1]
        float(wavelength)
        printed = " | ".join((line, one, two, percent))
        for column, cell in ((first, one), (second, two)):
            match = _UNIT.match(cell)
            if match is None:
                msg = f"unexpected Table II cell {cell!r}"
                raise ValueError(msg)
            column.append(Point(wavelength, wavelength, match.group(1), "", printed))
    return first, second


def extract() -> dict[str, list[Point]]:
    """Return Table I and the two orders of Table II."""
    first, second = _table_ii()
    return {
        f"{PREFIX}table_1.csv": _table_i(),
        f"{PREFIX}first_order.csv": first,
        f"{PREFIX}second_order.csv": second,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

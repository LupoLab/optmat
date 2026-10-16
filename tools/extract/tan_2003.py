# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tan, Lemon & French (2003): the band gaps of Suprasil 1 silica (Table I).

Table I (``tbl_1982492688af``, cells only, no CSV) lists, for three ranges of the
absorption coefficient, the band gap from a direct-gap fit (alpha^2 E^2 linear in E)
and from an indirect-gap fit (alpha^(1/2) linear in E), each with its fitting range in
cm^-1. Two data files, direct and indirect, one row per fitting range; the row label
(the dataset's ``x``) is the range as printed, rewritten in ASCII ("2.5e5-5.6e5 cm-1").
The dagger on the first row's gaps ("Calculated from ellipsometric spectroscopy") is
stripped from the value and kept in ``printed``.

Checked against the table image (2026-09-27): all cells agree. The first column prints
"Suprasil 1" on the first row and "Amorphous SiO2" on the second (the third is blank):
one sample, Suprasil 1, throughout (the caption).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Tan_2003_OpticalPropertiesLondon"
TABLE = "tbl_1982492688af"
PREFIX = "SiO2/tan_2003__"
EN_DASH = "\N{EN DASH}"
_POWER = re.compile(r"\s*\\times\s*10\^\{(\d+)\}")


def _range(cell: str) -> str:
    """Return a printed fitting range (two numbers and an en dash) in ASCII."""
    low, high = (_POWER.sub(r"e\1", part).strip() for part in cell.split(EN_DASH))
    float(low)
    float(high)
    return f"{low}-{high} cm-1"


def extract() -> dict[str, list[Point]]:
    """Return the direct and the indirect band gaps, one row per fitting range."""
    rows = json_table(DIRECTORY, TABLE)
    if rows[0][1] != "Direct gap (eV)" or rows[0][3] != "Indirect gap (eV)":
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    direct: list[Point] = []
    indirect: list[Point] = []
    for index, cells in enumerate(rows[1:], start=1):
        printed = " | ".join(cells)
        for points, gap, fit_range in ((direct, 1, 2), (indirect, 3, 4)):
            value = number(cells[gap].removesuffix("^{\\dagger}"))
            label = _range(cells[fit_range])
            points.append(Point(str(index), label, value, "", printed))
    return {f"{PREFIX}direct_gap.csv": direct, f"{PREFIX}indirect_gap.csv": indirect}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

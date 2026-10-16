# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Häusler and Kerl (1988): the molar refractivity A_1 of SF6 (Table I).

Only SF6 is a curated species; C2H4 and C2H6 are not extracted.

Table I (pp. 122-123) is split by the library into ``tbl_e56cd39cc6f5`` (T, p_max,
K, w/uw, A_1, A_2; no CSV rendering, read from its JSON cells) and the unlabelled
``tbl_2a63c5387037`` (A_3, A_4, alpha_0, -B, C), which has lost its row labels. Only
the first part is used: A_1 (cm^3/mol) is the molar refractivity of the free molecule,
and the second part's alpha_0 = 3 A_1/(4 pi N_A) repeats it. Each isotherm has a
weighted (w) and an unweighted (uw) row; the standard deviation of the last digits is
in parentheses, after a space ("10.9346 (32)"). The temperature cell of the isotherm
marked ``^a`` ("298.32^a (1)") was measured at 546.22 nm (the table's note a; the text,
p. 124, says 546.226 nm), the others at 632.99 nm. Rows and values checked against the
table image.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Husler_1988_MeanPolarizabilitiesSecond"
TABLE_I = "tbl_e56cd39cc6f5"
HE_NE = "632.99"  # nm, vacuum (p. 118)
HG = "546.22"  # nm, note a of Table I
_TEMPERATURE = re.compile(r"^(?P<value>[\d.]+)(?P<hg>\^a)? \(\d+\)$")
_CELL = re.compile(r"^(?P<value>[\d.]+) \((?P<sigma>\d+)\)$")


def extract() -> dict[str, list[Point]]:
    """Return one A_1 file per SF6 isotherm and weighting."""
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_I)
    start = next(i for i, cells in enumerate(rows) if cells[0] == "SF_6")
    name = wavelength = ""
    for cells in rows[start + 1 :]:
        if cells[0]:  # the weighted row opens an isotherm
            match = _TEMPERATURE.match(cells[0].strip())
            if match is None:
                msg = f"unexpected temperature {cells[0]!r}"
                raise ValueError(msg)
            wavelength = HG if match["hg"] else HE_NE
            suffix = "K_hg" if match["hg"] else "K"
            name = match["value"].replace(".", "_") + suffix
        weighting = cells[3].strip()
        if weighting not in {"w", "uw"}:
            msg = f"unexpected weighting {cells[3]!r}"
            raise ValueError(msg)
        match = _CELL.match(cells[4].strip())
        if match is None:
            msg = f"unexpected A_1 cell {cells[4]!r}"
            raise ValueError(msg)
        value = number(match["value"])
        decimals = len(value.partition(".")[2])
        sigma = f"{int(match['sigma']) * 10.0**-decimals:.{decimals}f}"
        point = Point(wavelength, wavelength, value, sigma, " | ".join(cells))
        files[f"SF6/hausler_1988__A_1_{weighting}_{name}_SF6.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

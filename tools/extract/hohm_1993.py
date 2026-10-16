# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hohm (1993), Mol. Phys. 78, 929: CH4 polarizabilities and Cauchy moments.

Only CH4 is a curated species; the other hydrocarbons are not extracted.

Table 1 (``tbl_90f137181bd7``, no CSV rendering: read from its JSON cells, checked
against the image, p. 935) gives alpha(omega_k) in atomic units at the five vacuum
wavelengths of its note a) (632.990, 611.971, 594.096, 543.516, 325.130 nm): row b
this work, with the uncertainty of the last digits in parentheses; row c [26], Thomas
and Meath's DOSD. The decimal point is printed "·".

Table 2 (``tbl_efe74204a8d8``, p. 937) prints, per gas, the fit parameters a_0-a_2 of
Eq. (13) (recorded as a formula, not here) and the Cauchy moments S(-2) ... S(-8) in
three stacked rows: this work (iterated), b [26] Thomas and Meath, e [13] Spackman's
TDHF. The library merges each stack into one cell ("63·2262·4150·07" for S(-4) =
63.22, 62.41, 50.07). The split below was read from the table image; the script checks
that the parts rejoin to the extracted cell exactly, so it cannot outlive a change in
the extraction.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hohm_1993_ExperimentalDeterminationDispersion"
TABLE_1 = "tbl_90f137181bd7"
TABLE_2 = "tbl_efe74204a8d8"
WAVELENGTHS = ("632.990", "611.971", "594.096", "543.516", "325.130")  # note a), nm
ORDERS = ("-2", "-4", "-6", "-8")
# CH4's merged S(-2k) cells, split as the image shows them (this work, [26], [13]);
# the S(-8) cell stacks only two values: Spackman [13] gives no S(-8).
SPLITS = {
    "17·257^{b}17\\cdot27^{e}15\\cdot895": ("17·257", "17·27", "15·895"),
    "63·2262·4150·07": ("63·22", "62·41", "50·07"),
    "293·8298·3199·11": ("293·8", "298·3", "199·11"),
    "15051714": ("1505", "1714"),
}
ROWS = ("this_work", "thomas_1977", "spackman_1991")
_CELL = re.compile(r"^(?:\^[a-z])?(?P<value>[\d.]+)(?:\((?P<sigma>\d+)\))?$")


def _split(cell: str) -> tuple[str, ...]:
    """Return the stacked values of a merged Table 2 cell, as printed."""
    parts = SPLITS[cell]
    rejoined = "".join(parts)
    plain = cell.replace("^{b}", "").replace("^{e}", "").replace("\\cdot", "·")
    if rejoined != plain:
        msg = f"Table 2 cell {cell!r} no longer splits as {parts}"
        raise ValueError(msg)
    return parts


def _table_1() -> dict[str, list[Point]]:
    """Return the CH4 rows of Table 1: this work and Thomas and Meath's DOSD."""
    rows = json_table(DIRECTORY, TABLE_1)
    measured, dosd = rows[1], rows[2]
    if measured[0] != "CH_4" or not measured[2].startswith("^b") or dosd[2][:2] != "^c":
        msg = "Table 1 no longer starts with the CH4 rows b and c"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for name, cells in (("alpha_CH4", measured), ("thomas_1977_CH4", dosd)):
        points = []
        printed = " | ".join(cells)
        for wavelength, cell in zip(WAVELENGTHS, cells[2:7], strict=True):
            match = _CELL.match(cell.replace("·", ".").strip())
            if match is None:
                msg = f"unexpected Table 1 cell {cell!r}"
                raise ValueError(msg)
            value = number(match["value"])
            sigma = ""
            if match["sigma"]:
                decimals = len(value.partition(".")[2])
                sigma = f"{int(match['sigma']) * 10.0**-decimals:.{decimals}f}"
            points.append(Point(wavelength, wavelength, value, sigma, printed))
        files[f"CH4/hohm_1993__{name}.csv"] = points
    return files


def _table_2() -> dict[str, list[Point]]:
    """Return CH4's Cauchy moments of Table 2: this work, [26] and [13]."""
    cells = json_table(DIRECTORY, TABLE_2)[1]
    if cells[0] != "CH_{4}":
        msg = "Table 2 no longer starts with CH4"
        raise ValueError(msg)
    printed = " | ".join(cells)
    points: dict[str, list[Point]] = {row: [] for row in ROWS}
    for order, cell in zip(ORDERS, cells[4:8], strict=True):
        for row, text in zip(ROWS, _split(cell), strict=False):
            value = number(text.replace("·", "."))
            points[row].append(Point(order, order, value, "", printed))
    return {
        "CH4/hohm_1993__moments_CH4.csv": points["this_work"],
        "CH4/hohm_1993__moments_thomas_1977_CH4.csv": points["thomas_1977"],
        "CH4/hohm_1993__moments_spackman_1991_CH4.csv": points["spackman_1991"],
    }


def extract() -> dict[str, list[Point]]:
    """Return the CH4 polarizabilities and moments."""
    return _table_1() | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

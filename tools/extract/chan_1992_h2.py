# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1992): H2 Lyman and Werner band oscillator strengths.

All four tables are cells files only (merged header cells); the images agree with the
cells (checked 2026-09-27).

- Table 1 (``tbl_aae134e5c906``, unlabelled in the extraction; the caption is "Absolute
  oscillator strengths for the vibronic transitions of the Lyman band"): f(v'' = 0 ->
  v') of B 1Sigma_u+ <- X 1Sigma_g+, v' = 0-22, one column per source (the present
  high-resolution dipole (e,e) work, Geiger and Schmoranzer, Fabian and Lewis, Hesser et
  al., Haddad et al., Allison and Dalgarno). Rows and the label abscissa are v'. The
  present values at v' = 14 and 20 are printed in parentheses with mark a,
  "(0.0101)^a", "Interpolated values": the parentheses and mark are stripped here and
  the rows are excluded in the record.
- Table 2 (``tbl_254f5e9a82b9``): the Werner band, C 1Pi_u <- X, v' = 0-6, as Table 1
  (four sources). The present v' = 3 value carries mark a ("0.0544a)": the underlying
  Lyman v' = 14 component subtracted), stripped here.
- Table 3 (``tbl_5acc9741b87e``): total integrated oscillator strengths of the two
  bands, one file per reference (rows ``lyman`` and ``werner``, label abscissa). The
  "(i) TDA", "(ii) RPA", "(i) dipole length", "(ii) dipole velocity" rows are variants
  of the reference named on the row above.

Table 4 (dipole strengths D_e(r_0), au) has no quantity and is not extracted (source
notes).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Chan_1992_AbsoluteOpticalOscillator_095875"
LYMAN = "tbl_aae134e5c906"
WERNER = "tbl_254f5e9a82b9"
TOTALS = "tbl_5acc9741b87e"
PREFIX = "H2/chan_1992_h2__"
# Table 1 and 2 column headers -> file stem.
BAND_SOURCES = {
    "present work": "present",
    "Geiger and Schmoranzer [57]": "geiger_1969",
    "Fabian and Lewis [44]": "fabian_1974",
    "Hesser et al. [43]": "hesser_1968",
    "Haddad et al. [41]": "haddad_1968",
    "Allison and Dalgarno [12,26]": "allison_1970",
    "theory Allison and Dalgarno [12,26]": "allison_1970",
}
# Table 3 reference -> file stem; the variant rows are appended to their parent.
TOTAL_SOURCES = {
    "Arrighini et al. [27]": "arrighini_1980",
    "Allison and Dalgarno [12,26]": "allison_1970",
    "Browne [24]": "browne_1969",
    "Rothenberg and Davidson [22]": "rothenberg_1967",
    "Miller and Krauss [21]": "miller_1967",
    "Peek and Lassettre [20]": "peek_1963",
    "Ehrenson and Phillipson [19]": "ehrenson_1961",
    "Shull [18]": "shull_1952",
    "Mulliken and Rieke [17]": "mulliken_1941",
    "present work (dipole (e, e))": "present",
    "Geiger and Schmoranzer [57] (electron impact)": "geiger_1969",
    "Hesser et al. [43] (optical: curve of growth)": "hesser_1968",
    "Hesser [49] (lifetimes)": "hesser_1968b",
}
VARIANTS = {
    "(i) TDA": "tda",
    "(ii) RPA": "rpa",
    "(i) dipole length": "length",
    "(ii) dipole velocity": "velocity",
}
_MARKED = re.compile(r"^\(?(?P<value>\d*\.\d+)\)?(?:\^?a\)?)?$")


def _band(table_id: str, band: str) -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, table_id)
    # The source names sit in header row 2, except Allison and Dalgarno in Table 1,
    # whose merged header cell spans rows 1-2 (anchored in row 1).
    headers = [rows[2][c] or rows[1][c] for c in range(len(rows[0]))]
    files: dict[str, list[Point]] = {}
    for cells in rows[3:]:
        level = cells[0].strip()
        printed = " | ".join(cells)
        for column in range(1, len(cells)):
            cell = cells[column].replace(" ", "")
            if not cell:
                continue
            match = _MARKED.match(cell)
            if match is None:
                msg = f"unreadable cell {cell!r}"
                raise ValueError(msg)
            stem = BAND_SOURCES[headers[column]]
            point = Point(level, level, number(match["value"]), "", printed)
            files.setdefault(f"{PREFIX}{band}_{stem}.csv", []).append(point)
    return files


def _totals() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    parent = ""
    for cells in json_table(DIRECTORY, TOTALS)[2:]:
        label = cells[1].strip()
        if label in VARIANTS:
            stem = f"{parent}_{VARIANTS[label]}"
        else:
            stem = parent = TOTAL_SOURCES[label]
        printed = " | ".join(cells)
        for column, band in ((2, "lyman"), (3, "werner")):
            if cells[column].strip():
                point = Point(band, band, number(cells[column]), "", printed)
                files.setdefault(f"{PREFIX}total_{stem}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the band oscillator strengths (Tables 1, 2) and totals (Table 3)."""
    return _band(LYMAN, "lyman") | _band(WERNER, "werner") | _totals()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

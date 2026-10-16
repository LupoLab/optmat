# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1993): O2 integrated oscillator strengths, 9.75-18 eV.

Both tables are unlabelled in the extraction (Tables 1 and 2 of the paper) and cells
files only; the images agree with the cells (checked 2026-09-27). The paper tabulates
no df/dE (its high-resolution spectrum is in Figs. 1-3 only).

- Table 1 (``tbl_64a7a5d7efc8``): the diffuse bands 9.75-12.07 eV. The experimental
  columns (present work, Huebner et al., Ogawa and Ogawa) are oscillator strengths
  integrated over the printed energy ranges: one ``interval_oscillator_strength`` file
  per column, row = the range as printed (en dash turned into "-"), x and x_upper its
  edges. Lewis et al.'s two values (isotopic O2, one averaged over 10.19-10.35 eV,
  mark b stripped) and the three theories' band values (Buenker et al., Yoshimine et
  al., Li et al.) belong to the assigned bands, not to the ranges: they go to
  ``oscillator_strength`` files with a label abscissa, rows ``E0``, ``E1``, ``E2`` (E
  3Sigma_u-, v' = 0-2) and ``Pi0``, ``Pi1``, ``Pi2`` (2 3Pi_u, v' = 0-2). The
  assignment of a row whose first cell is blank is the one above it only where the
  first cell is a bare v' ("1", "2"); otherwise the row is unassigned.
- Table 2 (``tbl_41c9c9856c23``): intervals 12.070-18.287 eV (present work, Huebner et
  al.), as Table 1's experimental columns.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Chan_1993_AbsoluteOpticalOscillator_374fef"
TABLE_1 = "tbl_64a7a5d7efc8"
TABLE_2 = "tbl_41c9c9856c23"
PREFIX = "O2/chan_1993_o2__"
INTERVAL_SOURCES = {
    "present work": "present",
    "Huebner et al. [28]": "huebner_1975",
    "Ogawa and Ogawa [14]": "ogawa_1975",
}
BAND_SOURCES = {
    "Lewis et al. [23]": "lewis_1988",
    "Buenker et al. [20]": "buenker_1976",
    "Yoshimine et al. [19]": "yoshimine_1976",
    "Li et al. [21]": "li_1992",
}
_STATE = re.compile(r"^(?P<state>E|2)\^\{3\}.*v' = (?P<v>\d)$")
_VALUE = re.compile(r"^(?P<value>\d*\.\d+)(?:\^\{b\)\})?$")


def _interval(printed: str) -> tuple[str, str, str]:
    interval = printed.replace("\N{EN DASH}", "-")
    low, high = interval.split("-")
    return interval, number(low), number(high)


def _band(label: str, previous: str) -> str:
    """Return the band row id of Table 1's assignment cell ('' for none)."""
    found = _STATE.match(label)
    if found:
        return ("E" if found["state"] == "E" else "Pi") + found["v"]
    if label.isdigit() and previous:
        return previous.rstrip("0123456789") + label
    return ""


def _table_1() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_1)
    headers = [rows[3][c] or rows[2][c] for c in range(len(rows[0]))]
    files: dict[str, list[Point]] = {}
    band = ""
    for cells in rows[4:]:
        band = _band(cells[0], band)
        interval, low, high = _interval(cells[1])
        printed = " | ".join(cells)
        for column in range(2, len(cells)):
            if not cells[column]:
                continue
            match = _VALUE.match(cells[column].replace(" ", ""))
            if match is None:
                msg = f"unreadable cell {cells[column]!r}"
                raise ValueError(msg)
            value = number(match["value"])
            if headers[column] in INTERVAL_SOURCES:
                stem = INTERVAL_SOURCES[headers[column]]
                point = Point(interval, low, value, "", printed, {"x_upper": high})
                files.setdefault(f"{PREFIX}intervals_{stem}.csv", []).append(point)
            else:
                stem = BAND_SOURCES[headers[column]]
                point = Point(band, band, value, "", printed)
                files.setdefault(f"{PREFIX}bands_{stem}.csv", []).append(point)
    return files


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        interval, low, high = _interval(cells[0])
        printed = " | ".join(cells)
        for column in (1, 2):
            if cells[column]:
                stem = INTERVAL_SOURCES[rows[1][column]]
                point = Point(
                    interval, low, number(cells[column]), "", printed, {"x_upper": high}
                )
                files.setdefault(f"{PREFIX}intervals_{stem}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1 and 2: interval and band oscillator strengths."""
    files = _table_1()
    for path, points in _table_2().items():
        files.setdefault(path, []).extend(points)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

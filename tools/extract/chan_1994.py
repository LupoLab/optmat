# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1994): N2O absolute oscillator strengths, 5.5-203 eV.

All three tables are unlabelled in the extraction (Tables 1-3 of the paper); images
checked against the cells (2026-09-27).

- Table 1 (``tbl_b536b17df805``, CSV): the low-resolution (1 eV FWHM) df/dE in
  ``10^-2 eV^-1`` at 7.40-203.00 eV, four column pairs read pair by pair. Rows are the
  energies as printed; two of them are misprinted (errata in the record).
- Table 2 (``tbl_c3cbc336c3cb``, cells file): the three dissociative continuum bands
  below 10 eV, each integrated over the energy range printed in its header (B 1Delta
  5.70-7.70, C 1Pi 7.70-8.97, D 1Sigma+ 8.97-10.23 eV). The present values and Lee
  and Suto's (integrated by Chan et al. over the same ranges from Lee and Suto's
  digitised figure, p. 9) go to ``interval_oscillator_strength`` files, row = the
  range (the present rows open the file of Table 3's present intervals, which they
  continue); the other sources' values are their own band strengths over ranges not
  stated here (Zelikoff et al.'s differ from these), so they go to
  ``oscillator_strength`` files with rows ``B``, ``C``, ``D``. The theory row prints
  its values beside "theory" and the reference, "Chutjian and Segal [5]", on the next
  (empty) row (image): the values belong to Chutjian and Segal.
- Table 3 (``tbl_989c4d1b6636``, cells file): integrated oscillator strengths over 35
  intervals, 10.23-16.40 eV, present work and "ref. [33]" (Huebner et al.).
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1994_DiscreteContinuumPhotoabsorption"
TABLE_1 = "tbl_b536b17df805"
TABLE_2 = "tbl_c3cbc336c3cb"
TABLE_3 = "tbl_989c4d1b6636"
PREFIX = "N2O/chan_1994__"
REFERENCES = {
    "1": "lee_1984",
    "3": "zelikoff_1953",
    "5": "chutjian_1972",
    "16": "rabalais_1971",
    "33": "huebner_1975",
}
# Table 2 sources whose values are integrals over the header's ranges.
SAME_RANGES = {"present", "lee_1984"}
_REFERENCE = re.compile(r"\[(\d+)\]")
_RANGE = re.compile(r"\((\d+\.\d+)-(\d+\.\d+) eV\)")


def _stem(label: str) -> str:
    if label.startswith("present"):
        return "present"
    found = _REFERENCE.search(label)
    if found is None:
        msg = f"no reference in {label!r}"
        raise ValueError(msg)
    return REFERENCES[found[1]]


def _interval_point(interval: str, value: str, printed: str) -> Point:
    low, high = interval.split("-")
    return Point(interval, low, number(value), "", printed, {"x_upper": high})


def _table_1() -> list[Point]:
    rows = table(DIRECTORY, TABLE_1)[1:]
    return [
        Point(
            number(cells[2 * k]),
            number(cells[2 * k]),
            number(cells[2 * k + 1]),
            "",
            f"{cells[2 * k].strip()} | {cells[2 * k + 1].strip()}",
        )
        for k in range(4)
        for cells in rows
        if cells[2 * k].strip()
    ]


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    intervals = []
    for cell in rows[2][1:]:
        found = _RANGE.search(cell.replace("\N{EN DASH}", "-"))
        assert found is not None, cell
        intervals.append(f"{found[1]}-{found[2]}")
    files: dict[str, list[Point]] = {}
    for index in range(3, len(rows)):
        cells = rows[index]
        if not any(cells[1:]):
            continue  # "Chutjian and Segal [5]" (read with the row above), "experiment"
        # The theory values sit on the "theory" row; the reference is on the next row.
        label = rows[index + 1][0] if cells[0] == "theory" else cells[0]
        stem = _stem(label)
        printed = " | ".join(cells) + ("" if label == cells[0] else f" || {label}")
        for interval, band, cell in zip(intervals, "BCD", cells[1:], strict=True):
            if not cell:
                continue
            if stem in SAME_RANGES:
                point = _interval_point(interval, cell, printed)
            else:
                point = Point(band, band, number(cell), "", printed)
            # The present bands continue Table 3's intervals: one dataset.
            name = "intervals" if stem == "present" else "bands"
            files.setdefault(f"{PREFIX}{name}_{stem}.csv", []).append(point)
    return files


def _table_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_3)
    for cells in rows[2:]:
        interval = cells[0].replace("\N{EN DASH}", "-")
        printed = " | ".join(cells)
        for column in (1, 2):
            if cells[column]:
                stem = _stem(rows[1][column])
                point = _interval_point(interval, cells[column], printed)
                files.setdefault(f"{PREFIX}intervals_{stem}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table 1 (df/dE), Table 2 (bands) and Table 3 (intervals)."""
    files = {f"{PREFIX}continuum.csv": _table_1()} | _table_2()
    for path, points in _table_3().items():
        files.setdefault(path, []).extend(points)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

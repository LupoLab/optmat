# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1993): CO2 absolute oscillator strengths, 6-203 eV.

Images checked against the cells (2026-09-27).

- Table 1 (``tbl_b45cfc94562b``, CSV, unlabelled): the low-resolution (1 eV FWHM) df/dE
  in ``10^-2 eV^-1`` at 8.53-203.00 eV, four column pairs read pair by pair.
- Table 2 (``tbl_914f5351589b``): the two weak bands below 10 eV (1Delta_u, 1Pi_g). The
  extraction labels this table "Table 3" and puts "Table 2" and the 7-10 eV caption
  in its footnotes (a detached caption, structural): it is the paper's Table 2 (image
  and text, p. 8). One file per row (source), rows ``Delta_u`` and ``Pi_g``, label
  abscissa.
- Table 3 (``tbl_ae9d1bb0325b``, unlabelled): the 1Sigma_u+ (10.30-11.28 eV) and 1Pi_u
  (11.28-11.49 eV) bands. The present values are integrals over the printed ranges
  (p. 9): an ``interval_oscillator_strength`` file with rows ``10.30-11.28`` and
  ``11.28-11.49``. The other rows (experiments and theories) give the states'
  strengths without saying over which range: ``oscillator_strength`` files with rows
  ``Sigma_u`` and ``Pi_u``.
- Table 4 (``tbl_c908d6e74402``, CSV, unlabelled): integrated oscillator strengths over
  16 intervals from 11.49 eV to the first ionisation threshold (13.78 eV).

The method in parentheses after each source ("(photoabsorption)") is not part of the
file name.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1993_ElectronicSpectrumCarbon"
TABLE_1 = "tbl_b45cfc94562b"
TABLE_2 = "tbl_914f5351589b"
TABLE_3 = "tbl_ae9d1bb0325b"
TABLE_4 = "tbl_c908d6e74402"
PREFIX = "CO2/chan_1993_co2__"
REFERENCES = {
    "9": "mccurdy_1974",
    "15": "padial_1981",
    "17": "wilkinson_1950",
    "18": "inn_1953",
    "25": "rabalais_1971",
    "39": "meyer_1965",
    "43": "klump_1978",
}
_REFERENCE = re.compile(r"\[(\d+)\]")
_RANGE = re.compile(r"\((\d+\.\d+)-(\d+\.\d+)\)")


def _stem(label: str) -> str:
    if label.startswith("present"):
        return "present"
    found = _REFERENCE.search(label)
    if found is None:
        msg = f"no reference in {label!r}"
        raise ValueError(msg)
    return REFERENCES[found[1]]


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
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_2)[2:]:
        printed = " | ".join(cells)
        for column, state in ((1, "Delta_u"), (2, "Pi_g")):
            if cells[column]:
                point = Point(state, state, number(cells[column]), "", printed)
                path = f"{PREFIX}bands_7_10_{_stem(cells[0])}.csv"
                files.setdefault(path, []).append(point)
    return files


def _table_3() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_3)
    ranges = [_RANGE.search(rows[1][column]) for column in (1, 2)]
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        if not any(cells[1:]):
            continue  # the "Theory" and "Experiment" headings
        stem, printed = _stem(cells[0]), " | ".join(cells)
        for column, state, found in zip(
            (1, 2), ("Sigma_u", "Pi_u"), ranges, strict=True
        ):
            if not cells[column]:
                continue
            value = number(cells[column])
            if stem == "present":
                assert found is not None
                low, high = found[1], found[2]
                point = Point(
                    f"{low}-{high}", low, value, "", printed, {"x_upper": high}
                )
                files.setdefault(f"{PREFIX}intervals_present.csv", []).append(point)
            else:
                point = Point(state, state, value, "", printed)
                path = f"{PREFIX}bands_10_11_{stem}.csv"
                files.setdefault(path, []).append(point)
    return files


def _table_4() -> list[Point]:
    points = []
    for cells in table(DIRECTORY, TABLE_4)[1:]:
        interval = cells[0].strip().replace("\N{EN DASH}", "-")
        low, high = interval.split("-")
        points.append(
            Point(
                interval,
                low,
                number(cells[1]),
                "",
                " | ".join(cells),
                {"x_upper": high},
            )
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return Table 1 (df/dE) and Tables 2-4 (bands and intervals)."""
    files = {f"{PREFIX}continuum.csv": _table_1()} | _table_2() | _table_3()
    files[f"{PREFIX}intervals_present.csv"].extend(_table_4())
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

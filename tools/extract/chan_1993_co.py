# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1993): CO absolute oscillator strengths, 7-200 eV.

Images checked against the cells (2026-09-27).

- Table 1 (``tbl_ca93596af562``, CSV): the low-resolution (1 eV FWHM) df/dE in
  ``10^-2 eV^-1`` at 7-200 eV, four column pairs read pair by pair (print order).
  Rows are the energies as printed (unique).
- Table 2 (``tbl_889679627bd2``, unlabelled in the extraction; its caption is the A 1Pi
  vibronic-band table): f(v'' = 0 -> v') of X 1Sigma+ -> A 1Pi, v' = 0-12, one file per
  column (source), rows and label v'. Rich's column header (a merged cell spanning two
  header rows) sits in row 1, the theory headers likewise.
- Table 3 (``tbl_05bec4dbd2ec``, unlabelled): the total A 1Pi oscillator strength, one
  file per reference, one row ``A``.
- Table 4 (``tbl_0aedd197bfd4``): B 1Sigma+, C 1Sigma+ and E 1Pi, v' = 0 and 1, rows
  ``B0`` ... ``E1``. The image prints the values of Nielsen et al., Padial et al.,
  Wood, Coughran et al. and Rose et al. centred between the v' = 0 and v' = 1 columns
  of each state: they are the state's summed strength (the text, p. 11: 'should be
  almost equal to the sum of the oscillator strengths for the v' = 0 and 1 bands');
  the cells put them under v' = 1. Structural fix: for those references each value
  goes to a ``*_sum`` file with rows ``B``, ``C``, ``E``.
- Table 5 (``tbl_d446b685ddb3``, unlabelled): integrated oscillator strengths over
  energy intervals, 12.13-16.98 eV (present work, Stark et al., Letzelter et al.); the
  row id is the printed interval, x its lower and x_upper its upper edge (eV). Two
  printed interval edges are misprinted and corrected by errata in the record.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1993_AbsoluteOpticalOscillator"
TABLE_1 = "tbl_ca93596af562"
TABLE_2 = "tbl_889679627bd2"
TABLE_3 = "tbl_05bec4dbd2ec"
TABLE_4 = "tbl_0aedd197bfd4"
TABLE_5 = "tbl_d446b685ddb3"
PREFIX = "CO/chan_1993_co__"
REFERENCES = {
    "7": "rose_1973",
    "8": "coughran_1973",
    "9": "wood_1974",
    "10": "nielsen_1980",
    "11": "lynch_1982",
    "12": "padial_1978",
    "14": "kirby_1989",
    "15": "chantranupong_1992",
    "35": "lee_1981",
    "36": "letzelter_1987",
    "39": "stark_1991",
    "40": "eidelsberg_1992",
    "43": "lassettre_1971",
    "44": "rich_1968",
    "45": "hesser_1968",
    "59": "field_1983",
}
SUMMED = {"nielsen_1980", "padial_1978", "wood_1974", "coughran_1973", "rose_1973"}
STATES = {1: "B0", 2: "B1", 3: "C0", 4: "C1", 5: "E0", 6: "E1"}
_REFERENCE = re.compile(r"\[(\d+)\]")


def _stem(label: str) -> str:
    """Return the file stem of a column or row label (the present work or a ref.)."""
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
    rows = json_table(DIRECTORY, TABLE_2)
    headers = [rows[2][c] or rows[1][c] for c in range(len(rows[0]))]
    files: dict[str, list[Point]] = {}
    for cells in rows[3:]:
        level, printed = cells[0], " | ".join(cells)
        for column in range(1, len(cells)):
            if cells[column]:
                stem = _stem(headers[column])
                point = Point(level, level, number(cells[column]), "", printed)
                files.setdefault(f"{PREFIX}a_state_{stem}.csv", []).append(point)
    return files


def _table_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_3)[1:]:
        if not cells[1]:
            continue  # the "Theory" and "Experiment" headings
        # The Lassettre and Skerbele label is garbled in the extraction
        # ("Lasscttce and Skcrbcle [43]"); its reference number is intact.
        stem = _stem(cells[0])
        point = Point("A", "A", number(cells[1]), "", " | ".join(cells))
        files[f"{PREFIX}a_total_{stem}.csv"] = [point]
    return files


def _table_4() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_4)[3:]:
        if not any(cells[1:]):
            continue  # the "Theory" and "Experiment" headings
        stem, printed = _stem(cells[0]), " | ".join(cells)
        for column, state in STATES.items():
            if not cells[column]:
                continue
            if stem in SUMMED:
                # Printed centred between v' = 0 and 1: the state's sum (docstring).
                point = Point(state[0], state[0], number(cells[column]), "", printed)
                files.setdefault(f"{PREFIX}bce_{stem}_sum.csv", []).append(point)
            else:
                point = Point(state, state, number(cells[column]), "", printed)
                files.setdefault(f"{PREFIX}bce_{stem}.csv", []).append(point)
    return files


def _table_5() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_5)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        interval = cells[0].replace("\N{EN DASH}", "-")
        low, high = interval.split("-")
        printed = " | ".join(cells)
        for column in range(1, 4):
            if cells[column]:
                stem = _stem(rows[1][column])
                point = Point(
                    interval,
                    number(low),
                    number(cells[column]),
                    "",
                    printed,
                    extra={"x_upper": number(high)},
                )
                files.setdefault(f"{PREFIX}intervals_{stem}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table 1 (df/dE) and Tables 2-5 (discrete and interval strengths)."""
    return (
        {f"{PREFIX}continuum.csv": _table_1()}
        | _table_2()
        | _table_3()
        | _table_4()
        | _table_5()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

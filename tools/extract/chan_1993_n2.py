# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper, Sodhi & Brion (1993): N2 oscillator strengths, 11-200 eV.

The tables are unlabelled in the extraction except Table 4; the paper's numbers are
used here. Images checked against the cells (2026-09-27).

- Table 1 (``tbl_0c46ead68b63``, CSV): the low-resolution (1 eV FWHM) df/dE in
  ``10^-2 eV^-1`` at 11-200 eV, four column pairs read pair by pair (print order).
- Table 2 (``tbl_cffab8566e21``, cells file): the high-resolution discrete transitions,
  12.500-14.860 eV, in two column blocks (energy, final state, v', integrated
  oscillator strength). One file, rows and abscissa the transition energy (eV). "-"
  (no value) is skipped. Four values are printed once beside a brace over several
  assigned transitions (image): 13.345 eV (b 1Pi_u v' = 9 and o 1Pi_u v' = 1), 13.760
  and 13.785 eV (b' 1Sigma_u+ v' = 10 and b 1Pi_u v' = 14), 14.230 eV (b' v' = 16, c
  1Pi_u v' = 5, c' 1Sigma_u+ v' = 5) and 14.720 eV (c v' = 7 and c' v' = 7). The cells
  put the braced labels and the energy on different rows; the rows are joined here
  (``_BLENDS``) and each blend is one point at the first printed energy of the group,
  row ``<energy>+`` (``13.760+`` for the two-energy brace, at 13.760 eV).
- Tables 3-7 (``tbl_70338eb74bd4`` b 1Pi_u, ``tbl_9e453019ce36`` b' 1Sigma_u+,
  ``tbl_ce19b48ffb29`` c 1Pi_u, ``tbl_84a3b2fe8549`` c' 1Sigma_u+, ``tbl_359396e23484``
  o 1Pi_u; CSV): per-state band strengths, one file per state and column, rows and
  label v' (some rows print no energy). Footnote marks ``^{b)}``, ``^{d)}`` stripped;
  "-" skipped.
- Table 8 (``tbl_c43a8a42aa78``, cells file): each state's total, one file per
  reference (rows b, bp, c, cp, o for b, b', c, c', o); "(a) Hartree-Fock", "(b) CI"
  are variants of Bielschowsky et al.
- Table 9 (``tbl_21883bc6ed1a``, CSV): integrated oscillator strengths over energy
  intervals, 14.92-16.91 eV.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1993_AbsoluteOpticalOscillator_6b574e"
TABLE_1 = "tbl_0c46ead68b63"
TABLE_2 = "tbl_cffab8566e21"
TABLE_8 = "tbl_c43a8a42aa78"
TABLE_9 = "tbl_21883bc6ed1a"
STATE_TABLES = {
    "b_state": "tbl_70338eb74bd4",
    "bp_state": "tbl_9e453019ce36",
    "c_state": "tbl_ce19b48ffb29",
    "cp_state": "tbl_84a3b2fe8549",
    "o_state": "tbl_359396e23484",
}
PREFIX = "N2/chan_1993_n2__"
REFERENCES = {
    "1": "zipf_1978",
    "3": "ajello_1989",
    "4": "james_1990",
    "6": "rescigno_1978",
    "7": "kosman_1985",
    "8": "rose_1973",
    "9": "hazi_1981",
    "10": "bielschowsky_1990",
    "13": "stahel_1983",
    "27": "lawrence_1968",
    "28": "carter_1972",
    "48": "chutjian_1977",
}
TOTAL_STATES = ("b", "bp", "c", "cp", "o")
TOTAL_VARIANTS = {"(a) Hartree-Fock": "hf", "(b) CI": "ci"}
# Table 2's braced groups: (block, first row, last row) of the cells grid, with the
# row that carries the value; the image prints one value for the group.
_BLENDS = {
    (0, 16): (16, 17),  # 13.345: b v' = 9 and o v' = 1
    (0, 27): (27, 28),  # 13.760, 13.785: b' v' = 10 and b v' = 14
    (1, 11): (11, 13),  # 14.230: b' v' = 16, c v' = 5, c' v' = 5
    (1, 24): (24, 25),  # 14.720: c v' = 7 and c' v' = 7
}
_REFERENCE = re.compile(r"\[(\d+)\]")
_MARK = re.compile(r"\^\{[a-d]\)\}$")


def _stem(label: str) -> str:
    if label.lower().startswith("present"):
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


def _table_2() -> list[Point]:
    grid = json_table(DIRECTORY, TABLE_2)
    points = []
    for block in (0, 1):
        columns = slice(4 * block, 4 * block + 4)
        skip: set[int] = set()
        for index in range(1, len(grid)):
            if index in skip:
                continue
            energy, state, level, value = grid[index][columns]
            if (block, index) in _BLENDS:
                first, last = _BLENDS[block, index]
                group = [grid[i][columns] for i in range(first, last + 1)]
                skip.update(range(first + 1, last + 1))
                energies = [cells[0] for cells in group if cells[0]]
                printed = " || ".join(" | ".join(cells) for cells in group)
                row = energies[0] + ("+" if len(energies) > 1 else "")
                points.append(
                    Point(row, number(energies[0]), number(value), "", printed)
                )
                continue
            if value in {"", "-"}:
                continue  # no value printed ("-", or b v' = 14 inside a brace)
            printed = " | ".join((energy, state, level, value))
            points.append(Point(energy, number(energy), number(value), "", printed))
    return sorted(points, key=lambda p: float(p.x))


def _states() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for name, table_id in STATE_TABLES.items():
        rows = table(DIRECTORY, table_id)
        headers = rows[0]
        for cells in rows[1:]:
            level, printed = cells[1].strip(), " | ".join(cells)
            for column in range(2, len(cells)):
                cell = _MARK.sub("", cells[column].strip())
                if cell in {"", "-"}:
                    continue
                stem = _stem(headers[column])
                point = Point(level, level, number(cell), "", printed)
                files.setdefault(f"{PREFIX}{name}_{stem}.csv", []).append(point)
    return files


def _table_8() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    parent = ""
    for cells in json_table(DIRECTORY, TABLE_8)[1:]:
        if not any(cells[1:]):
            if _REFERENCE.search(cells[0]):
                parent = _stem(cells[0])  # Bielschowsky et al., with variant rows
            continue  # or the "Experiment" and "Theory" headings
        if cells[0] in TOTAL_VARIANTS:
            stem = f"{parent}_{TOTAL_VARIANTS[cells[0]]}"
        else:
            stem = _stem(cells[0])
        printed = " | ".join(cells)
        for state, cell in zip(TOTAL_STATES, cells[1:], strict=True):
            if cell:
                point = Point(state, state, number(cell), "", printed)
                files.setdefault(f"{PREFIX}totals_{stem}.csv", []).append(point)
    return files


def _table_9() -> list[Point]:
    points = []
    for cells in table(DIRECTORY, TABLE_9)[1:]:
        interval = cells[0].strip().replace("\N{EN DASH}", "-")
        low, high = interval.split("-")
        points.append(
            Point(
                interval,
                number(low),
                number(cells[1]),
                "",
                " | ".join(cells),
                {"x_upper": number(high)},
            )
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return Tables 1-9: df/dE, discrete, per-state, totals and intervals."""
    return (
        {
            f"{PREFIX}continuum.csv": _table_1(),
            f"{PREFIX}discrete_present.csv": _table_2(),
            f"{PREFIX}intervals_present.csv": _table_9(),
        }
        | _states()
        | _table_8()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

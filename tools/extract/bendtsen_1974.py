# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bendtsen (1974), J. Raman Spectrosc. 2, 133: the Raman spectrum of 14N2.

All tables checked against their images. Only 14N2 (the material N2) is taken; the
14N15N and 15N2 tables (II, III, V, VI) and rows are in the record's notes.

* Table I (``tbl_0a37c986866c``, CSV): the pure rotational Raman lines of 14N2, |nu|
  (the mean of the Stokes and anti-Stokes shifts) against the lower level J, with the
  estimated deviation of each line (the ``sigma`` column). The table prints J; the row
  is written as the S(J) line (J -> J + 2) it is.
* Table IV (``tbl_1caab6610c5d``, CSV): the lines of the fundamental band, one dataset
  per branch (O, Q, S), rows as printed ("O(17)").
* Table VII (``tbl_77805abff7bd``, JSON cells only): B_0, D x 10^6, nu_0 and B_0 - B_1
  from the pure rotational analysis ("Rot.") and the unconstrained and constrained
  rotation-vibration analyses, "+- n" in units of the last decimal. The constrained
  row's B_0 and D are printed in parentheses: the pure rotational values held fixed,
  not results, and not taken.
* Table VIII (``tbl_fb40cdaa8814``, JSON): alpha from the unconstrained and the
  constrained B_0 - B_1; the alpha mu^(3/2) columns are in the notes.
* Table IX (``tbl_b9751e81a126``, CSV): B_e (r_e in the notes: no quantity).
* Table X (the unlabelled ``tbl_66df62a78ce2``, JSON; its caption "TABLE X Molecular
  constants for 14N2" is extracted as text before it): the other authors' constants,
  one dataset per constant with a row per author. The values in parentheses
  ("transferred from Stoicheff [1] and used in a constrained analysis") are Stoicheff's,
  held fixed in Lofthus's and Wilkinson's analyses, and are not taken; nor is the "This
  work" row (Table VII's values, rounded). Stoicheff's nu_0 carries the note mark "a"
  ("Q branch maximum"), dropped from the number.

The uncertainties "+- n" are converted to the printed unit by the number of decimals of
the value ("1.989574 +- 12" is 0.000012); the printed cells are kept.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Bendtsen_1974_RotationalRotationVibrational"
TABLE_1 = "tbl_0a37c986866c"
TABLE_4 = "tbl_1caab6610c5d"
TABLE_7 = "tbl_77805abff7bd"
TABLE_8 = "tbl_fb40cdaa8814"
TABLE_9 = "tbl_b9751e81a126"
TABLE_10 = "tbl_66df62a78ce2"
PREFIX = "N2/bendtsen_1974__"
_PLUS_MINUS = re.compile(r"^\s*([\d.]+)\s*(?:\\pm|±)\s*(\d+)\s*(?:\^a)?\s*$")


def _last_decimal(value: str, digits: str) -> str:
    """Return an uncertainty printed in units of the value's last decimal."""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) * 10**-decimals:.{decimals}f}"


def _with_sigma(cell: str) -> tuple[str, str]:
    r"""Return the value and the absolute uncertainty of a "1.989574 \pm 12" cell."""
    match = _PLUS_MINUS.match(cell)
    if match is None:
        msg = f"not a 'value +- n' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1))
    return value, _last_decimal(value, match.group(2))


def _table1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_1)
    if grid[0] != ["J", "|v|_{obs}", "est. dev.", "obs. - calc."]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[1:]:
        label = f"S({int(cells[0])})"
        points.append(
            Point(label, label, number(cells[1]), number(cells[2]), " | ".join(cells))
        )
    return {f"{PREFIX}table1_S_N2.csv": points}


def _table4() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_4)
    if grid[0] != ["assign.", "v_{obs}", "est. dev.", "obs. - calc."]:
        msg = f"unexpected Table IV header {grid[0]!r}"
        raise ValueError(msg)
    branches: dict[str, list[Point]] = {"O": [], "Q": [], "S": []}
    for cells in grid[1:]:
        label = cells[0].strip()
        branches[label[0]].append(
            Point(label, label, number(cells[1]), number(cells[2]), " | ".join(cells))
        )
    return {
        f"{PREFIX}table4_{branch}_N2.csv": points for branch, points in branches.items()
    }


def _point(label: str, cell: str, printed: str) -> Point:
    value, sigma = _with_sigma(cell)
    return Point(label, label, value, sigma, printed)


def _table7() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_7)
    if grid[0][2:] != ["B_0", "D \\times 10^6", "\\nu_0", "B_0 - B_1"]:
        msg = f"unexpected Table VII header {grid[0]!r}"
        raise ValueError(msg)
    # Rows 1-3 are 14N2 (the material cell spans them); the image agrees.
    if not grid[1][0].startswith("^{14}\\mathrm {N}_2"):
        msg = f"unexpected Table VII row {grid[1]!r}"
        raise ValueError(msg)
    rot, free, fixed = grid[1], grid[2], grid[3]
    labels = [row[1] for row in (rot, free, fixed)]
    if labels != ["Rot.", "Rot.-Vib. unconstrained", "Rot.-Vib. constrained"]:
        msg = f"unexpected Table VII analyses {labels!r}"
        raise ValueError(msg)

    def printed(row: list[str]) -> str:
        return " | ".join(row)

    return {
        f"{PREFIX}table7_B0_rot_N2.csv": [_point(rot[1], rot[2], printed(rot))],
        f"{PREFIX}table7_D_rot_N2.csv": [_point(rot[1], rot[3], printed(rot))],
        f"{PREFIX}table7_B0_rotvib_N2.csv": [_point(free[1], free[2], printed(free))],
        f"{PREFIX}table7_D_rotvib_N2.csv": [_point(free[1], free[3], printed(free))],
        f"{PREFIX}table7_nu0_unconstrained_N2.csv": [
            _point(free[1], free[4], printed(free))
        ],
        f"{PREFIX}table7_nu0_constrained_N2.csv": [
            _point(fixed[1], fixed[4], printed(fixed))
        ],
        f"{PREFIX}table7_dB_unconstrained_N2.csv": [
            _point(free[1], free[5], printed(free))
        ],
        f"{PREFIX}table7_dB_constrained_N2.csv": [
            _point(fixed[1], fixed[5], printed(fixed))
        ],
    }


def _table8() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_8)
    if grid[0][1::2] != ["unconstrained", "constrained"] or grid[1][1::2] != [
        "\\alpha",
        "\\alpha",
    ]:
        msg = f"unexpected Table VIII header {grid[:2]!r}"
        raise ValueError(msg)
    row = grid[2]
    if row[0] != "^{14}\\mathrm {N}_{2}":
        msg = f"unexpected Table VIII row {row!r}"
        raise ValueError(msg)
    printed = " | ".join(row)
    return {
        f"{PREFIX}table8_alpha_unconstrained_N2.csv": [
            _point("unconstrained", row[1], printed)
        ],
        f"{PREFIX}table8_alpha_constrained_N2.csv": [
            _point("constrained", row[3], printed)
        ],
    }


def _table9() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_9)
    row = grid[1]
    if grid[0][1] != "B_e" or not row[0].startswith("^{14}\\text{N}_2"):
        msg = f"unexpected Table IX layout {grid[:2]!r}"
        raise ValueError(msg)
    return {f"{PREFIX}table9_Be_N2.csv": [_point("B_e", row[1], " | ".join(row))]}


# Table X: (dataset, column) and, per author row, the row label written to the file.
_TABLE_10_COLUMNS = {"B0": 2, "D": 3, "alpha": 4, "nu0": 5}
_TABLE_10_AUTHORS = {
    "Stoicheff [1]": "Stoicheff 1954",
    "Butcher et al. [2]": "Butcher 1971",
    "Lofthus [3]": "Lofthus 1956",
    "Wilkinson [4]": "Wilkinson 1957",
    "Vanderslice [5]": "Vanderslice 1964",
}


def _table10() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_10)
    if grid[0][2:] != ["B_0", "D \\times 10^6", "\\alpha", "v_0", "r_e(\\text{\\AA})"]:
        msg = f"unexpected Table X header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for name, column in _TABLE_10_COLUMNS.items():
        points = []
        for row in grid[1:]:
            if row[1] == "This work":
                continue  # Table VII's values, rounded
            label = _TABLE_10_AUTHORS[row[1]]
            cell = row[column].strip()
            if not cell or cell.startswith("("):
                continue  # none printed, or Stoicheff's value held fixed
            printed = " | ".join(row)
            if "pm" in cell or "±" in cell:
                points.append(_point(label, cell, printed))
            else:
                points.append(Point(label, label, number(cell), "", printed))
        files[f"{PREFIX}table10_{name}_N2.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the 14N2 files of Tables I, IV and VII-X."""
    return _table1() | _table4() | _table7() | _table8() | _table9() | _table10()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

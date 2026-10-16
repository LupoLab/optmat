# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sabbagh and Sadeghi (1977): transition probabilities of Xe I lines.

J. Quant. Spectrosc. Radiat. Transfer 17, 297-301 (1977) = PDF pages 1-5 (journal page
= PDF page + 296). Table 2 (PDF p. 4, journal p. 300; library table
``tbl_a841626f3146``, cells only, every cell checked against the table image
2026-10-09): 'Transition probabilities, 10^6 sec^-1', column 'Present work', 25 lines,
each with its Racah label and its air wavelength in Angstrom. Taken, the row named by
the transition in the paper's notation (``6p[5/2]3-6s[3/2]2``) and ``x`` the printed
wavelength:

* the 19 lines to the metastable 6s[3/2]2 measured against 8819 Angstrom (directly, or
  the blue and violet ones through 4671 Angstrom);
* the 4 lines to the resonant 6s[3/2]1 measured against 8280 Angstrom;
* the two reference lines, 8819 and 8280 Angstrom, whose values are the means of the
  literature values of Table 1 (p. 298), not measurements of this paper.

The library cell of the 8280 Angstrom row reads "6p[1/2]_6-6s[3/2]_1"; the image prints
6p[1/2]0 (a J = 6 level cannot decay to J = 1 by a dipole line), so the row label is
6p[1/2]0 and the cell is kept as extracted in ``printed``. The comparison columns
(Chen and Garstang; Aymar-Fehrenbach; Malakhov; Karimov and Klimkin) are not taken.

**Checks** (:func:`checks`; the script fails if one is missed): the abstract's counts
(19 lines to 6s[3/2]2 and four to 6s[3/2]1, besides the references); Table 2's
reference values equal Table 1's mean values ('30 ± 0.5', '34 ± 2.5',
``tbl_730574616f21``); and Gornik et al.'s f values (gornik_1981, Table I: 0.12, 0.15,
0.30) are the paper's Eq. (3) applied to its 8280, 8952.2 and 9923.2 Angstrom values,
to their two digits.
"""

import math
import re

from _json_table import json_table
from _library import Point, number, write

from optmat.constants import E_CHARGE, EPS0, M_E, C

DIRECTORY = "Sabbagh_1977_ExperimentalTransitionProbabilities"
SOURCE = "sabbagh_1977"
TABLE_1 = "tbl_730574616f21"
TABLE_2 = "tbl_a841626f3146"
REFERENCES = ("8819", "8280")
_LABEL = re.compile(r"^(\d[a-z]'?\[\d/2\])_(\d)-(6s\[3/2\])_(\d)$")
_VALUE = re.compile(r"^([\d.]+) ± ([\d.]+)$")
# Lower-level J and file of each group.
FILES = {"2": "Xe/sabbagh_1977__A_6s2_Xe.csv", "1": "Xe/sabbagh_1977__A_6s1_Xe.csv"}
REFERENCE_FILE = "Xe/sabbagh_1977__A_reference_Xe.csv"
# Gornik et al.'s f of 6s[3/2]1 -> 6p (their Table I), by Sabbagh's line; g_u / g_l.
GORNIK = {"8280": ("0.12", 1 / 3), "8952.2": ("0.15", 5 / 3), "9923.2": ("0.30", 5 / 3)}


def _rows() -> list[tuple[str, str, str, str, str]]:
    """Return (label, lower J, wavelength, value, sigma) of Table 2's 'Present work'."""
    grid = json_table(DIRECTORY, TABLE_2)
    if grid[1][2] != "Present work":
        msg = f"unexpected header {grid[1]!r}"
        raise ValueError(msg)
    rows = []
    for cells in grid[4:]:
        match = _LABEL.match(cells[0])
        value = _VALUE.match(cells[2])
        if match is None or value is None:
            msg = f"unexpected row {cells!r}"
            raise ValueError(msg)
        upper, j_upper, lower, j_lower = match.groups()
        if cells[1] == "8280" and j_upper == "6":
            j_upper = "0"  # the image prints 6p[1/2]0 (module notes)
        label = f"{upper}{j_upper}-{lower}{j_lower}"
        rows.append((label, j_lower, number(cells[1]), *value.groups()))
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the measured lines by lower level, and the two reference lines."""
    grid = json_table(DIRECTORY, TABLE_2)
    printed = {cells[1]: " | ".join(cells[:3]) for cells in grid[4:]}
    files: dict[str, list[Point]] = {
        name: [] for name in (*FILES.values(), REFERENCE_FILE)
    }
    for label, j_lower, wavelength, value, sigma in _rows():
        target = REFERENCE_FILE if wavelength in REFERENCES else FILES[j_lower]
        point = Point(label, wavelength, value, sigma, printed[wavelength])
        files[target].append(point)
    return files


def checks() -> None:
    """Check the counts, the references against Table 1 and Gornik et al.'s f."""
    files = extract()
    counts = {name: len(points) for name, points in files.items()}
    if [counts[FILES["2"]], counts[FILES["1"]], counts[REFERENCE_FILE]] != [19, 4, 2]:
        msg = f"counts {counts}"
        raise ValueError(msg)
    means = {cells[1]: cells[-1] for cells in json_table(DIRECTORY, TABLE_1)[2:]}
    for point in files[REFERENCE_FILE]:
        if means.get(point.x) != f"{point.value} ± {point.sigma}":
            msg = f"reference {point.x}: Table 1 {means.get(point.x)!r}"
            raise ValueError(msg)
    factor = EPS0 * M_E * C / (2 * math.pi * E_CHARGE**2)
    lines = {p.x: p for points in files.values() for p in points}
    for wavelength, (printed, weights) in GORNIK.items():
        a = float(lines[wavelength].value) * 1e6
        f = factor * (float(wavelength) * 1e-10) ** 2 * weights * a
        if f"{f:.2f}" != printed:
            msg = f"{wavelength}: f = {f:.4f} against Gornik et al.'s {printed}"
            raise ValueError(msg)


def main() -> None:
    """Check the table and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

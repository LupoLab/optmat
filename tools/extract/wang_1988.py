# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wang et al. (1988): the ionisation potentials of CO2 (He I photoelectrons).

Table 1 (``tbl_e1db285ab087``, image checked, p. 171): adiabatic ionisation potentials
(eV) of CO2, COS and CS2 to the X, A, B and C states. Only the CO2 row is taken (COS and
CS2 are not curated); its X cell carries the footnote mark "b" ("The uncertainties are
+-0.002 eV"), removed from the value and kept in ``printed``.

Table 3 (``tbl_25887e8a6993``, cells only; image checked, p. 173): the peak positions
(eV) of the X 2Pi_g band of CO2+, each assigned to a vibrational level (v1 v2 v3) of
the Omega = 3/2 or Omega = 1/2 spin-orbit component, or to a hot band. The digits in
parentheses are the precision of the relative positions (footnote a: 'The absolute
peak positions can only be determined within +-0.002 eV'), so they are kept in
``printed`` only, and the dataset carries the absolute +-0.002 eV. The first row's
footnote mark "^a" is removed from the value.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Wang_1988_HighResolutionUV"
TABLE_1 = "tbl_e1db285ab087"
TABLE_3 = "tbl_25887e8a6993"
STATES = {
    "\\tilde{X} ^{2}\\Pi_{\\text{g}}": "X 2Pi_g",
    "\\tilde{A} ^{2}\\Pi_{\\text{u}}": "A 2Pi_u",
    "\\tilde{B} ^{2}\\Sigma_{\\text{u}}^{+}": "B 2Sigma_u+",
    "\\tilde{C} ^{2}\\Sigma_{\\text{g}}^{+}": "C 2Sigma_g+",
}
_PEAK = re.compile(r"^([\d.]+)\((\d+)\)(?:\^a)?$")


def _table_1() -> list[Point]:
    grid = table(DIRECTORY, TABLE_1)
    states = [STATES[cell.strip()] for cell in grid[0][1:]]
    (cells,) = (row for row in grid[1:] if row[0].strip() == "CO_{2}")
    printed = " | ".join(cells)
    # The X cell's footnote mark "b" (the +-0.002 eV) is not part of the number.
    values = [number(cell.strip().removesuffix("b")) for cell in cells[1:]]
    return [
        Point(state.split()[0], state, value, "", printed)
        for state, value in zip(states, values, strict=True)
    ]


def _table_3() -> list[Point]:
    points = []
    for cells in json_table(DIRECTORY, TABLE_3)[2:]:
        match = _PEAK.match(cells[0].strip())
        if match is None:
            msg = f"unreadable Table 3 peak {cells[0]!r}"
            raise ValueError(msg)
        three_halves, one_half = cells[1].strip(), cells[2].strip()
        if bool(three_halves) == bool(one_half):
            msg = f"Table 3 row {cells!r} is not assigned to one component"
            raise ValueError(msg)
        if three_halves == "2_1^1(H.B.)":
            row, label = "hot_band", "2_1^1 (hot band), Omega = 3/2"
        else:
            omega = "3/2" if three_halves else "1/2"
            level = (three_halves or one_half).replace(" ", "")
            row, label = f"{level}_{omega}", f"({level}) Omega = {omega}"
        points.append(Point(row, label, number(match.group(1)), "", " | ".join(cells)))
    return points


def extract() -> dict[str, list[Point]]:
    """Return the adiabatic IPs of CO2 (Table 1) and the X-band peaks (Table 3)."""
    return {
        "CO2/wang_1988__adiabatic_CO2.csv": _table_1(),
        "CO2/wang_1988__x_band_CO2.csv": _table_3(),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

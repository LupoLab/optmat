# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract De et al. (2011): the calculated KERs of O2+ a 4Pi_u (Table II), TC-962.

Table II (p. 6, ``tbl_a5247d8a13ea``, read from its cells, ``.json``; checked against
the table image 2026-10-09) prints, for selected vibrational levels v of four states,
the "calculated oscillation frequency" nu_n (THz) and the kinetic energy release of
the fragments, the latter "calculated on the assumption of single-photon absorption"
of an 800 nm photon (caption). Only the O2+ a 4Pi_u block (v = 0-15, the levels from
Marian et al. 1982, their ref. [38]) is taken: its KER of O2+ -> O+ + O through the
f 4Pi_g state. The CO+ and CO++ blocks belong to materials optmat does not register;
the nu_n column has no quantity and stays in the ``printed`` cells.

Each row's KER is the level's energy plus one photon above the O+ + O limit, so the
KERs of adjacent levels v and v + 1 differ by their spacing. The script checks every
step against h times the mean of the two rows' printed nu_n, the frequency midway
between them on the anharmonic ladder (Eq. 2 of the paper is linear in v); the
largest misfit is 0.9 meV. Against h nu_v alone the steps are 0.9-2.4 meV short, so
the printed nu_n of row v is not the v -> v + 1 spacing of the KER column exactly.
"""

import itertools

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "De_2011_FollowingDynamicNuclear"
TABLE = "tbl_a5247d8a13ea"
HEADER = ["State", "\\nu", "\\nu_{\\text{n}} (THz)", "KER (eV)", "References"]
STATE = "O_{2}^{+}a^{4}\\Pi_{\\text{u}}"
LEVELS = 16
# Planck's constant in eV per THz (CODATA 2018, exact h and e).
H_EV_PER_THZ = 6.62607015e-34 / 1.602176634e-19 * 1e12
# The KERs are printed to 1 meV, so a difference of two carries up to 1 meV of
# rounding, and the frequencies to 0.01 THz (0.04 meV); 1.5 meV covers both and still
# catches a misread digit or a dropped minus sign (a step is 85-128 meV).
TOLERANCE_EV = 0.0015


def extract() -> dict[str, list[Point]]:
    """Return the O2+ a 4Pi_u block of Table II."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[0] != HEADER:
        msg = f"unexpected header of {TABLE}: {grid[0]!r}"
        raise ValueError(msg)
    if grid[1][0] != STATE or grid[1][4] != "[38]":
        msg = f"Table II does not open with O2+ a 4Pi_u [38]: {grid[1]!r}"
        raise ValueError(msg)
    block = grid[1 : 1 + LEVELS]
    if grid[1 + LEVELS][0] == "" or any(cells[0] for cells in block[1:]):
        msg = "the O2+ a 4Pi_u block is not 16 rows"
        raise ValueError(msg)
    points = []
    for level, cells in enumerate(block):
        if cells[1] != str(level):
            msg = f"row {level} of the O2+ block prints v = {cells[1]!r}"
            raise ValueError(msg)
        points.append(
            Point(cells[1], cells[1], number(cells[3]), "", " | ".join(cells))
        )
    rows = zip(points, block, strict=True)
    for (lower, low), (upper, high) in itertools.pairwise(rows):
        spacing = float(upper.value) - float(lower.value)
        mean = (float(number(low[2])) + float(number(high[2]))) / 2
        expected = mean * H_EV_PER_THZ
        if abs(spacing - expected) > TOLERANCE_EV:
            msg = f"v = {lower.row}: KER step {spacing:.4f} eV, h nu {expected:.4f} eV"
            raise ValueError(msg)
    return {"O2+/de_2011__ker_a4Piu_O2+.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

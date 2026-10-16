# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Balakrishnan, Smith and Stoicheff (1992): the dissociation energy of H2.

Phys. Rev. Lett. 68, 2149-2152 (1992); journal p. 2149 is PDF page 1. The Letter that
first reports the second dissociation limit of H2 by XUV fluorescence-excitation
spectroscopy near 84.5 nm with delayed quenching of H(2s); the full account, with D2
and a revised uncertainty budget, is ``balakrishnan_1994``.

Table I (``tbl_8f13575d8bbe``, p. 2151; the library has its cells only as JSON, read
with ``_json_table``; image checked) prints D0(H2) in cm^-1 by reference. Only the
'Present' row is taken (one file); the other rows are the values that
``balakrishnan_1994``'s Table III reproduces.

:func:`checks` tests the text's arithmetic (p. 2151, p. 2152): D_L = 118 258.57 +
118.49 = 118 377.06; D0 = D_L - 82 258.95 = 36 118.11, the value of the 'Present' row;
and D0(H2+) = D0 + I(H) - I(H2) with I(H2) the average of the two values the text
quotes, 21 379.36.
"""

import re
from typing import Final

from _json_table import json_table
from _library import Point, number, write

DIRECTORY: Final = "Balakrishnan_1992_DissociationEnergyHydrogen"
TABLE_I: Final = "tbl_8f13575d8bbe"
HEADER: Final = ["", "Reference", "D_0(H_2) (cm^{-1})"]
_CELL: Final = re.compile(r"^([\d .]+?)\s*(?:±\s*([\d.]+))?$")

DL_J1: Final = 118258.57
"""The threshold above J'' = 1, +-0.04 (p. 2151)."""
F1: Final = 118.49
"""The rotational spacing between J'' = 0 and 1 (p. 2151, their Ref. [9])."""
DL: Final = 118377.06
"""The second dissociation limit (abstract; p. 2151)."""
INTERVAL: Final = 82258.95
"""The 2s - 1s separation of H (p. 2151, their Ref. [8])."""
I_H: Final = 109678.764
"""IP(H) (p. 2152, their Ref. [10])."""
I_H2: Final = (124417.501, 124417.524)
"""The two IP(H2) values averaged (p. 2152, their Refs. [13, 14])."""
D0_ION: Final = 21379.36
"""D0(H2+) (abstract; p. 2152)."""


def _present() -> tuple[str, str, list[str]]:
    """Return Table I's 'Present' value, uncertainty and row."""
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    rows = [row for row in grid if row[1] == "Present"]
    if len(rows) != 1:
        msg = f"Table I: {len(rows)} 'Present' rows"
        raise ValueError(msg)
    match = _CELL.match(rows[0][2])
    if match is None or match.group(2) is None:
        msg = f"Table I: cannot read {rows[0][2]!r}"
        raise ValueError(msg)
    return number(match.group(1).replace(" ", "")), number(match.group(2)), rows[0]


def extract() -> dict[str, list[Point]]:
    """Return the Letter's D0(H2), one file."""
    value, sigma, row = _present()
    point = Point("Present", "H2", value, sigma, " | ".join(row))
    return {"H2/balakrishnan_1992__d0_H2.csv": [point]}


def checks() -> None:
    """Check the Letter's arithmetic against Table I and the text; raise if not."""
    d0 = float(_present()[0])
    for name, got, printed in (
        ("D_L", DL_J1 + F1, DL),
        ("D0(H2)", DL - INTERVAL, d0),
        ("D0(H2+)", d0 + I_H - sum(I_H2) / 2, D0_ION),
    ):
        if abs(got - printed) > 0.005:
            msg = f"{name}: recomputed {got:.4f}, printed {printed:.4f}"
            raise ValueError(msg)


def main() -> None:
    """Check the Letter, then write the data file."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

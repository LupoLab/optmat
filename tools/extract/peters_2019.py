# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Peters et al. (2019): the N2+ dissociative recombination coefficient.

Table I (``tbl_7c2f4f3bf473``, cells only; image checked, p. 243301-3): the model's
reactions and rate coefficients. Of the dissociative recombinations, only R12, e- +
N2+ -> 2N, concerns a registered ion; it is printed "2.8 x 10^-7 sqrt(300/T_e)" cm^3
s^-1, refs. 23 (Kossyi et al. 1992) and 40, not tuned. The coefficient's mantissa is
written (the record's scale holds 10^-7); the printed form is checked here.

Register: TC-977.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Peters_2019_KineticsModelFemtosecond"
TABLE_I = "tbl_7c2f4f3bf473"
_R12 = re.compile(
    r"^([\d.]+) \\times 10\^\{-7\} \\sqrt\{\\frac\{300\}\{T_\{\\mathrm\{e\}\}\}\}$"
)


def extract() -> dict[str, list[Point]]:
    """Return R12's coefficient at 300 K."""
    for cells in json_table(DIRECTORY, TABLE_I):
        if cells[0].startswith("R12"):
            match = _R12.match(cells[2])
            if match is None:
                msg = f"R12 printed as {cells[2]!r}"
                raise ValueError(msg)
            point = Point("R12", "", match.group(1), "", " | ".join(cells))
            return {"N2+/peters_2019__dr_N2+.csv": [point]}
    msg = "no row R12 in Table I"
    raise ValueError(msg)


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

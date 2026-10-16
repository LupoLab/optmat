# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mahla & Mehnen (2025): vertical ionisation potentials of O2 and O3.

Table I (``tbl_0127e8228e7c``, cells only; image checked, p. 5): per target and final
ionic state, the R-matrix value (this work) and the experimental values, each marked by
its footnote (a: PES, ref. 74; b: PFI, ref. 75; c: PES, ref. 76; d: PES, refs. 35, 41,
42; e: PES, refs. 28-30; f: PFI-ZEKE, ref. 77). One file per target and column or
footnote, rows labelled by the ionic state (in ASCII). "12.75^{e,f}" belongs to both e
and f.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Mahla_2025_PhotoionizationDynamicsO2"
TABLE = "tbl_0127e8228e7c"
TARGETS = {"O_{2}": "O2", "O_{3}": "O3"}
STATES = {
    "X^{2}\\Pi_{g}": "X 2Pi_g",
    "a^{4}\\Pi_{u}": "a 4Pi_u",
    "A^{2}\\Pi_{u}": "A 2Pi_u",
    "b^{4}\\Sigma_{g}^{-}": "b 4Sigma_g-",
    "B^{2}\\Sigma_{g}^{-}": "B 2Sigma_g-",
    "3^{2}\\Pi_{u}": "3 2Pi_u",
    "^{2}A_{1}(6a_{1}^{-1})": "2A1 (6a1^-1)",
    "^{2}B_{2}(4b_{2}^{-1})": "2B2 (4b2^-1)",
    "^{2}A_{2}(1a_{2}^{-1})": "2A2 (1a2^-1)",
}
MARKS = {
    "a": "baltzer",
    "b": "tonkyn",
    "c": "edqvist",
    "d": "katsumata",
    "e": "brundle",
    "f": "willitsch",
}
_ENTRY = re.compile(r"^([\d.]+)\^\{([a-f](?:,[a-f])*)\}$")


def extract() -> dict[str, list[Point]]:
    """Return the R-matrix and experimental vertical I_p of O2 and O3."""
    files: dict[str, list[Point]] = {}
    target = ""
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        target = TARGETS.get(cells[0].strip(), target)
        state = STATES[cells[1].strip()]
        printed = " | ".join(cells)
        files.setdefault(f"{target}/mahla_2025__rmatrix_{target}.csv", []).append(
            Point(state, state, number(cells[2]), "", printed)
        )
        for entry in (e.strip() for e in cells[3].split(", ")):
            match = _ENTRY.match(entry)
            if match is None:
                msg = f"unreadable experimental entry {entry!r}"
                raise ValueError(msg)
            for mark in match.group(2).split(","):
                path = f"{target}/mahla_2025__expt_{MARKS[mark]}_{target}.csv"
                files.setdefault(path, []).append(
                    Point(state, state, number(match.group(1)), "", printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

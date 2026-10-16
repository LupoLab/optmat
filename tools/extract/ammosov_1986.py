# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ammosov, Delone & Krainov (1986): the coefficients C_n*l* of Table II.

Table II (``tbl_06845899aeaa``, cells file only; image checked, p. 1193) is printed in
two halves side by side, each with the columns atom or ion, n*, "A(Z/n*)^(n*-0.5)
(A-from handbook^4)" and "From Eq. (11)". Per atom or ion, the handbook's numerical
coefficient converted to Eq. (3)'s definition, with its "+-" error
(``handbook_<target>``), and the analytic value of Eq. (11) (``eq11_<target>``); n* is
written as the per-row column ``n_star``. The page prints decimal commas
("1,99+-0,007"); the cells file has points except in the Kr+ cell "2,30=0.12", whose
"=" stands where every other cell of the column prints "+-" (the page prints "=" too,
a print defect): it is read as 2.30 +- 0.12. Every target of the table is a curated
atom or ion.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Ammosov_1986_TunnelIonizationComplex"
TABLE_II = "tbl_06845899aeaa"
TARGETS = {
    "He": "He",
    "Ne": "Ne",
    "Ne^+": "Ne+",
    "Ar": "Ar",
    "Ar^+": "Ar+",
    "Kr": "Kr",
    "Kr^+": "Kr+",
    "Xe": "Xe",
    "Xe^+": "Xe+",
}
_HANDBOOK = re.compile(r"^(?P<value>\d[.,]\d+)\s*[±=]\s*(?P<error>\d[.,]\d+)$")


def _decimal(text: str) -> str:
    """Return a printed number with its decimal comma as a point."""
    return number(text.replace(",", "."))


def extract() -> dict[str, list[Point]]:
    """Return the handbook and Eq. (11) coefficients of each atom and ion."""
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_II)
    for cells in grid[2:]:
        for half in (cells[:4], cells[4:]):
            label = half[0].strip()
            if not label:
                continue
            target = TARGETS[label]
            printed = " | ".join(half)
            n_star = {"n_star": _decimal(half[1])}
            match = _HANDBOOK.match(half[2].strip())
            if match is None:
                msg = f"unreadable handbook cell {half[2]!r}"
                raise ValueError(msg)
            files[f"{target}/ammosov_1986__handbook_{target}.csv"] = [
                Point(
                    "C",
                    "",
                    _decimal(match["value"]),
                    _decimal(match["error"]),
                    printed,
                    n_star,
                )
            ]
            files[f"{target}/ammosov_1986__eq11_{target}.csv"] = [
                Point("C", "", _decimal(half[3]), "", printed, n_star)
            ]
    if len(files) != 2 * len(TARGETS):
        msg = f"Table II lost a row: {sorted(files)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

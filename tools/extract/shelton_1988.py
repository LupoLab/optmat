# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton & Ulivi (1988): the vibrational hyperpolarizability of SF6.

Table III (``tbl_5ed5e11e4508``, JSON cells; image checked) gives gamma^v in 1e-62
C^4 m^4 J^-3 for the static limit, dc Kerr, ESHG and THG at nu = 20000 cm^-1, split into
the mu^4, mu^2 alpha + mu alpha mu, mu beta and alpha^2 terms, and their total. One file
per process for the total (``gammav_<process>``) and one for the four terms
(``gammav_terms_<process>``). Some cells carry their own power of ten ("-1.9 \times
10^{-3}"): unlike a column factor it differs from cell to cell, so it is kept in the
number as E notation ("-1.9e-3"), not stripped.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Shelton_1988_VibrationalHyperpolarizabilitySF6"
TABLE = "tbl_5ed5e11e4508"
PROCESSES = {1: ("static", "inf"), 2: ("dc_kerr", "20000"), 3: ("eshg", "20000")}
PROCESSES[4] = ("thg", "20000")
TERMS = {
    "\\mu^{4}": "mu^4",
    "\\mu^{2}\\alpha + \\mu\\alpha\\mu": "mu^2 alpha + mu alpha mu",
    "\\mu\\beta": "mu beta",
    "\\alpha^{2}": "alpha^2",
}
_POWER = re.compile(r"^\s*(-?[\d.]+)\s*\\times\s*10\^\{(-?\d+)\}\s*$")


def value(cell: str) -> str:
    """Return the printed number, with a per-cell power of ten as E notation."""
    match = _POWER.match(cell)
    if match:
        return f"{match.group(1)}e{match.group(2)}"
    return number(cell)


def extract() -> dict[str, list[Point]]:
    """Return one total and one term file per process."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[1][1:] != ["Static", "dc Kerr", "ESHG", "THG"]:
        msg = f"Table III columns changed: {grid[1]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for column, (process, wavenumber) in PROCESSES.items():
        terms = []
        for cells in grid[2:]:
            printed = " | ".join(cells)
            if cells[0] == "Total":
                files[f"SF6/shelton_1988__gammav_{process}.csv"] = [
                    Point("total", wavenumber, value(cells[column]), "", printed)
                ]
            else:
                label = TERMS[cells[0]]
                terms.append(
                    Point(label, wavenumber, value(cells[column]), "", printed)
                )
        files[f"SF6/shelton_1988__gammav_terms_{process}.csv"] = terms
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

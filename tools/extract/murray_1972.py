# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Murray & Javan (1972), J. Mol. Spectrosc. 42, 1: H2, D2 Q widths and shifts.

Both tables exist only as JSON cells; both were checked against their images. Values are
in 10^-3 cm-1 per amagat (per amagat^2 for b) at 25 C; cells print "v \pm s" or "v".

* Table II (``tbl_cab7e16264ae``, p. 11): density broadening coefficients of the
  Q_1(n) lines, "This work" with its rms deviation delta (kept in ``printed``), May et
  al. (Ref. 11) and the theory of Ref. 32, for H2 (n = 0-3) and D2 (n = 1-3; May et al.
  and the theory give none). The gas names head their sections in rows of their own;
  delta's header sits in a row between the H2 header and the first line.
* Table III (``tbl_ca64e120ea87``, p. 11): the line-shift coefficients, this work's
  linear fit (a, delta), its quadratic fit (a, b, delta; H2 Q_1(0) and Q_1(1) only) and
  May et al.'s quadratic fit (a, b; H2 n = 0-2). **Structural fix:** the library merged
  each gas's section header into its first line's cell ("HydrogenQ_{1}(n)n = 0",
  "DeuteriumQ_{1}(n)n = 1"); the image prints the header on lines of its own. The
  script reads the gas and n from that cell.

Rows are the lines, "Q1(n)". One file per table, column group and gas.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Murray_1972_EffectsCollisionsRaman"
TABLE_II = "tbl_cab7e16264ae"
TABLE_III = "tbl_ca64e120ea87"
GASES = {"Hydrogen": "H2", "Deuterium": "D2"}
_LINE = re.compile(r"^(Hydrogen|Deuterium)?(?:Q_\{1\}\(n\))?(?:n = )?(\d)$")


def _split(cell: str) -> tuple[str, str]:
    r"""Return (value, sigma) of a "v \pm s" or "v" cell."""
    value, _, sigma = cell.partition("\\pm")
    return number(value), number(sigma) if sigma else ""


def _rows(grid: list[list[str]]) -> list[tuple[str, str, list[str]]]:
    """Return (gas, line label, cells) for each data row of a table."""
    rows = []
    gas = ""
    for cells in grid:
        header = cells[0].split("Q_{1}(n)")[0].strip()
        if header in GASES and not any(cells[1:]):
            gas = GASES[header]  # a section header on its own row
            continue
        match = _LINE.match(cells[0])
        if match is None:
            continue  # a column-header row
        if match.group(1):
            gas = GASES[match.group(1)]  # the header merged into the row (Table III)
        if not gas:
            msg = f"a line before any gas: {cells!r}"
            raise ValueError(msg)
        rows.append((gas, f"Q1({match.group(2)})", cells))
    return rows


def _add(
    files: dict[str, list[Point]],
    name: str,
    gas: str,
    line: str,
    cell: str,
    printed: str,
) -> None:
    if not cell:
        return
    value, sigma = _split(cell)
    files.setdefault(f"{gas}/murray_1972__{name}_{gas}.csv", []).append(
        Point(line, line, value, sigma, printed)
    )


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    if grid[0] != ["Line", "This work", "", "May, et al. (11)", "Theoretical (32)"]:
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for gas, line, cells in _rows(grid[1:]):
        printed = " | ".join(cells)
        _add(files, "table2", gas, line, cells[1], printed)
        _add(files, "table2_may", gas, line, cells[3], printed)
        _add(files, "table2_theory", gas, line, cells[4], printed)
    return files


def _table_iii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    if grid[1] != ["", "a", "\\delta", "a", "b", "\\delta", "a", "b"]:
        msg = f"unexpected Table III header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    columns = {
        1: "table3_linear",
        3: "table3_quadratic_a",
        4: "table3_quadratic_b",
        6: "table3_may_a",
        7: "table3_may_b",
    }
    for gas, line, cells in _rows(grid[2:]):
        printed = " | ".join(cells)
        for index, name in columns.items():
            _add(files, name, gas, line, cells[index], printed)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Tables II and III."""
    return _table_ii() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

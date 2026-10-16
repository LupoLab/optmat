# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rouille et al. (1992), J. Mol. Spectrosc. 154, 372: O2 lines, constants.

All tables were checked against their images.

* Table I (``tbl_02cb213c6e14``, JSON, p. 376): the v = 0 constants of this work and of
  Zink and Mizushima (Ref. 5), "v(u) \times 10^{-k}" cells (u: one standard deviation
  on the last digits). One file per source and constant; the row is v.
* Table II (``tbl_be0dee4f0c09``, CSV, p. 377): the ^SS(N, J) and ^OO(N, J) triplet
  lines of the 1-0 band (the Obs. - calc. columns, 10^-4 cm-1, stay in ``printed``).
  Rows "SS(N,J)", "OO(N,J)" with N and J of the lower level.
* Table III (``tbl_b87d755b4ebf``, CSV, p. 378): the v = 1 constants of this work and of
  Refs. 6, 7 and 16, and the band origin nu_10. **Structural fix:** the library labels
  the primed v = 1 rows "\lambda_i", "\lambda_i''", "\mu_i'", "\mu_i''", "\mu_i''''";
  the image prints lambda_1', lambda_1'', mu_1', mu_1'', mu_1''' (the subscript 1 read
  as i, one prime lost on lambda_1'). The script maps the row labels to the image's.
  Values marked a (fixed to v = 0) or b (fixed) are not written (the record's notes).
* Table IV (``tbl_592310ecda73``, CSV, p. 379): the calculated Q-branch triplet
  components, this work and Dang-Nhu et al. (Ref. 14); the observed column repeats
  Table V.
* Table V (``tbl_9c9d20888448``, CSV, p. 381): the unresolved Q(N) lines of the
  fundamental and the first and second hot bands, with their precision. The column
  header prints J; the lines are labelled by N (Table IV's Q(N); J = N when the spin is
  ignored, Eq. 11).
* Table VI (``tbl_796b518eb4e9``, CSV, p. 381): the v = 2 and 3 constants; the band
  origins and B_v (note a, fitted) and the extrapolated D, lambda, lambda', mu, mu'
  are written; the fixed H, lambda'', mu'', mu''' are not.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Rouille_1992_HighResolutionStimulated"
TABLE_I = "tbl_02cb213c6e14"
TABLE_II = "tbl_be0dee4f0c09"
TABLE_III = "tbl_b87d755b4ebf"
TABLE_IV = "tbl_592310ecda73"
TABLE_V = "tbl_9c9d20888448"
TABLE_VI = "tbl_796b518eb4e9"
_CELL = re.compile(
    r"^(?P<value>-?\d+\.?\d*)(?:\((?P<digits>\d+)\))?"
    r"(?:\\times ?10\^?\{?-(?P<power>\d+)\}?)?(?P<note>[ab]?)$"
)
# Row labels of Tables I, III and VI (as the library reads them) and the file stems.
CONSTANTS = {
    "B": "B",
    "D": "D",
    "H": "H",
    "\\lambda": "lambda",
    "\\lambda'": "lambda_D",
    "\\lambda''": "lambda_H",
    "\\mu": "mu",
    "\\mu'": "mu_D",
    "\\mu''": "mu_H",
    "\\mu'''": "mu_3",
}
# Table III's misread labels of the v = 1 primed constants (see the docstring).
TABLE_III_LABELS = {
    "\\lambda_i": "\\lambda'",
    "\\lambda_i''": "\\lambda''",
    "\\mu_i'": "\\mu'",
    "\\mu_i''": "\\mu''",
    "\\mu_i''''": "\\mu'''",
}
# Table VI's H, lambda'', mu'' and mu''' are the v = 0 (or zero) values, fixed without
# a mark (p. 379): notes only. Tables I and III mark their fixed values (a, b).
FIXED_VI = {"H", "lambda_H", "mu_H", "mu_3"}


def _parse(cell: str) -> tuple[str, str, str, str]:
    """Return (value, sigma, power of ten, note) of a constant's cell."""
    match = _CELL.match(cell.replace(" ", ""))
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value = number(match["value"])
    sigma = ""
    if match["digits"]:
        decimals = len(value.partition(".")[2])
        sigma = f"{int(match['digits']) / 10**decimals:.{decimals}f}"
    return value, sigma, match["power"] or "0", match["note"]


def _label(cell: str) -> tuple[str, str]:
    r"""Split a row label such as "\lambda_0'" into ("\lambda'", "0")."""
    match = re.fullmatch(r"(\\?[A-Za-z]+)_?(\d|v)?('*)", cell)
    if match is None:
        msg = f"unexpected row label {cell!r}"
        raise ValueError(msg)
    return match[1] + match[3], match[2] or ""


def _constants(
    grid: list[list[str]],
    stem: str,
    sources: dict[int, str],
    labels: dict[str, str],
    rows_v: list[str],
) -> dict[str, list[Point]]:
    """Read a table of constants: one file per source and constant, rows v."""
    files: dict[str, list[Point]] = {}
    powers: dict[str, set[str]] = {}
    for cells in grid[1:]:
        raw = labels.get(cells[0], cells[0])
        if raw.startswith("^") or raw.startswith("v_") or raw.startswith("\\nu"):
            continue  # a note line, or a band origin (read elsewhere)
        symbol, _ = _label(raw)
        name = CONSTANTS[symbol]
        for column, source in sources.items():
            cell = cells[column]
            if not cell:
                continue
            value, sigma, power, note = _parse(cell)
            # In Tables I and III, notes a and b mark fixed values; in Table VI, note a
            # marks the fitted ones and the fixed ones are unmarked (FIXED_VI).
            fixed = name in FIXED_VI if stem == "table6" else bool(note)
            if fixed or name == "mu_3":
                continue  # fixed values (and mu''' = 0): the record's notes
            v = rows_v[0] if len(rows_v) == 1 else rows_v[column - 1]
            path = f"O2/rouille_1992__{stem}_{source}_{name}.csv"
            powers.setdefault(path, set()).add(power)
            files.setdefault(path, []).append(
                Point(v, v, value, sigma, f"{cells[0]} | {cell}")
            )
    # One column factor per file: the record's scale is 10^-power.
    mixed = {path: found for path, found in powers.items() if len(found) > 1}
    if mixed:
        msg = f"column factors differ within a file: {mixed!r}"
        raise ValueError(msg)
    return files


def _table_i() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0] != ["", "This work", "(5)"]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    return _constants(grid, "table1", {1: "this_work", 2: "zink"}, {}, ["0"])


def _table_iii() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_III)
    if grid[0] != ["", "This work", "(6)", "(7)", "(16)"]:
        msg = f"unexpected Table III header {grid[0]!r}"
        raise ValueError(msg)
    sources = {1: "this_work", 2: "endo", 3: "loete", 4: "albritton"}
    files = _constants(grid, "table3", sources, TABLE_III_LABELS, ["1"])
    origin = grid[1]
    if origin[0] != "v_{10}":
        msg = f"unexpected Table III origin row {origin!r}"
        raise ValueError(msg)
    for column, source in sources.items():
        if origin[column]:
            value, sigma, _, _ = _parse(origin[column])
            files[f"O2/rouille_1992__table3_{source}_nu.csv"] = [
                Point("1-0", "1-0", value, sigma, " | ".join(origin))
            ]
    return files


def _table_vi() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_VI)
    if grid[0] != ["", "v = 2", "v = 3"]:
        msg = f"unexpected Table VI header {grid[0]!r}"
        raise ValueError(msg)
    files = _constants(grid, "table6", {1: "this_work", 2: "this_work"}, {}, ["2", "3"])
    origin = grid[1]
    if origin[0] != "v_{vv-1}":
        msg = f"unexpected Table VI origin row {origin!r}"
        raise ValueError(msg)
    files["O2/rouille_1992__table6_this_work_nu.csv"] = [
        Point(band, band, *_parse(cell)[:2], " | ".join(origin))
        for band, cell in (("2-1", origin[1]), ("3-2", origin[2]))
    ]
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_II)
    if grid[0][:3] != ["N", "J", "^{S}S(N, J)"] or grid[0][4] != "^{O}O(N, J)":
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        n, j = cells[0], cells[1]
        printed = " | ".join(cells)
        for column, branch in ((2, "SS"), (4, "OO")):
            if cells[column]:
                row = f"{branch}({n},{j})"
                files.setdefault(f"O2/rouille_1992__table2_{branch}.csv", []).append(
                    Point(row, row, number(cells[column]), "", printed)
                )
    return files


def _table_iv() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_IV)
    if grid[0] != ["N", "J", "(14) (calc.)", "This work (calc.)", "Q(N) (obs.)"]:
        msg = f"unexpected Table IV header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        row = f"QQ({cells[0]},{cells[1]})"
        printed = " | ".join(cells)
        for column, name in ((2, "table4_dang_nhu"), (3, "table4_this_work")):
            if cells[column]:
                files.setdefault(f"O2/rouille_1992__{name}.csv", []).append(
                    Point(row, row, number(cells[column]), "", printed)
                )
    return files


def _table_v() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_V)
    if grid[0] != [
        "J",
        "Fundamental band",
        "J",
        "First hot-band",
        "J",
        "Second hot-band",
    ]:
        msg = f"unexpected Table V header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        printed = " | ".join(cells)
        for column, band in ((0, "1-0"), (2, "2-1"), (4, "3-2")):
            if not cells[column]:
                continue
            row = f"Q({cells[column]})"
            value, sigma, _, _ = _parse(cells[column + 1])
            name = f"O2/rouille_1992__table5_Q_{band.replace('-', '')}.csv"
            files.setdefault(name, []).append(Point(row, row, value, sigma, printed))
    return files


def extract() -> dict[str, list[Point]]:
    """Return the O2 files of Tables I-VI."""
    return (
        _table_i() | _table_ii() | _table_iii() | _table_iv() | _table_v() | _table_vi()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

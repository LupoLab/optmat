# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Smyth, Rosasco & Hurst (1987), J. Chem. Phys. 87, 1001: D2 Q-branch widths.

Both tables have CSVs, checked against their images:

* Table I (``tbl_bedfcf2a5c4b``, p. 1005): the HWHM line-broadening coefficients
  (10^-3 cm-1/amagat) of D2 Q(0)-Q(5) by D2, He, Ar, H2 and CH4; one file per
  collision partner, one row per line.
* Table III (``tbl_01b17bbca40c``, p. 1008): pure D2 against J, this work (Table I's
  first row with its errors), the paper's theory column, and four other measurements;
  one file per column. Cells print "v ± s"; blank cells are lines a source did not give.

Table II (the rate-law parameters) is recorded as formulas in the source record.
"""

from _library import Point, number, table, write

DIRECTORY = "Smyth_1987_MeasurementRateLaw"
TABLE_I = "tbl_bedfcf2a5c4b"
TABLE_III = "tbl_01b17bbca40c"
PERTURBERS = {
    "D_{2}:D_{2}": "self",
    "D_{2}:He": "He",
    "D_{2}:Ar": "Ar",
    "D_{2}:H_{2}": "H2",
    "D_{2}:CH_{4}": "CH4",
}
TABLE_III_COLUMNS = {
    1: "table3_this_work",
    2: "table3_theory",
    3: "table3_murray_javan",
    4: "table3_robertson",
    5: "table3_gustafson_byer",
    6: "table3_roh",
}


def _split(cell: str) -> tuple[str, str]:
    """Return (value, sigma) of a "v ± s" or "v" cell."""
    value, _, sigma = cell.replace(" ", "").partition("±")
    return number(value), number(sigma) if sigma else ""


def _table_i() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != ["Gas mixture", "Q(0)", "Q(1)", "Q(2)", "Q(3)", "Q(4)", "Q(5)"]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        name = f"D2/smyth_1987__table1_{PERTURBERS[cells[0]]}.csv"
        files[name] = [
            Point(line, line, number(cell), "", " | ".join(cells))
            for line, cell in zip(grid[0][1:], cells[1:], strict=True)
        ]
    return files


def _table_iii() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_III)
    header = [
        "Rotational state J",
        "This work",
        "Theory^a",
        "Murray and Javan (Ref. 10)",
        "Robertson et al. (Ref. 31)",
        "Gustafson and Byer (Ref. 32)^b",
        "Roh (Ref. 33)",
    ]
    if grid[0] != header:
        msg = f"unexpected Table III header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        line = f"Q({int(cells[0])})"
        for column, name in TABLE_III_COLUMNS.items():
            if not cells[column]:  # no value from that source
                continue
            value, sigma = _split(cells[column])
            files.setdefault(f"D2/smyth_1987__{name}.csv", []).append(
                Point(line, line, value, sigma, " | ".join(cells))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I and Table III files."""
    return _table_i() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

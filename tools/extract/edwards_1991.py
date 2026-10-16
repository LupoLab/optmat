# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Edwards, Long & Sherwood (1991): H2 and D2 pure rotational Raman broadening.

J. Raman Spectrosc. 22, 607. Broadening coefficients 10^3 sigma in cm-1 bar-1 of the
S_0(J) lines, from the cells files (Table 1 also as CSV: a summary without values, in
the record's notes); every table checked against its image:

* Table 2 (``tbl_10ca63d1c797``): self-broadening of 1H2 and 2H2 at 293 K. The D2
  S_0(4) cell prints a bullet where "±" is meant ("1.14 ● 0.10": the image shows the
  bullet; every other cell has "±"); it is read as "±".
* Table 3 (``tbl_9facf32ff8c9``): self-broadening of 1H2 at 90 K and 2H2 at 90, 200
  and 375 K (the merged header cells run the molecule into the temperature,
  "^1H_290 K").
* Tables 4 and 5 (``tbl_e5739c8fc9bb``, ``tbl_5006f1b3c92e``): 1H2 and 2H2 broadened
  by He, Ar, O2, N2, CO and HCl at 293 K ("±" printed as ``\pm``). HCl (not a curated
  material) is left in the record's notes.
* Table 6 (``tbl_c97b652acc1c``): 2H2 S_0(0) broadened by N2 and CO, 100-293 K.
* Tables 7-9 (``tbl_2a3deb14c621``, ``tbl_3b59e5b50be9``, ``tbl_27adee9b3c93``): the
  comparisons with earlier work (no uncertainties printed).

"—" is no value.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Edwards_1991_LineWidthsPure"
TABLE_2 = "tbl_10ca63d1c797"
TABLE_3 = "tbl_9facf32ff8c9"
TABLE_4 = "tbl_e5739c8fc9bb"
TABLE_5 = "tbl_5006f1b3c92e"
TABLE_6 = "tbl_c97b652acc1c"
TABLE_7 = "tbl_2a3deb14c621"
TABLE_8 = "tbl_3b59e5b50be9"
TABLE_9 = "tbl_27adee9b3c93"
DASH = "\N{EM DASH}"
BULLET = "\N{BLACK CIRCLE}"
PERTURBERS = {
    "He": "He",
    "Ar": "Ar",
    "O_2": "O2",
    "O_{2}": "O2",
    "N_2": "N2",
    "N_{2}": "N2",
    "CO": "CO",
}
# HCl ("HCI" in Table 4's cells) is not a curated material: its rows stay in the notes.
SKIPPED = {"HCI", "HCl"}


def _split(cell: str) -> tuple[str, str]:
    r"""Return (value, sigma) of a "v ± s", "v \pm s" or "v" cell."""
    text = cell.replace("\\pm", "±").replace(BULLET, "±").replace(" ", "")
    value, _, sigma = text.partition("±")
    return number(value), number(sigma) if sigma else ""


def _line(cell: str) -> str:
    """Return the S_0(J) row name of a "0" or "S_0(0)" cell."""
    j = cell.removeprefix("S_0(").removesuffix(")")
    return f"S0({number(j)})"


def _check(grid: list[list[str]], cells: dict[tuple[int, int], str]) -> None:
    """Fail loudly if header cells are not where the image puts them."""
    for (row, column), text in cells.items():
        if grid[row][column] != text:
            msg = f"unexpected header cell ({row}, {column}) {grid[row][column]!r}"
            raise ValueError(msg)


def _by_line(
    grid: list[list[str]], first: int, columns: dict[int, str]
) -> dict[str, list[Point]]:
    """Return one file per column of S_0(J) rows."""
    files: dict[str, list[Point]] = {}
    for cells in grid[first:]:
        row = _line(cells[0])
        for column, name in columns.items():
            if cells[column] == DASH:
                continue
            value, sigma = _split(cells[column])
            files.setdefault(name, []).append(
                Point(row, row, value, sigma, " | ".join(cells))
            )
    return files


def _by_perturber(
    grid: list[list[str]], material: str, table: str
) -> dict[str, list[Point]]:
    """Return one file per perturber row of S_0(J) columns."""
    files: dict[str, list[Point]] = {}
    lines = [_line(cell.replace("{", "").replace("}", "")) for cell in grid[1][1:]]
    for cells in grid[2:]:
        if cells[0] in SKIPPED:
            continue
        perturber = PERTURBERS[cells[0]]
        points = []
        for line, cell in zip(lines, cells[1:], strict=True):
            if cell == DASH:
                continue
            value, sigma = _split(cell)
            points.append(Point(line, line, value, sigma, " | ".join(cells)))
        files[f"{material}/edwards_1991__{table}_{material}_{perturber}.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 2-9 files."""
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_2)
    _check(grid, {(1, 1): "^1H_2", (1, 2): "^2H_2"})
    files |= _by_line(
        grid,
        2,
        {1: "H2/edwards_1991__table2_H2.csv", 2: "D2/edwards_1991__table2_D2.csv"},
    )
    grid = json_table(DIRECTORY, TABLE_3)
    _check(
        grid,
        {(1, 1): "^1H_290 K", (1, 2): "90 K", (1, 3): "^2H_2200 K", (1, 4): "375 K"},
    )
    files |= _by_line(
        grid,
        2,
        {
            1: "H2/edwards_1991__table3_H2_90K.csv",
            2: "D2/edwards_1991__table3_D2_90K.csv",
            3: "D2/edwards_1991__table3_D2_200K.csv",
            4: "D2/edwards_1991__table3_D2_375K.csv",
        },
    )
    files |= _by_perturber(json_table(DIRECTORY, TABLE_4), "H2", "table4")
    files |= _by_perturber(json_table(DIRECTORY, TABLE_5), "D2", "table5")
    grid = json_table(DIRECTORY, TABLE_6)
    _check(grid, {(1, 1): "100 K", (1, 2): "150 K", (1, 3): "200 K", (1, 4): "293 K"})
    temperatures = [cell.removesuffix(" K") for cell in grid[1][1:]]
    for cells in grid[2:]:
        perturber = PERTURBERS[cells[0]]
        points = []
        for temperature, cell in zip(temperatures, cells[1:], strict=True):
            value, sigma = _split(cell)
            points.append(
                Point(f"{temperature}K", temperature, value, sigma, " | ".join(cells))
            )
        files[f"D2/edwards_1991__table6_D2_S0_{perturber}.csv"] = points
    grid = json_table(DIRECTORY, TABLE_7)
    _check(
        grid,
        {
            (1, 1): "This work",
            (1, 2): "Ref. 4",
            (1, 3): "Ref. 5",
            (1, 4): "Ref. 7",
            (1, 5): "Ref. 10",
        },
    )
    files |= _by_line(
        grid,
        2,
        {
            1: "H2/edwards_1991__table7_this_work.csv",
            2: "H2/edwards_1991__table7_may.csv",
            3: "H2/edwards_1991__table7_cooper_1968.csv",
            4: "H2/edwards_1991__table7_cooper_1970.csv",
            5: "H2/edwards_1991__table7_knaap.csv",
        },
    )
    grid = json_table(DIRECTORY, TABLE_8)
    _check(grid, {(1, 1): "This work", (1, 2): "Ref. 10"})
    files |= _by_line(
        grid,
        2,
        {
            1: "D2/edwards_1991__table8_this_work.csv",
            2: "D2/edwards_1991__table8_knaap.csv",
        },
    )
    grid = json_table(DIRECTORY, TABLE_9)
    _check(
        grid,
        {
            (1, 1): "This work",
            (1, 3): "Ref. 4",
            (1, 4): "Ref. 11",
            (2, 1): "He",
            (2, 2): "N_2",
            (2, 3): "He",
            (2, 4): "N_2",
        },
    )
    files |= _by_line(
        grid,
        3,
        {
            1: "H2/edwards_1991__table9_this_work_He.csv",
            2: "H2/edwards_1991__table9_this_work_N2.csv",
            3: "H2/edwards_1991__table9_may_He.csv",
            4: "H2/edwards_1991__table9_dayan_N2.csv",
        },
    )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

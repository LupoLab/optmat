# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bischel & Dyer (1986): H2 Q(0) and Q(1) Raman linewidths and line shifts.

Phys. Rev. A 33, 3113. The tables exist as cells files only (no CSV); every cell used
was checked against its table image.

* Table I (``tbl_64d6921655d1``): the coefficients A (MHz amagat) and B (MHz/amagat) of
  Eq. (1), Delta nu_FWHM = A/rho + B rho, for Q(J) at 81, 298 and 474 K and ortho:para
  ratios 3:1 (normal) and 1:7.7 (para-enriched), with the earlier measurements in the
  "Previous work" columns. The extraction runs each "Previous work" cell's stacked
  values together ("66.0±5.0^a70^b83.0±3.0^c"); the image shows one value per line,
  each with its reference mark, and they are split at the marks.
* Table II (``tbl_9778e6276a73``): the line-shift coefficients C (MHz/amagat) and D
  (MHz/amagat^2) of Eq. (3). The image's superscript reference marks are extracted as
  "^a" or as a bare trailing letter ("0.16±0.02a"); both are dropped. The stacked
  previous values sit in rows of their own (the row's T, J and R cells empty); they
  belong to the row above.
* Table III (``tbl_ac118c75b5a0``): the Doppler FWHM Delta nu_D (MHz) for forward
  scattering; D_0 and rho_c are in the record's notes (no quantity holds them).
* Table IV (``tbl_b44c3eba642d``): van Kranendonk's (Ref. 3) and Hunt et al.'s (Ref. 8)
  calculated widths, and the experiment; the "<2" cells at 81 K are not numbers and
  are left in the notes.

Table V (the coefficients of Eqs. (8) and (9)) is recorded as formulas.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bischel_1986_TemperatureDependenceRaman"
TABLE_I = "tbl_64d6921655d1"
TABLE_II = "tbl_9778e6276a73"
TABLE_III = "tbl_ac118c75b5a0"
TABLE_IV = "tbl_b44c3eba642d"

# One stacked "previous work" value: a number, an optional "±" error, and its
# reference mark ("^a", or a bare letter where the extraction lost the caret).
_PREVIOUS = re.compile(r"(-?[\d.]+)(?:±([\d.]+))?\^?([a-f])")


def _split(cell: str) -> tuple[str, str]:
    """Return (value, sigma) of a "v±s" or "v" cell, spaces and marks removed."""
    text = cell.replace(" ", "").rstrip("abcdef").removesuffix("^")
    value, _, sigma = text.partition("±")
    return number(value), number(sigma) if sigma else ""


def _previous(cell: str) -> dict[str, tuple[str, str]]:
    """Return the stacked previous-work values of a cell, by reference mark."""
    found = {
        m[2]: (number(m[0]), number(m[1]) if m[1] else "")
        for m in _PREVIOUS.findall(cell)
    }
    if not found and cell:
        msg = f"unreadable previous-work cell {cell!r}"
        raise ValueError(msg)
    return found


def _expect(cells: list[str], expected: list[str]) -> None:
    """Fail loudly if a row's leading cells are not what the image shows."""
    if cells[: len(expected)] != expected:
        msg = f"unexpected row {cells!r}, expected it to start {expected!r}"
        raise ValueError(msg)


def _table_i() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0][3] != "A (MHz amagat)" or grid[0][5] != "B (MHz/amagat)":
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    # (row index, T, J, R) as the image shows them; blank T cells repeat the one above.
    rows = [
        (2, "81", "0", "3:1"),
        (3, "81", "0", "1:7.7"),
        (4, "81", "1", "3:1"),
        (5, "298", "0", "3:1"),
        (6, "298", "1", "3:1"),
        (7, "474", "1", "3:1"),
    ]
    series = {
        ("0", "3:1"): "Q0_normal",
        ("0", "1:7.7"): "Q0_para",
        ("1", "3:1"): "Q1_normal",
    }
    files: dict[str, list[Point]] = {}
    previous: dict[tuple[str, str], list[Point]] = {}
    for index, temperature, j, ratio in rows:
        cells = grid[index]
        _expect(cells[1:3], [j, ratio])
        if cells[0] not in {temperature, ""}:
            msg = f"Table I row {index}: T {cells[0]!r}, expected {temperature}"
            raise ValueError(msg)
        printed = " | ".join(cells)
        name = series[j, ratio]
        for column, coefficient in ((3, "A"), (5, "B")):
            value, sigma = _split(cells[column])
            files.setdefault(
                f"H2/bischel_1986__table1_{coefficient}_{name}.csv", []
            ).append(Point(f"{temperature}K", temperature, value, sigma, printed))
        for column, coefficient in ((4, "A"), (6, "B")):
            for mark, (value, sigma) in _previous(cells[column]).items():
                previous.setdefault((coefficient, mark), []).append(
                    Point(f"Q({j})", f"Q({j})", value, sigma, printed)
                )
    # Reference marks of Table I's notes: a Ref. 22, b Ref. 6, c Ref. 8, d Ref. 9,
    # e Ref. 5, f Ref. 7.
    names = {
        "a": "murray_javan",
        "b": "allin",
        "c": "hunt",
        "d": "owyoung",
        "e": "foltz",
        "f": "lallemand",
    }
    for (coefficient, mark), points in previous.items():
        files[f"H2/bischel_1986__table1_{coefficient}_{names[mark]}.csv"] = points
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    if grid[0][3] != "C (MHz/amagat)" or grid[0][5] != "D (MHz/amagat2)":
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    # (row index, T, J, R) of the rows that start a condition; the rows after 5 and 7
    # hold only further previous-work values of the row above (image).
    rows = [
        (2, "81", "0", "3:1"),
        (3, "81", "0", "1:7.7"),
        (4, "81", "1", "3:1"),
        (5, "298", "0", "3:1"),
        (6, "298", "0", "3:1"),
        (7, "298", "1", "3:1"),
        (8, "298", "1", "3:1"),
        (9, "298", "1", "3:1"),
        (10, "474", "1", "3:1"),
    ]
    series = {
        ("0", "3:1"): "Q0_normal",
        ("0", "1:7.7"): "Q0_para",
        ("1", "3:1"): "Q1_normal",
    }
    names = {"a": "may", "b": "murray_javan", "c": "lallemand"}
    files: dict[str, list[Point]] = {}
    for index, temperature, j, ratio in rows:
        cells = grid[index]
        continuation = cells[:3] == ["", "", ""]
        if not continuation:
            _expect(cells[1:3], [j, ratio])
            if cells[0] not in {temperature, ""}:
                msg = f"Table II row {index}: T {cells[0]!r}, expected {temperature}"
                raise ValueError(msg)
        printed = " | ".join(cells)
        name = series[j, ratio]
        for column, coefficient in ((3, "C"), (5, "D")):
            if cells[column]:
                value, sigma = _split(cells[column])
                files.setdefault(
                    f"H2/bischel_1986__table2_{coefficient}_{name}.csv", []
                ).append(Point(f"{temperature}K", temperature, value, sigma, printed))
        for column, coefficient in ((4, "C"), (6, "D")):
            for mark, (value, sigma) in _previous(cells[column]).items():
                files.setdefault(
                    f"H2/bischel_1986__table2_{coefficient}_{names[mark]}.csv", []
                ).append(
                    Point(
                        f"Q({j}) {temperature} K",
                        f"Q({j}) {temperature} K",
                        value,
                        sigma,
                        printed,
                        {"temperature": temperature},
                    )
                )
    return files


def _table_iii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    if grid[0][3] != r"\Delta\nu_D(MHz)":
        msg = f"unexpected Table III header {grid[0]!r}"
        raise ValueError(msg)
    temperature = ""
    points = []
    for cells in grid[1:]:
        temperature = cells[0] or temperature
        row = f"Q({cells[1]}) {temperature} K {cells[2]}"
        points.append(
            Point(
                row,
                row,
                number(cells[3]),
                "",
                " | ".join(cells),
                {"temperature": temperature},
            )
        )
    return {"H2/bischel_1986__table3_doppler.csv": points}


def _table_iv() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_IV)
    if grid[2][3:5] != ["Ref. 3", "Ref. 8"] or grid[1][5] != "Experiment (this work)":
        msg = f"unexpected Table IV header {grid[:3]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    temperature = ""
    for cells in grid[3:]:
        temperature = cells[0] or temperature
        j = cells[1]
        # The two 81 K Q(0) rows differ in ortho:para ratio (n_J 0.25 for 3:1, 0.89
        # for 1:7.7, as in Table I); the row names say which.
        ratio = {"0.25": " 3:1", "0.89": " 1:7.7"}.get(cells[2], "")
        row = f"Q({j}) {temperature} K{ratio}"
        printed = " | ".join(cells)
        extra = {"temperature": temperature}
        for column, name in ((3, "van_kranendonk"), (4, "hunt")):
            if cells[column].startswith("<"):
                continue  # "<2" at 81 K: an upper bound, in the notes
            files.setdefault(f"H2/bischel_1986__table4_theory_{name}.csv", []).append(
                Point(row, row, number(cells[column]), "", printed, extra)
            )
        averaged = cells[5].endswith("^a")  # "Average of values in Refs. 6, 8, and 22"
        value, sigma = _split(cells[5].removesuffix("^a"))
        name = "average" if averaged else "experiment"
        files.setdefault(f"H2/bischel_1986__table4_{name}.csv", []).append(
            Point(row, row, value, sigma, printed, extra)
        )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I-IV files."""
    return _table_i() | _table_ii() | _table_iii() | _table_iv()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

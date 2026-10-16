# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Herring, Dyer & Bischel (1986): N2 and H2 rotational Raman widths, shifts.

Phys. Rev. A 34, 1944. Tables I-VI exist as cells files only (Table VII also as CSV,
recorded as formulas in the source record); every cell used was checked against its
table image. The merged header cells are extracted out of place (Table I's "295 K"
sits beside the unit, Table V's "295 K" and "80 K" run into "Ref. 20" and "This
work"); the columns are addressed by position, as the images show them:

* Table I (``tbl_efca0dbdd612``): N2 S(J) self-broadening B (MHz/amagat) at 295, 195
  and 80 K.
* Table II (``tbl_f0f74c352efe``): the room-temperature comparison: present work, Ref. 3
  (Jammu et al.), Ref. 6 (a 1964 Opt. Spectrosc. paper), Ref. 7 (Edwards et al.).
* Table III (``tbl_27787ea84e18``): O2 broadening of N2 S(J) at 295 and 195 K. The
  cells print "±" as ``\pm``.
* Table IV (``tbl_0e34a8d11d00``): H2 S(0), S(1) parameters: nu_R (cm-1), the Doppler
  FWHM (MHz) and A (MHz amagat); D_0 and rho_c are in the record's notes.
* Table V (``tbl_bc667320a0be``): H2 self-broadening B at 295 and 80 K, this work and
  Refs. 20, 21.
* Table VI (``tbl_e8a79c044c46``): H2 density shifts at 295 and 80 K, this work and
  Refs. 21, 20.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Herring_1986_TemperatureDensityDependence"
TABLE_I = "tbl_efca0dbdd612"
TABLE_II = "tbl_f0f74c352efe"
TABLE_III = "tbl_27787ea84e18"
TABLE_IV = "tbl_0e34a8d11d00"
TABLE_V = "tbl_bc667320a0be"
TABLE_VI = "tbl_e8a79c044c46"


def _split(cell: str) -> tuple[str, str]:
    r"""Return (value, sigma) of a "v±s", "v \pm s" or "v" cell."""
    value, _, sigma = cell.replace("\\pm", "±").replace(" ", "").partition("±")
    return number(value), number(sigma) if sigma else ""


def _columns(
    grid: list[list[str]], first: int, columns: dict[int, str], material: str
) -> dict[str, list[Point]]:
    """Return one file per column of S(J) rows (the row label in column 0)."""
    files: dict[str, list[Point]] = {}
    for cells in grid[first:]:
        row = f"S({number(cells[0])})"
        for column, name in columns.items():
            if not cells[column]:
                continue
            value, sigma = _split(cells[column])
            files.setdefault(f"{material}/herring_1986__{name}.csv", []).append(
                Point(row, row, value, sigma, " | ".join(cells))
            )
    return files


def _check_header(grid: list[list[str]], cells: dict[tuple[int, int], str]) -> None:
    """Fail loudly if a table's header cells are not where the image puts them."""
    for (row, column), text in cells.items():
        if grid[row][column] != text:
            msg = f"unexpected header cell ({row}, {column}) {grid[row][column]!r}"
            raise ValueError(msg)


def _table_i() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_I)
    _check_header(grid, {(0, 1): "295 K", (1, 2): "195 K", (1, 3): "80 K"})
    return _columns(
        grid, 2, {1: "table1_295K", 2: "table1_195K", 3: "table1_80K"}, "N2"
    )


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    _check_header(
        grid,
        {(0, 1): "Present work", (1, 2): "Ref. 3", (1, 3): "Ref. 6", (1, 4): "Ref. 7"},
    )
    names = {
        1: "table2_present",
        2: "table2_jammu",
        3: "table2_pinter",
        4: "table2_edwards",
    }
    return _columns(grid, 2, names, "N2")


def _table_iii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_III)
    _check_header(grid, {(1, 1): "295 K", (1, 2): "195 K"})
    return _columns(grid, 2, {1: "table3_O2_295K", 2: "table3_O2_195K"}, "N2")


def _table_iv() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_IV)
    _check_header(
        grid,
        {(0, 2): "\\nu_R(cm-1)", (0, 3): "\\Delta\\nu_D(MHz)", (0, 5): "A(MHz amagat)"},
    )
    files: dict[str, list[Point]] = {}
    temperature = ""
    for cells in grid[1:]:
        temperature = cells[0] or temperature
        line = f"S({cells[1]})"
        row = f"{line} {temperature} K"
        printed = " | ".join(cells)
        extra = {"temperature": temperature}
        for column, name in ((3, "table4_doppler"), (5, "table4_A")):
            files.setdefault(f"H2/herring_1986__{name}.csv", []).append(
                Point(row, row, number(cells[column]), "", printed, extra)
            )
        if temperature == "295":  # nu_R is printed for each temperature, the same
            files.setdefault("H2/herring_1986__table4_shift.csv", []).append(
                Point(line, line, number(cells[2]), "", printed)
            )
    return files


def _table_v() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_V)
    _check_header(
        grid, {(0, 1): "This work", (1, 2): "295 KRef. 20", (1, 3): "Ref. 21"}
    )
    _check_header(grid, {(1, 4): "80 KThis work", (1, 5): "Ref. 20"})
    names = {
        1: "table5_295K",
        2: "table5_van_den_hout_295K",
        3: "table5_cooper_295K",
        4: "table5_80K",
        5: "table5_van_den_hout_80K",
    }
    return _columns(grid, 2, names, "H2")


def _table_vi() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_VI)
    _check_header(
        grid,
        {(1, 1): "295 K", (2, 1): "This work", (2, 2): "Ref. 21", (2, 3): "This work"},
    )
    _check_header(grid, {(0, 3): "(MHz/amagat)80 K", (2, 4): "Ref. 20"})
    names = {
        1: "table6_295K",
        2: "table6_cooper_295K",
        3: "table6_80K",
        4: "table6_van_den_hout_80K",
    }
    return _columns(grid, 3, names, "H2")


def extract() -> dict[str, list[Point]]:
    """Return the Table I-VI files."""
    return (
        _table_i() | _table_ii() | _table_iii() | _table_iv() | _table_v() | _table_vi()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mikhailov (1959): N2 rotational and Q-branch Raman widths against pressure.

Sov. Phys. JETP 9, 974 (Zh. Eksp. Teor. Fiz. 36, 1368). Both tables checked against
their images.

* Table I (``tbl_e9551d96728b``, CSV): the widths delta omega (rotational lines,
  averaged over J = 0-15) and delta nu (the vibrational Q branch) in cm-1 at 27 C
  against pressure in atm, printed as two side-by-side blocks of columns (p, delta
  omega, delta nu), read here block by block. "±" is extracted as ``\pm``; the 3 atm
  cell, "~0,5*" in the image (approximate, "Measured by V. A. Iozenas"), is kept as 0.5
  and excluded in the record. "—" is no value.
* Table II (``tbl_57a9c368f39b``, cells only): the calculated rotational-structure,
  collisional and total widths of the Q branch and the measured width, at 200 and 27 C
  (the liquid row, -190 C, is not a gas value and is left in the notes). The image
  prints decimal commas in some cells ("2,34", "1,75"); they are read as points.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Mikhailov_1959_InfluenceTemperaturePressure"
TABLE_I = "tbl_e9551d96728b"
TABLE_II = "tbl_57a9c368f39b"
DASH = "\N{EM DASH}"


def _split(cell: str) -> tuple[str, str]:
    r"""Return (value, sigma) of a "v \pm s" cell, or of "\sim v^{*}"."""
    text = cell.replace("\\sim", "").replace("^{*}", "").replace("\\pm", "±")
    value, _, sigma = text.replace(" ", "").replace(",", ".").partition("±")
    return number(value), number(sigma) if sigma else ""


def _table_i() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_I)
    header = ["p, atmos", "\\delta\\omega ,cm^{-1}", "\\delta\\nu ,cm^{-1}"]
    if grid[0] != header * 2:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for block in (0, 3):
        for cells in grid[1:]:
            pressure = cells[block]
            if pressure == DASH:
                continue
            printed = " | ".join(cells[block : block + 3])
            for offset, name in ((1, "table1_rotational"), (2, "table1_vibrational")):
                cell = cells[block + offset]
                if cell == DASH:
                    continue
                value, sigma = _split(cell)
                files.setdefault(f"N2/mikhailov_1959__{name}.csv", []).append(
                    Point(f"{pressure}atm", number(pressure), value, sigma, printed)
                )
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    if grid[0][2] != "Calculated for the gas" or grid[1][5] != "\\delta\\nu, cm^{-1}":
        msg = f"unexpected Table II header {grid[:2]!r}"
        raise ValueError(msg)
    names = {
        2: "table2_rotational_structure",
        3: "table2_collisional",
        4: "table2_total",
        5: "table2_experimental",
    }
    files: dict[str, list[Point]] = {}
    for cells in grid[2:4]:  # the gas rows, 200 and 27 C
        if cells[1] != "50 atmos":
            msg = f"unexpected Table II row {cells!r}"
            raise ValueError(msg)
        row = f"{cells[0]} C"
        for column, name in names.items():
            value, sigma = _split(cells[column])
            files.setdefault(f"N2/mikhailov_1959__{name}.csv", []).append(
                Point(
                    row, row, value, sigma, " | ".join(cells), {"temperature": cells[0]}
                )
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I and II files."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

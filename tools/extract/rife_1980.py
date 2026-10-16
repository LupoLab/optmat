# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rife & Osantowski (1980): n and k of two SiO2-based glasses, 68-308 A.

Tables I (glass A, Cer-Vit C101, ``tbl_5e4d1f912c2e``) and II (glass B, Corning 7971
ULE, ``tbl_33f026d123dd``), cells only (no CSV): each prints two blocks of columns side
by side, wavelength (A), n and k, 68-188 A on the left and 198-308 A on the right. Some
cells carry an error estimate, "0.9854 ± 0.002": the number goes to ``value`` and the
estimate to ``sigma``. Four data files, n and k of each glass, in wavelength order (the
left block, then the right).

Every cell checked against the table images (2026-09-27): the cells agree. Table III
(reflectances, measured and calculated) has no vocabulary quantity and is not extracted.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rife_1980_ExtremeUltravioletOptical"
TABLES = {
    "Cer-Vit-C101/rife_1980__glass_a": "tbl_5e4d1f912c2e",
    "SiO2-TiO2/rife_1980__glass_b": "tbl_33f026d123dd",
}
HEADER = ["", "n", "k", "", "n", "k"]


def _split(cell: str) -> tuple[str, str]:
    """Return a cell's value and its error estimate ("" when none is printed)."""
    value, _, sigma = cell.partition("±")
    return number(value), number(sigma) if sigma else ""


def _glass(table_id: str) -> tuple[list[Point], list[Point]]:
    """Return the n and k points of one table, in wavelength order."""
    rows = json_table(DIRECTORY, table_id)
    if rows[1] != HEADER or rows[0][0] != "Wavelength (Å)":
        msg = f"{table_id}: unexpected header {rows[:2]!r}"
        raise ValueError(msg)
    n: list[Point] = []
    k: list[Point] = []
    for block in (slice(0, 3), slice(3, 6)):
        for cells in rows[2:]:
            wavelength, n_cell, k_cell = cells[block]
            if not wavelength:
                continue
            printed = " | ".join(cells[block])
            x = number(wavelength)
            n.append(Point(x, x, *_split(n_cell), printed))
            k.append(Point(x, x, *_split(k_cell), printed))
    return n, k


def extract() -> dict[str, list[Point]]:
    """Return n and k of glass A and glass B."""
    files: dict[str, list[Point]] = {}
    for stem, table_id in TABLES.items():
        n, k = _glass(table_id)
        files[f"{stem}_n.csv"] = n
        files[f"{stem}_k.csv"] = k
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

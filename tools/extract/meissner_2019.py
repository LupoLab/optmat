# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Meißner, Hölzer & Seeger (2019): N2 S-branch Raman widths, pure and in air.

Appl. Opt. 58, C47. Both tables as cells files only, checked against their images:
FWHM widths Gamma_J'' of the N2 S(J'') lines, "value \pm sd", one column per
temperature (the second row), headed "Gamma_J''(FWHM)/10^-3 cm" (so printed; see the
record). "-" is no value (in Table 2 the J'' lines that coincide with O2 lines are
blank in every column).

* Table 1 (``tbl_a526fa9af339``): pure N2, 295-1900 K.
* Table 2 (``tbl_2fa66b541f5b``): N2 in air ("N2-O2"), 295-1900 K.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Meissner_2019_DeterminationBranchRaman"
TABLES = {"tbl_a526fa9af339": "table1_self", "tbl_2fa66b541f5b": "table2_air"}


def _split(cell: str) -> tuple[str, str]:
    value, _, sigma = cell.replace("\\pm", "±").replace(" ", "").partition("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return one file per table and temperature."""
    files: dict[str, list[Point]] = {}
    for table_id, name in TABLES.items():
        grid = json_table(DIRECTORY, table_id)
        if not grid[0][0].startswith("\\Gamma_{J''}(FWHM)/10^{-3}"):
            msg = f"unexpected unit header {grid[0][0]!r} in {table_id}"
            raise ValueError(msg)
        header = grid[1]
        if header[0] != "J''" or not all(cell.endswith(" K") for cell in header[1:]):
            msg = f"unexpected header {header!r} in {table_id}"
            raise ValueError(msg)
        for cells in grid[2:]:
            row = f"S({number(cells[0])})"
            for head, cell in zip(header[1:], cells[1:], strict=True):
                if cell == "-":
                    continue
                value, sigma = _split(cell)
                temperature = head.removesuffix(" K")
                files.setdefault(
                    f"N2/meissner_2019__{name}_{temperature}K.csv", []
                ).append(Point(row, row, value, sigma, " | ".join(cells)))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

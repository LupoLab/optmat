# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hölzer, Meißner & Seeger (2021): O2 S-branch Raman widths, self and in air.

Appl. Opt. 60, C76. Both tables as CSV, checked against their images: FWHM widths
Gamma_N'' in 10^-3 cm-1 atm-1 of the S(N'') lines, "value ± sd", one column per
temperature. The unit header is split over the first row's cells ("0N00(FWH",
"M)/10-3 cm-1", "atm-1": Gamma_N''(FWHM)/10^-3 cm^-1 atm^-1 in the image) and the row
label N'' is extracted as "N00"; the temperatures are the second row.

* Table 1 (``tbl_2ed5b19b32b8``): pure O2, 295-1900 K.
* Table 2 (``tbl_2e8c929f073a``): O2 in air, 295-1900 K.
"""

from _library import Point, number, table, write

DIRECTORY = "Holzer_2021_ImprovementCoherentModel"
TABLES = {"tbl_2ed5b19b32b8": "table1_self", "tbl_2e8c929f073a": "table2_air"}


def _split(cell: str) -> tuple[str, str]:
    value, _, sigma = cell.replace(" ", "").partition("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return one file per table and temperature."""
    files: dict[str, list[Point]] = {}
    for table_id, name in TABLES.items():
        grid = table(DIRECTORY, table_id)
        header = grid[1]
        if header[0] != "N00" or not all(cell.endswith("K") for cell in header[1:]):
            msg = f"unexpected header {header!r} in {table_id}"
            raise ValueError(msg)
        for cells in grid[2:]:
            row = f"S({number(cells[0])})"
            for head, cell in zip(header[1:], cells[1:], strict=True):
                if not cell:
                    continue
                value, sigma = _split(cell)
                files.setdefault(f"O2/holzer_2021__{name}_{head}.csv", []).append(
                    Point(row, row, value, sigma, " | ".join(cells))
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

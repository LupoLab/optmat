# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tarazkar, Romanov & Levis (2016): static n2 of the neutral noble gases.

Table II (``tbl_53f14a8f3c08``, CSV; image checked, every cell agrees): n2 in 10^-9
cm^2/TW for He, Ne, Ar, Kr and their ions, in the static limit and at a species-specific
shortest wavelength (not stated numerically). Only the neutral atoms' static column is
written (the ions are not curated species; the lambda_min column has no wavelength).
"""

from _library import Point, table, write

DIRECTORY = "Tarazkar_2016_NonlinearOpticalResponses"
TABLE_II = "tbl_53f14a8f3c08"
NEUTRALS = ("He", "Ne", "Ar", "Kr")


def extract() -> dict[str, list[Point]]:
    """Return the static n2 of each neutral atom."""
    grid = table(DIRECTORY, TABLE_II)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        if cells[0] not in NEUTRALS:
            continue
        value = cells[1]
        float(value)
        files[f"{cells[0]}/tarazkar_2016__n2_static_{cells[0]}.csv"] = [
            Point(cells[0], "inf", value, "", " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

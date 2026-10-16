# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Abjean, Mehu & Johannin-Gilles (1971): Xe and Kr n - 1 at 180.7-253.7 nm.

Tables 1 (xenon) and 2 (krypton), checked on the table images: the measured column and
the column calculated with the paper's formulas, ``(n - 1) x 10^6``; the difference
column is not extracted. Table 2's last calculated cell is printed ``507.7``, as the
image shows: a misprint corrected by an erratum in the source record.

Tables 3 (xenon) and 4 (krypton) have no CSV rendering (merged header) and are read
from their JSON cells. Their "present work" column is the paper's formula evaluated at
round wavelengths (it reproduces the formulas to 0.05e-6). The other columns reproduce
Kronjager (ref. 1), Koch (ref. 2), Chashchina and Shreider (ref. 3), Kingston (ref. 9)
and Langhoff and Karplus (ref. 10). The table note says "a) Experimental values": the
cells of refs. 1-3 without the mark are the other authors' formulas extrapolated, which
are not extracted (as in ``bideau_mehu_1981``). Refs. 9 and 10 are calculations, taken
whole; "-" is no value. Table 3's heading prints ``(n-1) . 10^3`` (the image shows the
3); the values are ``x 10^6`` as in Table 1 (see the source record).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Abjean_1971_InterferometricDeterminationXenon"
MEASURED = {"Xe": "tbl_84f9d5968312", "Kr": "tbl_00feebf8a634"}
COMPARED = {"Xe": "tbl_30db947f05b2", "Kr": "tbl_9e525b4dedbd"}

# Column index of Tables 3 and 4 -> dataset stem, and whether only the cells marked
# "a)" (experimental) are taken.
COMPARED_COLUMNS = {
    1: ("formula_table", False),
    2: ("kronjager_1935", True),
    3: ("koch_1949", True),
    4: ("chashchina_shreider_1969", True),
    5: ("kingston_1964", False),
    6: ("langhoff_karplus_1969", False),
}
MARK = "a)"


def extract() -> dict[str, list[Point]]:
    """Return the measured, calculated and reproduced columns per gas."""
    files: dict[str, list[Point]] = {}
    for gas, table_id in MEASURED.items():
        rows = table(DIRECTORY, table_id)[1:]
        for column, stem in ((1, "experiment"), (2, "fit")):
            files[f"{gas}/abjean_1971__{stem}_{gas}.csv"] = [
                Point(
                    cells[0].strip(),
                    cells[0].strip(),
                    number(cells[column]),
                    "",
                    " | ".join(cells),
                )
                for cells in rows
            ]
    for gas, table_id in COMPARED.items():
        rows = json_table(DIRECTORY, table_id)[2:]  # caption row, column heads
        for column, (stem, marked_only) in COMPARED_COLUMNS.items():
            points = []
            for cells in rows:
                wavelength, cell = cells[0].strip(), cells[column].strip()
                marked = cell.endswith(MARK)
                cell = cell.removesuffix(MARK).strip()
                if cell in {"", "-"} or (marked_only and not marked):
                    continue
                points.append(
                    Point(wavelength, wavelength, number(cell), "", " | ".join(cells))
                )
            files[f"{gas}/abjean_1971__{stem}_{gas}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

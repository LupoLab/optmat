# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Santran et al. (2004): n2 of ten fused silicas at 800 and 400 nm (Table 2).

Table 2 is published only as JSON cells (merged header cells). Its rows are the samples
(one grade each); its two value columns are the nonlinear index at 800 and at 400 nm.
Each grade is a dataset of its own (the grade is a dataset parameter).

Structural fix, read from the table image (``tables/tbl_ece90cb0d86a.jpg``): the unit
cells, extracted as "(10216 cm2/W)", print (10^-16 cm^2/W) (the pdf's font maps the
minus to "2"); the record's scale is 1e-16. The header rows (the first three) carry no
values and are skipped.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Santran_2004_PreciseAbsoluteMeasurements"
TABLE = "tbl_ece90cb0d86a"
HEADER_ROWS = 3
UNIT_CELL = "(10216 cm2/W)"
COLUMNS = {1: "800", 2: "400"}

# Sample, as printed -> dataset id (the grade parameter keeps the printed name).
SAMPLES = {
    "Heraeus S300": "heraeus_s300",
    "Herasil": "herasil",
    "Suprasil": "suprasil",
    "Heraeus Homosil": "heraeus_homosil",
    "Heraeus S1": "heraeus_s1",
    "Herasil S1V": "herasil_s1v",
    "Suprasil EN1027A": "suprasil_en1027a",
    "F851053": "f851053",
    "Heraeus H1": "heraeus_h1",
    "Schott SQ1": "schott_sq1",
}


def extract() -> dict[str, list[Point]]:
    """Return one file per fused-silica sample."""
    rows = json_table(DIRECTORY, TABLE)
    if rows[HEADER_ROWS - 1][1:] != [UNIT_CELL, UNIT_CELL]:
        msg = f"unexpected unit cells {rows[HEADER_ROWS - 1]!r}"
        raise ValueError(msg)
    files = {}
    for cells in rows[HEADER_ROWS:]:
        name = SAMPLES[cells[0]]
        files[f"SiO2/santran_2004__{name}.csv"] = [
            Point(wavelength, wavelength, number(cells[column]), "", " | ".join(cells))
            for column, wavelength in COLUMNS.items()
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

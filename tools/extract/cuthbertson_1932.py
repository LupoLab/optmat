# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract C. and M. Cuthbertson (1932): neon and helium n - 1 at 275-546 nm.

Table II (``tbl_4ee024d8daf8``, neon) and Table III (``tbl_01cd89e4e3aa``, helium)
have no CSV in the library and are read from their cells files. Each gives, per
wavelength (``lambda x 10^8``, i.e. Angstrom: "in vacuo" for neon, "four figures" for
helium), the observed refractivity ``(mu - 1) x 10^8``, the value calculated from the
paper's one-term formula and their difference in the ninth decimal place. The observed
column is the measurement; the calculated one is extracted as a check of the formula.
Neon's second column (relative numbers of bands) is not extracted.

The tables print the decimal point as a raised dot (``6725·0``); it is read as a
decimal point. The table images agree with the extracted cells.

Not extracted: the unlabelled table ``tbl_9c520fc15186`` (four series of the neon
refractivity at 5462.23 A whose mean, 6725e-8, is Table II's first row) and the three
helium values in the text (p. 45), whose mean is Table III's first row; and Tables I
and IV (dispersion-formula constants), which are formulas.
"""

import json

from _library import LIBRARY, Point, number, write

DIRECTORY = "Cuthbertson_1932_RefractionDispersionNeon"
TABLES = {"Ne": "tbl_4ee024d8daf8", "He": "tbl_01cd89e4e3aa"}
# gas -> (observed column, calculated column) in the cells grid.
COLUMNS = {"Ne": (2, 3), "He": (1, 2)}
HEADER_ROWS = 2  # a spanning header and the "Observed. / Calculated." row


def _cells_grid(directory: str, table_id: str) -> list[list[str]]:
    """Return a table's rows from its cells file, for tables without a CSV."""
    path = LIBRARY / directory / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    n_rows = 1 + max(c["row"] for c in cells)
    n_columns = max(c["column"] + c["column_span"] for c in cells)
    grid = [[""] * n_columns for _ in range(n_rows)]
    for cell in cells:
        grid[cell["row"]][cell["column"]] = cell["text"]
    return grid


def _decimal(cell: str) -> str:
    """Read the printed raised decimal point as a decimal point."""
    return number(cell.replace("\N{MIDDLE DOT}", "."))


def extract() -> dict[str, list[Point]]:
    """Return the observed and calculated columns of Tables II and III."""
    files = {}
    for gas, table_id in TABLES.items():
        observed_column, calculated_column = COLUMNS[gas]
        observed, calculated = [], []
        for cells in _cells_grid(DIRECTORY, table_id)[HEADER_ROWS:]:
            wavelength = _decimal(cells[0])
            printed = " | ".join(cells)
            value = _decimal(cells[observed_column])
            fitted = _decimal(cells[calculated_column])
            observed.append(Point(wavelength, wavelength, value, "", printed))
            calculated.append(Point(wavelength, wavelength, fitted, "", printed))
        files[f"{gas}/cuthbertson_1932__experiment_{gas}.csv"] = observed
        files[f"{gas}/cuthbertson_1932__fit_{gas}.csv"] = calculated
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

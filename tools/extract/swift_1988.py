# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Swift, Schlie & Rathge (1988): polarizabilities at 458-1315 nm (Table I).

Table I (``tbl_98326ad22a2d``, cells file only; image checked) gives, per gas and
wavelength (A: 4579, 4880, 5145, 6328, 13150), the polarizability in 1e-24 cm^3 per
molecule followed in parentheses by the Gladstone-Dale constant beta = 2 pi alpha N_s
(1e-4), and a last column of "maximum uncertainty" in percent. The library's header
row is split over two rows (the wavelengths are the second); the numbers are intact.
Only alpha is extracted (beta is alpha times a constant, which holds for every cell of
the curated gases); the gases not curated (H2O, Cl2, ICl, CF3I, I2, n-C3F7I) are
skipped.
"""

import json
import re

from _library import LIBRARY, Point, number, write

DIRECTORY = "Swift_1988_DispersionGasesAtomic"
TABLE = "tbl_98326ad22a2d"
GASES = {
    "He": "He",
    "N_2": "N2",
    "O_2": "O2",
    "Ar": "Ar",
    "CO_2": "CO2",
    "O_3": "O3",
    "SF_6": "SF6",
}
_CELL = re.compile(r"^([\d.]+)\s*\(([\d.]+)\)$")


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


def extract() -> dict[str, list[Point]]:
    """Return alpha per curated gas at the five wavelengths."""
    _, wavelengths, *rows = _cells_grid(DIRECTORY, TABLE)
    files = {}
    for cells in rows:
        gas = GASES.get(cells[0].strip())
        if gas is None:
            continue
        points = []
        for column in range(1, 6):
            match = _CELL.match(cells[column].strip())
            if match is None:
                msg = f"cannot read {cells[column]!r} for {gas}"
                raise ValueError(msg)
            wavelength = wavelengths[column].strip()
            points.append(
                Point(
                    wavelength,
                    wavelength,
                    number(match.group(1)),
                    "",
                    " | ".join(cells),
                )
            )
        files[f"{gas}/swift_1988__experiment_{gas}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

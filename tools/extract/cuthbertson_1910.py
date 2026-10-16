# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract C. and M. Cuthbertson (1910): Ar, He, Ne, Kr, Xe dispersion at 480-671 nm.

All the paper's own values are *doubled* refractivities, 2 (mu - 1), "in order to
facilitate comparison with the diatomic gases" (p. 13); the datasets' scales halve
them. Wavelengths (``lambda x 10^8``, Angstrom) are the air wavelengths of the Cd, Hg,
Ag and Li lines. The tables print the decimal point as a raised dot, read here as a
decimal point (and the argon table's ``6438 ·5`` loses its space).

- Table I (``tbl_7e2ecfe6276e``), argon: observed, calculated, difference. The
  library's header reads ``(mu - 1) 10^9``; the table image prints ``10^8``, which the
  values confirm (56460 x 1e-8 / 2 = 2.8230e-4, the text's 1.00028230). The 5769.5 A
  value, which the authors left out of their fit as "evidently affected by error",
  goes to an excluded dataset.
- Table II (``tbl_c0066726ecc8``, cells file only), helium, neon, krypton, xenon,
  observed and calculated, headed ``(mu - 1) x 10^8``. The xenon columns are
  ``x 10^7`` (14109.8 x 1e-7 / 2 = 7.05e-4; the paper's xenon formula gives the
  printed calculated column only on that scale): their scale says so. The helium value
  printed at 5769.5 A, 6963.8, is 19.7 below the calculated 6983.5 and below the
  6438.5 A value; it goes to an excluded dataset (see the source record).
- The unlabelled table of earlier argon work (``tbl_f8dd6a8e65a1``, cells file only,
  ``(mu - 1) 10^7``, not doubled): Burton's and Ahrberg's values with a wavelength are
  extracted as reproduced datasets; the white-light values (Rayleigh, Ramsay and
  Travers) and Ahrberg's "5790-70" (the unresolved yellow mercury pair) are not.

The constants tables (``tbl_70d50384af0d``, ``tbl_f14deb352ff6``) are formulas,
recorded in the source record.
"""

import json

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Cuthbertson_1910_RefractionDispersionArgon"
ARGON = "tbl_7e2ecfe6276e"
OTHERS = "tbl_c0066726ecc8"
EARLIER = "tbl_f8dd6a8e65a1"
# Table II: gas -> (observed column, calculated column).
TABLE_II_COLUMNS = {"He": (1, 2), "Ne": (3, 4), "Kr": (5, 6), "Xe": (7, 8)}
# Rows set aside: (gas, wavelength as printed).
EXCLUDED = {("Ar", "5769.5"), ("He", "5769.5")}
# Earlier argon work: observer (first row of their block) -> dataset, and the
# wavelengths taken.
EARLIER_ROWS = {
    "Burton": ("burton_1908_Ar", {"6563", "5896"}),
    "Ahrberg": ("ahrberg_1909_Ar", {"6439", "5461", "4359"}),
}


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
    """Read the printed raised decimal point (and a stray space) as a decimal point."""
    return number(cell.replace("\N{MIDDLE DOT}", ".").replace(" ", ""))


def _columns(
    gas: str, rows: list[list[str]], observed_column: int, calculated_column: int
) -> dict[str, list[Point]]:
    observed, excluded, calculated = [], [], []
    for cells in rows:
        if cells[observed_column].strip() in {"—", "-", ""}:
            continue
        wavelength = _decimal(cells[0])
        printed = " | ".join(cells)
        value = _decimal(cells[observed_column])
        point = Point(wavelength, wavelength, value, "", printed)
        (excluded if (gas, wavelength) in EXCLUDED else observed).append(point)
        fitted = _decimal(cells[calculated_column])
        calculated.append(Point(wavelength, wavelength, fitted, "", printed))
    prefix = f"{gas}/cuthbertson_1910__"
    files = {
        f"{prefix}experiment_{gas}.csv": observed,
        f"{prefix}fit_{gas}.csv": calculated,
    }
    if excluded:
        files[f"{prefix}experiment_{gas}_excluded.csv"] = excluded
    return files


def _earlier() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    observer = ""
    for cells in _cells_grid(DIRECTORY, EARLIER)[1:]:
        # The observer is named on the first row of their block only.
        if cells[0].strip():
            observer = cells[0].split()[0].rstrip("*†‡§")
        if observer not in EARLIER_ROWS:
            continue
        dataset, wavelengths = EARLIER_ROWS[observer]
        wavelength = cells[1].strip()
        if wavelength not in wavelengths:
            continue
        files.setdefault(f"Ar/cuthbertson_1910__{dataset}.csv", []).append(
            Point(wavelength, wavelength, number(cells[2]), "", " | ".join(cells))
        )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I and II and the earlier argon values as data files."""
    files = _columns("Ar", table(DIRECTORY, ARGON)[1:], 1, 2)
    rows = _cells_grid(DIRECTORY, OTHERS)[2:]  # two header rows
    for gas, (observed_column, calculated_column) in TABLE_II_COLUMNS.items():
        files |= _columns(gas, rows, observed_column, calculated_column)
    files |= _earlier()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

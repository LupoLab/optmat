# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mansfield & Peck (1969): helium n - 1 at 480-2059 nm and comparison tables.

Table I (``tbl_aae92dff687e``, no CSV in the library: read from its cells file) gives
the measured refractivity ``(n - 1) x 10^9`` at 760 torr, 0 degrees C, at fifteen
vacuum wavelengths (Angstrom), with the deviation from the paper's Eq. (4). The table
image prints two values in parentheses with a footnote ``a``, "Values not included
while fitting dispersion formulas" (9660.434 and 15300.153 A); those two go to a
separate, excluded dataset. The deviation column is not extracted (it is a residual,
and it does not reproduce from the printed Eq. (4); see the source record).

Table II (``tbl_330bf663e80e``) lists oscillator strengths and wavenumbers (in
um^-1) of the 1s -> np, (2s, np) and (3s, np) transitions, all taken from other
authors or estimated; the footnote letters are stripped. The image shows the row
labels of the three series; rows get ids naming the upper state.

Tables III, IV and VI are one-value-per-author comparisons at 5462.258 A (III, IV)
and in the static limit (VI, the polarizability column); each author becomes its own
dataset. "This work" in Table III repeats Table I's 5462.258 A value and is not
extracted again. Tables V, VII and VIII (effective F values, Becquerel constants and
Verdet constants derived from dispersion formulas) are not data and are not extracted.
"""

import json
import re

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Mansfield_1969_DispersionHelium"
TABLE_I = "tbl_aae92dff687e"
TABLE_II = "tbl_330bf663e80e"
TABLE_III = "tbl_cd9d9ebee3ff"
TABLE_IV = "tbl_acc0ad617db8"
TABLE_VI = "tbl_a09b6733c084"
PREFIX = "He/mansfield_peck_1969__"
GREEN = "5462.258"

# Table II row ids, in the printed order, read from the image: the 1s np series
# (n = 2..10, then the sum from n = 11 lumped at 505.5 A), the doubly excited
# (2s, np) series (n = 2..6 and the sum from 7) and (3s, np) (n = 2..4, sum from 5).
TABLE_II_ROWS = (
    *(f"1s{n}p" for n in range(2, 11)),
    "1snp_11+",
    *(f"2s{n}p" for n in range(2, 7)),
    "2snp_7+",
    *(f"3s{n}p" for n in range(2, 5)),
    "3snp_5+",
)

# Table III rows (author -> dataset); "This work" is Table I's value again.
TABLE_III_ROWS = {"Koch": "koch_1913", "Cuthbertson and Cuthbertson": "cuthbertson"}
# Table IV rows (author -> dataset, wavelength in A).
TABLE_IV_ROWS = {
    "Korwin": ("korwin_1934", GREEN),
    "Present calc.": ("theory", GREEN),
    "Musulin and Epstein": ("musulin_epstein_1964", GREEN),
    "Karplus and Kolker": ("karplus_kolker_1963", GREEN),
    "Migneron and Levinger (5500 Å)": ("migneron_levinger_1965", "5500"),
}
# Table VI rows (author -> dataset); the polarizability column.
TABLE_VI_ROWS = {
    "Cuthbertson and Cuthbertson": "alpha_cuthbertson",
    "Present expt.": "alpha_experiment",
    "Present theor. work": "alpha_theory",
    "Chung and Hurst": "alpha_chung_hurst_1966",
    "Maryott and Buckley": "alpha_maryott_buckley_1953",
    "Chung^{\\text{c}}": "alpha_chung_1968",
}

_FOOTNOTE = re.compile(r"\^\{?(\\text\{)?[a-e]\}?\}?$")


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


def _strip_mark(cell: str) -> str:
    r"""Drop a trailing footnote letter (``^a``, ``^{\text{b}}``) from a cell."""
    return _FOOTNOTE.sub("", cell.strip())


def _table_i() -> dict[str, list[Point]]:
    fitted, excluded = [], []
    for cells in _cells_grid(DIRECTORY, TABLE_I)[1:]:
        wavelength, value = cells[1].strip(), cells[2].strip()
        printed = " | ".join(cells)
        # "(34759.27)^a": parenthesised, footnote a = not used in the fits.
        if value.startswith("("):
            number_ = number(value.removesuffix("^a").strip("()"))
            excluded.append(Point(wavelength, wavelength, number_, "", printed))
        else:
            fitted.append(Point(wavelength, wavelength, number(value), "", printed))
    return {
        f"{PREFIX}experiment.csv": fitted,
        f"{PREFIX}experiment_excluded.csv": excluded,
    }


def _table_ii() -> list[Point]:
    rows = table(DIRECTORY, TABLE_II)[1:]
    if len(rows) != len(TABLE_II_ROWS):
        msg = f"Table II has {len(rows)} rows, expected {len(TABLE_II_ROWS)}"
        raise ValueError(msg)
    return [
        Point(
            row, number(cells[1]), number(_strip_mark(cells[2])), "", " | ".join(cells)
        )
        for row, cells in zip(TABLE_II_ROWS, rows, strict=True)
    ]


def _one_per_author(
    table_id: str, authors: dict[str, str], column: int, wavelength: str
) -> dict[str, list[Point]]:
    files = {}
    for cells in table(DIRECTORY, table_id)[1:]:
        author = cells[0].strip()
        if author not in authors:
            continue
        # "3.4890 \pm 0.05 (?)": the value is the first number; the stated error
        # (and the paper's query mark) stay in the printed cells.
        value = number(_strip_mark(cells[column].split("\\pm")[0]))
        files[f"{PREFIX}{authors[author]}.csv"] = [
            Point("1", wavelength, value, "", " | ".join(cells))
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I, II, III, IV and VI as data files."""
    files = _table_i()
    files[f"{PREFIX}oscillator_strengths.csv"] = _table_ii()
    files |= _one_per_author(TABLE_III, TABLE_III_ROWS, 1, GREEN)
    for cells in table(DIRECTORY, TABLE_IV)[1:]:
        dataset, wavelength = TABLE_IV_ROWS[cells[0].strip()]
        files[f"{PREFIX}{dataset}.csv"] = [
            Point("1", wavelength, number(cells[1]), "", " | ".join(cells))
        ]
    files |= _one_per_author(TABLE_VI, TABLE_VI_ROWS, 2, "inf")
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

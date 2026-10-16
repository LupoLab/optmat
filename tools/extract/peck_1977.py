# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Peck & Huang (1977): H2 refractivity at 0.40-1.69 um (Tables I, II and IV).

Table I (``tbl_7beee2e55930``, checked against the image) gives 17 vacuum wavelengths
(um, printed with a thin space after the fourth decimal), the measured
``10^6 (n - 1)`` and the value of the paper's Eq. (4). The image prints every decimal
point as a comma-like low dot; the library extraction converted all but one cell
(``137,746``), which is read here as ``137.746``. Two rows are printed in parentheses,
the opening one before the wavelength and the closing one after the last column: they
were excluded from the paper's fits (p. 1551) and are written to a separate file.
The Eq. (4) column's ``135.305`` at 1.3507884 um is a misprint for 136.305 (an
erratum in the source record).

Table II (``tbl_f3f6edf512d1``) reproduces older refractivities at 0.5462252 um from
Cuthbertson & Cuthbertson (1933); its last row, 'This work', repeats Eq. (3) and is
skipped. Table IV (``tbl_c46a760f2d3f``, JSON cells only) gives polarizabilities in
units of a0^3 at infinite wavelength and at 0.5462252 um, one column per origin; the
'...' cell (no Kolos-Wolniewicz value at 0.5462 um) is skipped.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Peck_1977_RefractivityDispersionHydrogen"
TABLE_I = "tbl_7beee2e55930"
TABLE_II = "tbl_f3f6edf512d1"
TABLE_IV = "tbl_c46a760f2d3f"

TABLE_IV_COLUMNS = {
    1: "alpha_this_work",
    2: "alpha_kirn",
    3: "alpha_kolos_wolniewicz",
    4: "alpha_ford_browne",
}


def _cell(text: str) -> str:
    """Return a Table I number: parentheses, thin spaces and the decimal comma gone."""
    return number(text.strip("()").replace(" ", "").replace(",", "."))


def _table_i() -> dict[str, list[Point]]:
    measured: list[Point] = []
    excluded: list[Point] = []
    fit: list[Point] = []
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        wavelength = _cell(cells[0])
        printed = " | ".join(cells)
        point = Point(wavelength, wavelength, _cell(cells[1]), "", printed)
        (excluded if cells[0].startswith("(") else measured).append(point)
        fit.append(Point(wavelength, wavelength, _cell(cells[2]), "", printed))
    return {
        "H2/peck_1977__experiment.csv": measured,
        "H2/peck_1977__experiment_excluded.csv": excluded,
        "H2/peck_1977__eq4.csv": fit,
    }


def _table_ii() -> dict[str, list[Point]]:
    points = [
        Point(
            f"{cells[1]} {cells[0]}",
            "0.5462252",
            number(cells[2]),
            "",
            " | ".join(cells),
        )
        for cells in table(DIRECTORY, TABLE_II)[1:]
        if cells[1] != "This work"
    ]
    return {"H2/peck_1977__literature_546.csv": points}


def _table_iv() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_IV)[2:]
    for column, dataset in TABLE_IV_COLUMNS.items():
        points = []
        for cells in rows:
            wavelength = "inf" if cells[0] == "\N{INFINITY}" else _cell(cells[0])
            if cells[column] in {"", "...", "\N{HORIZONTAL ELLIPSIS}"}:
                continue
            points.append(
                Point(
                    wavelength, wavelength, number(cells[column]), "", " | ".join(cells)
                )
            )
        files[f"H2/peck_1977__{dataset}.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I, II and IV data files."""
    return _table_i() | _table_ii() | _table_iv()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

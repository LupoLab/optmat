# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Weaver, Krafka, Lynch & Koch (1981): the aluminium and gold tables.

The library holds the book's front matter (``original/source_01``) and the Ag, Al and
Au chapters and the bibliography (``source_02`` to ``source_05``), whose tables are in
``supplement_01`` (Ag), ``supplement_02`` (Al), ``supplement_03`` (Au) and
``supplement_04`` (bibliography, no tables). Every table prints, against the photon
energy (eV), epsilon_1, epsilon_2, n, k, Im(-1/epsilon) and R(phi = 0); n and k are
extracted as the primary data and epsilon_1, epsilon_2 as checks, and the loss
function and R stay in ``printed`` with the row's other cells.

**Aluminium** (pp. 79-86, "publication by E. Shiles, T. Sasaki, M. Inokuti, and D.Y.
Smith in Phys. Rev. B 22, 1612 (1980) based on the following tabulation"): eight
tables, one per page, in page order (``ALUMINIUM``), 0.040-1551.000 eV, 447 rows, three
decimals (R four). The pages continue one another (0.190/0.195, 1.460/1.480,
14.400/14.500, 25.500/26.000, 72.400/72.500, 111.000/112.000, 340.000/350.000 eV) and
the energies increase throughout. Every page was checked against its image, and every
row with the redundant columns (epsilon = N^2, the loss function from epsilon, R from
N) to the printed rounding; the three cells that failed were read against the images:

- 1.700 eV, R: the library reads ".8977", the image prints ".8972" (structural fix);
- 27.000 eV, epsilon_1: the library reads "0.700", the image prints "0.706" (fix);
- 134.000 eV, Im(-1/epsilon): the library reads "0.036", the image "0.038" (fix);
- 17.500 eV, k: the image prints "0.036", as read; a misprint corrected by an erratum
  in the record (not here).

**Gold** (pp. 64-69, "Unpublished reflectance data for electropolished Au(110) by J.H.
Weaver (0.2 <= hv <= 4 eV) and C.G. Olson and D.W. Lynch (~4 <= hv <= 30 eV)"): six
tables in page order (``GOLD``), 0.10-30.00 eV, 284 printed rows, two decimals (R
three). Checked cell by cell against an independent transcription of the images (a
scratch file) and with the redundant columns; the library misreads four cells:

- 2.20 eV, R: ".980" for ".880";
- 2.50 eV, epsilon_1: "-1.96" for "-1.86";
- 19.40 eV, k: "0.76" for "0.78";
- 28.80 eV, epsilon_2: "0.65" for "0.85".

The book prints the 8.40 eV row twice (p. 66; the record excludes the repeat) and
"6.50" twice, the second for 6.60 (p. 66); their second rows get the ids "8.40
(repeat)" and "6.50 (repeat)" so that row ids stay unique; "26.50" (p. 68) is printed
for 26.60. The record corrects the energies by errata.

**The workbook's gold** (``optical-constantsCuAgAu.xlsx``, attached to Babar & Weaver's
paper in the library, sheet ``JHW-Au-ZAED``, read with :mod:`_xlsx`): Babar & Weaver's
copy of the gold table, 217
rows, 0.1-30 eV, E (eV), R, n, k, e1, e2 (datasets ``sheet_*``, checks). Its e1, e2
are the book's, its n, k sqrt(e1 + i e2) of them to five figures; it differs from the
book at 12.2, 17.6 and 21.4 eV (epsilon) and 11.4 eV (R) (see the record).

The silver table (pp. 52-53, Hagemann, Gudat & Kunz's) is not extracted.
"""

from _library import Point, attachment, number, table, write
from _xlsx import sheet

DIRECTORY = "Weaver_1981_OpticalPropertiesMetals"
HEADER_START = "Energy (eV)"

# Book page -> library table, in page order.
ALUMINIUM = {
    79: "tbl_41f659cbf54e",
    80: "tbl_41f03578c473",
    81: "tbl_550bc2056f85",
    82: "tbl_6b7f728e7af0",
    83: "tbl_5937497d4313",
    84: "tbl_7ae05093ee33",
    85: "tbl_b1243e5e0609",
    86: "tbl_bd839bc39e50",
}
GOLD = {
    64: "tbl_e1aa4a8fd3ac",
    65: "tbl_0099904bfe57",
    66: "tbl_3b6da8494ea6",
    67: "tbl_c59b9057cc82",
    68: "tbl_22fe363af672",
    69: "tbl_e323eaae89b9",
}
ROWS = {"Al": 447, "Au": 284}

# Library misreadings, checked against the table images: (metal, energy as printed,
# column) -> (the library's cell, the image's). Columns: 0 E, 1 epsilon_1, 2 epsilon_2,
# 3 n, 4 k, 5 Im(-1/epsilon), 6 R.
FIXES = {
    ("Al", "1.700", 6): (".8977", ".8972"),
    ("Al", "27.000", 1): ("0.700", "0.706"),
    ("Al", "134.000", 5): ("0.036", "0.038"),
    ("Au", "2.20", 6): (".980", ".880"),
    ("Au", "2.50", 1): ("-1.96", "-1.86"),
    ("Au", "19.40", 4): ("0.76", "0.78"),
    ("Au", "28.80", 2): ("0.65", "0.85"),
}

# Dataset suffix -> column of the book's tables.
COLUMNS = {"n": 3, "k": 4, "eps1": 1, "eps2": 2}

WORKBOOK_PAPER = "Babar_2015_OpticalConstantsCu"
WORKBOOK = "optical-constantsCuAgAu.xlsx"
SHEET = "JHW-Au-ZAED"
SHEET_HEADER = ["E (eV)", "R", "n", "k", "e1", "e2"]
SHEET_ROWS = 217
# Dataset suffix -> column of the sheet.
SHEET_COLUMNS = {"n": 2, "k": 3, "eps1": 4, "eps2": 5}


def _book(metal: str, pages: dict[int, str], supplement: str) -> list[list[str]]:
    """Return a metal's table rows in page order, with the misreadings fixed."""
    rows: list[list[str]] = []
    for page, table_id in pages.items():
        cells = table(f"{DIRECTORY}/{supplement}", table_id)
        if cells[0][0] != HEADER_START or len(cells[0]) != 7:
            msg = f"{metal} p. {page} ({table_id}): unexpected header {cells[0]}"
            raise ValueError(msg)
        rows += [[cell.strip() for cell in row] for row in cells[1:]]
    if len(rows) != ROWS[metal] or any(len(row) != 7 for row in rows):
        msg = f"{metal}: {len(rows)} rows, not {ROWS[metal]} of seven cells"
        raise ValueError(msg)
    for (fix_metal, energy, column), (library, image) in FIXES.items():
        if fix_metal != metal:
            continue
        (row,) = [row for row in rows if row[0] == energy]
        if row[column] != library:
            msg = f"{metal} {energy}: expected the misreading {library!r}"
            raise ValueError(msg)
        row[column] = image
    return rows


def _points(metal: str, rows: list[list[str]]) -> dict[str, list[Point]]:
    """Return the book's datasets; a repeated printed energy gets a suffixed row id."""
    ids: list[str] = []
    for row in rows:
        ids.append(f"{row[0]} (repeat)" if row[0] in ids else row[0])
    return {
        f"{metal}/weaver_1981__{name}_{metal}.csv": [
            Point(row_id, number(row[0]), number(row[column]), "", " | ".join(row))
            for row_id, row in zip(ids, rows, strict=True)
        ]
        for name, column in COLUMNS.items()
    }


def _sheet() -> dict[str, list[Point]]:
    """Return the workbook's copy of the gold table."""
    rows = sheet(attachment(WORKBOOK_PAPER, WORKBOOK), SHEET)
    if rows[0] != SHEET_HEADER:
        msg = f"{SHEET}: unexpected header {rows[0]}"
        raise ValueError(msg)
    body = [cells for cells in rows[1:] if cells]
    if len(body) != SHEET_ROWS or any(len(c) != 6 for c in body):
        msg = f"{SHEET}: {len(body)} rows, not {SHEET_ROWS} of six cells"
        raise ValueError(msg)
    if len({c[0] for c in body}) != SHEET_ROWS:
        msg = f"{SHEET}: repeated energies"
        raise ValueError(msg)
    return {
        f"Au/weaver_1981__sheet_{name}_Au.csv": [
            Point(c[0], number(c[0]), number(c[column]), "", " | ".join(c))
            for c in body
        ]
        for name, column in SHEET_COLUMNS.items()
    }


def extract() -> dict[str, list[Point]]:
    """Return aluminium's and gold's n, k and epsilon, and the workbook's gold."""
    aluminium = _book("Al", ALUMINIUM, "supplement_02")
    gold = _book("Au", GOLD, "supplement_03")
    return _points("Al", aluminium) | _points("Au", gold) | _sheet()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Babar & Weaver (2015): Cu, Ag and Au, from Table 1 and the website file.

Table 1 (``tbl_8cb2a14ef702``, p. 479, 0.10-2.00 eV) and its continuation
(``tbl_76629faf19ce``, p. 480, 2.10-6.00 eV) print, against the photon energy E (eV),
three columns per metal: the reflectance R, epsilon_1 and epsilon_2. Both parts were
checked against their images: the header (a merged "Cu | Ag | Au" row over "E (eV) | R |
epsilon_1 | epsilon_2 ...") is the library's first two CSV rows, and every cell agrees.
R has no quantity: it stays in ``printed``, with the row's other cells. The printed R
agrees with |(sqrt(eps) - 1)/(sqrt(eps) + 1)|^2 from the same row to at most 0.0013
(Ag at 3.4 eV, where epsilon_2 = 0.25 is printed to two decimals), a check on the
transcription (2026-09-27).

The supplement (``original/source_02``, 17 pages) holds figures only (n, k, epsilon and
R against E to 6 and 30 eV); nothing is extracted from it.

The website data: the paper's "Weaver website" file, ``optical-constantsCuAgAu.xlsx``
(attached to the paper in the library, with its URL and retrieval date), sheet
``This work``, read with :mod:`_xlsx`. Row 1 names the metals over three blocks of six
columns (a merged cell each), row 2 the columns ``E | R | n | k | e1 | e2``; each metal
has its own energy column (the grids agree to 6.10 eV, then Cu and Ag run 6.42, 6.92
... 29.92 eV, and Ag to 30.42, while Au runs 6.51, 7.01 ... 30.01 eV). Cu has 169 rows,
Ag 170 and Au 168. The numbers are the cells' stored values in their shortest exact
form (the sheet displays them rounded to 0.00 or 0.000; that is presentation). n and k
are extracted as the primary data and epsilon_1, epsilon_2 as checks; R stays in
``printed`` with the row's other cells. Checked (2026-09-28): in every row of the
three metals, epsilon_1 = n^2 - k^2 and epsilon_2 = 2nk hold within the bound set by
the five printed digits of each number (worst: Cu at 10.92 eV, at 0.96 of the bound),
and R = |(N - 1)/(N + 1)|^2 to 2.3e-5 (Ag, 5.95 eV); at the 69 energies of Table 1
(0.10-6.00 eV) the sheet's R, epsilon_1 and epsilon_2 round to Table 1's printed values
in every row of every metal.
"""

from _library import Point, attachment, number, table, write
from _xlsx import sheet

DIRECTORY = "Babar_2015_OpticalConstantsCu"
PARTS = ("tbl_8cb2a14ef702", "tbl_76629faf19ce")

# The metal's first column (R) in each row; epsilon_1 and epsilon_2 follow it.
METALS = {"Cu": 1, "Ag": 4, "Au": 7}
HEADER = ["E (eV)", "R", "ε1", "ε2", "R", "ε1", "ε2", "R", "ε1", "ε2"]

WORKBOOK_PAPER = "Babar_2015_OpticalConstantsCu"
WORKBOOK = "optical-constantsCuAgAu.xlsx"
# The metal's first column (its own E) in the sheet "This work".
WEB_METALS = {"Cu": 0, "Ag": 6, "Au": 12}
WEB_HEADER = ["E", "R", "n", "k", "e1", "e2"]
# Dataset suffix -> offset from the metal's E column.
WEB_COLUMNS = {"n": 2, "k": 3, "eps1": 4, "eps2": 5}
WEB_ROWS = {"Cu": 169, "Ag": 170, "Au": 168}


def extract() -> dict[str, list[Point]]:
    """Return Table 1's epsilon and the website's n, k and epsilon of each metal."""
    files: dict[str, list[Point]] = {
        f"{metal}/babar_2015__{part}_{metal}.csv": []
        for metal in METALS
        for part in ("eps1", "eps2")
    }
    for table_id in PARTS:
        rows = table(DIRECTORY, table_id)
        if rows[0][2::3] != ["Cu", "Ag", "Au"] or rows[1] != HEADER:
            msg = f"{table_id}: unexpected header {rows[:2]}"
            raise ValueError(msg)
        for cells in rows[2:]:
            energy = number(cells[0])
            printed = " | ".join(cells)
            for metal, first in METALS.items():
                eps1, eps2 = number(cells[first + 1]), number(cells[first + 2])
                files[f"{metal}/babar_2015__eps1_{metal}.csv"].append(
                    Point(energy, energy, eps1, "", printed)
                )
                files[f"{metal}/babar_2015__eps2_{metal}.csv"].append(
                    Point(energy, energy, eps2, "", printed)
                )
    return files | _website()


def _website() -> dict[str, list[Point]]:
    """Return n, k, epsilon_1 and epsilon_2 of each metal from the website's sheet."""
    rows = sheet(attachment(WORKBOOK_PAPER, WORKBOOK), "This work")
    if [rows[0][i] for i in WEB_METALS.values()] != list(WEB_METALS):
        msg = f"This work: unexpected metals row {rows[0]}"
        raise ValueError(msg)
    if rows[1] != WEB_HEADER * 3:
        msg = f"This work: unexpected header {rows[1]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for metal, first in WEB_METALS.items():
        block = [cells[first : first + 6] for cells in rows[2:]]
        block = [cells for cells in block if cells and cells[0]]
        if len(block) != WEB_ROWS[metal] or any(len(c) != 6 for c in block):
            msg = f"This work: {metal} has {len(block)} rows"
            raise ValueError(msg)
        for name, offset in WEB_COLUMNS.items():
            files[f"{metal}/babar_2015__web_{name}_{metal}.csv"] = [
                Point(c[0], number(c[0]), number(c[offset]), "", " | ".join(c))
                for c in block
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

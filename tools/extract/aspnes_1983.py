# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Aspnes & Studna (1983): the optical properties of Si, Ge, GaAs, InP and GaSb.

Tables II-IX (one per material, pp. 6-21 of the pdf) print, against the photon energy
(eV, 1.500-6.000 in 0.100 eV steps, 46 rows), the pseudodielectric function
<epsilon_1> and <epsilon_2> measured by spectroscopic ellipsometry, and n, k, the
normal-incidence reflectance R and the absorption coefficient 10^-3 alpha (cm^-1)
calculated from it. n and k are extracted as the primary data, epsilon_1 and
epsilon_2 as checks; R and alpha stay in ``printed`` with the row's other cells.

Only the materials registered in optmat are extracted: Si (Table II), Ge (III), GaAs
(V), GaSb (VI) and InP (VII). GaP (IV), InAs (VIII) and InSb (IX) are not registered.

Every table was compared with its image (Si cell by cell; the others through the
redundant columns): in every row of the five tables epsilon = N^2 holds within the
printed rounding, R from N within its rounding and alpha = 4 pi k E / (h c) within
the rounding of epsilon_2, except one cell of Si, 3.100 eV <epsilon_2>, which the
image prints "4.321" (as read); the record corrects it to 4.312 by an erratum.
"""

from _library import Point, number, table, write

DIRECTORY = "Aspnes_1983_DielectricFunctionsOptical"
HEADER = [
    "E (eV)",
    r"\langle \epsilon_1 \rangle",
    r"\langle \epsilon_2 \rangle",
    "n",
    "k",
    "R",
    r"10^3\alpha (cm-1)",
]

# Material -> library table (the paper's table number in the comment).
TABLES = {
    "Si": "tbl_f343b1b2f612",  # Table II
    "Ge": "tbl_24578d7291a9",  # Table III
    "GaAs": "tbl_fdeda94e867c",  # Table V
    "GaSb": "tbl_e0d0f2a98631",  # Table VI
    "InP": "tbl_f0b0b06653cf",  # Table VII
}

# Dataset prefix -> column of the tables.
COLUMNS = {"n": 3, "k": 4, "eps1": 1, "eps2": 2}

ENERGIES = [f"{0.1 * step:.3f}" for step in range(15, 61)]


def _rows(material: str) -> list[list[str]]:
    """Return a material's table rows, checking the header and the energy grid."""
    rows = table(DIRECTORY, TABLES[material])
    if rows[0] != HEADER:
        msg = f"{material}: unexpected header {rows[0]}"
        raise ValueError(msg)
    body = [[cell.strip() for cell in row] for row in rows[1:]]
    if [row[0] for row in body] != ENERGIES or any(len(row) != 7 for row in body):
        msg = f"{material}: the rows are not 1.500-6.000 eV in 0.1 eV steps"
        raise ValueError(msg)
    return body


def extract() -> dict[str, list[Point]]:
    """Return n, k, epsilon_1 and epsilon_2 of the five materials."""
    files = {}
    for material in TABLES:
        rows = _rows(material)
        for name, column in COLUMNS.items():
            files[f"{material}/aspnes_1983__{name}_{material}.csv"] = [
                Point(row[0], row[0], number(row[column]), "", " | ".join(row))
                for row in rows
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Green (2008): alpha, n and k of intrinsic silicon at 300 K (Table 1).

Table 1 ("Optical properties of intrinsic silicon at 300 K including normalised
temperature coefficients", pp. 1306-1307) runs over three library tables:
``tbl_85bbdac968fa`` (p. 2 of the PDF, 0.25-0.82 um), ``tbl_282abc185c1f`` (p. 3,
0.83-1.41 um, under a "Table 1 (continued)" line) and ``tbl_76272313742f`` (p. 4,
1.42-1.45 um): 121 rows at 0.01 um steps. The columns are the free-space wavelength
lambda (um), alpha (/cm), n, k, C_n (10^-4/K, (1/n)(dn/dT)) and C_k,alpha (10^-4/K,
(1/k)(dk/dT), which the note says alpha shares). The library renders the Greek
letters and the micro sign as Latin ("l (mm)", "a (/cm)"); the images print
"lambda (um)" and "alpha (/cm)".

Each part was checked against its image, and every cell against the PDF's text layer
(``pdftotext -layout``), on 2026-10-14: they agree. k is printed only to 1.00 um; from
1.01 um the image prints an en dash ("-" in the CSV), and those rows are left out of
the k file (alpha and n run to 1.45 um). The temperature coefficients have no quantity
(they are normalised, d ln p / dT); they stay in each row's ``printed`` cells.
"""

from _library import Point, number, table, write

DIRECTORY = "Green_2008_SelfConsistentOptical"
PARTS = ("tbl_85bbdac968fa", "tbl_282abc185c1f", "tbl_76272313742f")
HEADER = [
    "l (mm)",
    "a (/cm)",
    "n",
    "k",
    "Cn (10-4/K) (1/n) (dn/dt)",
    "Ck,a (10-4/K) (1/k)(dk/dt)",
]
CONTINUED = ["Table 1 (continued)", "", "", "", "", ""]
ROWS = 121
PREFIX = "Si/green_2008__"


def rows() -> list[list[str]]:
    """Return Table 1's data rows from its three parts, headers checked and removed."""
    found: list[list[str]] = []
    for table_id in PARTS:
        cells = table(DIRECTORY, table_id)
        if cells[0] == CONTINUED:
            cells = cells[1:]
        if cells[0] != HEADER:
            msg = f"{table_id}: unexpected header {cells[0]}"
            raise ValueError(msg)
        found += cells[1:]
    wavelengths = [number(c[0]) for c in found]
    expected = [f"{i / 100:g}" for i in range(25, 146)]
    if len(found) != ROWS or wavelengths != expected:
        msg = f"Table 1: {len(found)} rows, not the 0.25-1.45 um grid"
        raise ValueError(msg)
    return found


def extract() -> dict[str, list[Point]]:
    """Return alpha, n and k at 300 K."""
    alpha, n, k = [], [], []
    for cells in rows():
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        alpha.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        n.append(Point(wavelength, wavelength, number(cells[2]), "", printed))
        if cells[3] != "-":
            k.append(Point(wavelength, wavelength, number(cells[3]), "", printed))
    return {
        f"{PREFIX}alpha_300K.csv": alpha,
        f"{PREFIX}n_300K.csv": n,
        f"{PREFIX}k_300K.csv": k,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

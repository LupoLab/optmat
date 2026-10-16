# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ordal et al. (1985): n, k, -epsilon_1 and epsilon_2 of Cu, 180-19400 cm-1.

Table II (``tbl_08c9b2542049``, p. 4495, PDF page 3) prints omega (cm-1), lambda (um),
-epsilon_1, epsilon_2, n and k for 56 wavenumbers: the authors' own Kramers-Kronig
analysis of the reflectance of copper. The library's cells agree with the image; they
keep the typeset spacing of the exponents ("1.80E + 2"), which the script removes to
store the number ("1.80E+2"), the cells staying as printed in ``printed``. The fourth
column is -epsilon_1: stored as printed, the dataset's scale of -1 makes it epsilon_1.

The library's note on this table, "Table III. Nickel", is the caption of the table
beside it on the page, which the library did not extract (Ni is not transcribed).
Every row satisfies epsilon_1 = n^2 - k^2, epsilon_2 = 2nk and lambda = 1e4 / omega to
the printed rounding (worst: epsilon_2 at 1900 cm-1, 239 printed, 238.05 from n and k;
checked 2026-09-28).
"""

from _library import Point, number, table, write

DIRECTORY = "Ordal_1985_OpticalPropertiesFourteen"
TABLE_II = "tbl_08c9b2542049"
HEADER = [
    "\\omega(cm^{-1})",
    "\\lambda(\\mu m)",
    "-\\epsilon_{1}",
    "\\epsilon_{2}",
    "n",
    "k",
]
# Dataset prefix -> column (the third column is -epsilon_1).
COLUMNS = {"n": 4, "k": 5, "eps1": 2, "eps2": 3}
ROWS = 56


def _number(cell: str) -> str:
    """Return the number of a cell printed with spaced exponents ("1.80E + 2")."""
    return number(cell.replace(" ", ""))


def extract() -> dict[str, list[Point]]:
    """Return copper's n, k, -epsilon_1 and epsilon_2 against the wavenumber."""
    rows = table(DIRECTORY, TABLE_II)
    if rows[0] != HEADER or len(rows) - 1 != ROWS:
        msg = f"{TABLE_II}: unexpected header {rows[0]} or {len(rows) - 1} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"Cu/ordal_1985__{name}_Cu.csv": [] for name in COLUMNS
    }
    for cells in rows[1:]:
        wavenumber = _number(cells[0])
        printed = " | ".join(cells)
        for name, column in COLUMNS.items():
            files[f"Cu/ordal_1985__{name}_Cu.csv"].append(
                Point(wavenumber, wavenumber, _number(cells[column]), "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

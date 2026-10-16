# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ripoche et al. (1997), Opt. Commun. 135, 310: the air parameters of Table 1.

Table 1 (``tbl_d9854060f4cb``, CSV; its caption, "Values of the non-linear refractive
index parameters for O2 and N2 [13]", is lost in the extraction; checked against the
page image, p. 313): per gas, the concentration in air (not written), the rotational
constant B (cm-1) and the polarizability anisotropy beta(omega) (C^2 m^2/J, which is C
m^2/V, the SI polarizability), all from Ref. [13] (Lin et al. 1976). The cell
"1.27 \times 10^{-40}" is written as 1.27 with the dataset scale 1e-40. No wavelength is
stated for beta: the abscissa cell is "not stated".
"""

from _library import Point, number, table, write

DIRECTORY = "Ripoche_1997_DeterminationTimeDependence"
TABLE_1 = "tbl_d9854060f4cb"
GASES = {"O_{2}": "O2", "N_{2}": "N2"}


def extract() -> dict[str, list[Point]]:
    """Return one B file and one anisotropy file per gas."""
    grid = table(DIRECTORY, TABLE_1)
    header = [
        "Gas",
        "Concentration in air",
        "B (cm^{-1})",
        "\\beta(\\omega) (C^{2} m^{2}/J)",
    ]
    if grid[0] != header:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        files[f"{gas}/ripoche_1997__table1_B_{gas}.csv"] = [
            Point("B", "B", number(cells[2]), "", printed)
        ]
        files[f"{gas}/ripoche_1997__table1_anisotropy_{gas}.csv"] = [
            Point(gas, "not stated", number(cells[3]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

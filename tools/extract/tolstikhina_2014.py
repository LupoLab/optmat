# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tolstikhina, Morishita & Tolstikhin (2014): the coefficients of Table I.

Table I (``tbl_bb65aa194b35``, CSV; image checked, p. 053413-4) prints, per atom with
not more than one electron above closed shells, the term, the ionising orbital, the
Hartree-Fock orbital exponent kappa (column headed "x" in the CSV), the occupation
factor p, the relaxation coefficient R_N, the orbital's asymptotic coefficient g_N
(Eq. 12) and the coefficient g = sqrt(p) R_N g_N of the rate Gamma = |g|^2 W(F) (Eqs. 1,
17). Only the curated species are taken: He, Ne, Ar, Kr and Xe (H's row, sqrt(2) by
definition, and Li, Be, B, Na, Mg, Al are not taken). kappa is printed with its digits
grouped by spaces ("1.354 958"); the grouping is removed.
"""

from _library import Point, number, table, write

DIRECTORY = "Tolstikhina_2014_ApplicationManyElectron"
TABLE_I = "tbl_bb65aa194b35"
HEADER = ["Atom", "^{2S+1}L", "n_Nl_N", "x", "p", "R_N", "g_N", "g"]
ATOMS = {"He": "1s", "Ne": "2p", "Ar": "3p", "Kr": "4p", "Xe": "5p"}
ORBITAL_L = {"s": "0", "p": "1"}


def _grouped(cell: str) -> str:
    """Remove the digit grouping of a printed number."""
    return number(cell.replace(" ", "").replace("\N{THIN SPACE}", ""))


def extract() -> dict[str, list[Point]]:
    """Return kappa, R_N, g_N and g of He, Ne, Ar, Kr and Xe."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    rows = {cells[0].strip(): cells for cells in grid[1:]}
    files: dict[str, list[Point]] = {}
    for gas, orbital in ATOMS.items():
        cells = rows[gas]
        if cells[2].strip() != orbital or cells[4].strip() != "2":
            msg = f"unexpected Table I row {cells!r}"
            raise ValueError(msg)
        printed = " | ".join(cells)
        wave = ORBITAL_L[orbital[-1]]
        stem = f"{gas}/tolstikhina_2014__"
        files[f"{stem}kappa_{gas}.csv"] = [
            Point("kappa", "", _grouped(cells[3]), "", printed)
        ]
        files[f"{stem}relaxation_{gas}.csv"] = [
            Point("R_N", "", number(cells[5]), "", printed)
        ]
        files[f"{stem}g_orbital_{gas}.csv"] = [
            Point(wave, wave, number(cells[6]), "", printed)
        ]
        files[f"{stem}g_{gas}.csv"] = [Point(wave, wave, number(cells[7]), "", printed)]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

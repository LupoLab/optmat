# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Loriot et al. (2011): the ionisation potentials of Table II.

Table II (``tbl_701638991ad8``, CSV; image checked, preprint p. 14): for N2, O2 and Ar,
IP_neutral and IP_ion (eV) and the ratio chi(3)_ion/chi(3)_neutral of Eq. (6). The
neutral's value is written as the molecule's or atom's ionisation energy
(``ip_<gas>``), the ion's as the cation's (``ip_<gas>+``). The ratio row, an
estimate from Eq. (6), is not written.
"""

from _library import Point, number, table, write

DIRECTORY = "Loriot_2010_NegativeHigherOrder"
TABLE = "tbl_701638991ad8"
SPECIES = {"N_{2}": "N2", "O_{2}": "O2", "Ar": "Ar"}
ROWS = {"IP_{neutral}(eV)": "", "IP_{ion}(eV)": "+"}


def extract() -> dict[str, list[Point]]:
    """Return the ionisation energies of N2, O2, Ar and their cations."""
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE)
    if grid[0][0] != "Species":
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    species = [SPECIES[cell.strip()] for cell in grid[0][1:]]
    for cells in grid[1:]:
        suffix = ROWS.get(cells[0].strip())
        if suffix is None:  # the chi(3) ratio row
            continue
        for gas, cell in zip(species, cells[1:], strict=True):
            material = gas + suffix
            printed = f"{gas} | {cells[0]} | {cell}"
            files[f"{material}/loriot_2011__ip_{material}.csv"] = [
                Point("IP", "", number(cell), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

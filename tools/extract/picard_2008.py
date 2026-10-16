# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Picard, Davis, Gläser & Fujii (2008): the CIPM-2007 composition of dry air.

Table 1 (``tbl_b2ed47a88f5d``, CSV; image checked digit by digit, no differences) gives
per constituent its molar mass M_i (10^-3 kg/mol), its mole fraction x_i and the
contribution x_i M_i. Only the mole fraction is extracted, one file per constituent
(material = the constituent; all twelve are registered gases), row = the constituent as
printed, abscissa none. Two structural fixes, both from the image: the digits of the
major fractions are printed in groups of three separated by a thin space ("0.780 848"),
removed here; the minor ones are printed "18.2 \times 10^{-6}", whose factor the helper
strips (their datasets carry scale = 1e-6).
"""

from _library import Point, number, table, write

DIRECTORY = "Picard_2008_RevisedFormulaDensity"
TABLE_1 = "tbl_b2ed47a88f5d"
MATERIALS = {
    "N_2": "N2",
    "O_2": "O2",
    "Ar": "Ar",
    "CO_2": "CO2",
    "Ne": "Ne",
    "He": "He",
    "CH_4": "CH4",
    "Kr": "Kr",
    "H_2": "H2",
    "N_2O": "N2O",
    "CO": "CO",
    "Xe": "Xe",
}


def extract() -> dict[str, list[Point]]:
    """Return Table 1's mole fractions, one file per constituent."""
    grid = table(DIRECTORY, TABLE_1)
    if grid[0][0] != "Constituent" or "fraction" not in grid[0][2]:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        material = MATERIALS[cells[0]]
        cell = cells[2]
        if "\\times" not in cell:
            cell = cell.replace(" ", "")  # thin-space digit groups (image)
        point = Point(material, "", number(cell), "", " | ".join(cells))
        files[f"{material}/picard_2008__x_{material}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Werner et al. (2009): epsilon of Cu, Ag and Au from DFT and from REELS.

Tables 11-26 give, for each metal but Te, epsilon_1, epsilon_2 and the energy-loss
function (ELF = Im(-1/epsilon)) from DFT and from REELS on a common grid of energies
(0.50-4.75 eV in 0.25 eV steps, 5.0-70.5 eV in 0.5 eV steps: 150 rows), with the REELS
single-scattering distributions w_b (DIIMFP, meV^-1) and w_s (DSEP, arbitrary units)
in the last two columns. Each table runs over three pages. The library publishes these
tables as JSON and HTML only (merged header cells), so they are read from the JSON cells
(``_json_table``). The parts and their order, established from the images (each part is
one page, headed with the metal and "DFT | REELS", and its energies continue the
previous page's):

- Cu, Table 16: ``tbl_bbde835d361c`` (p. 1059, 0.50-21.5 eV), ``tbl_77acafd389e3``
  (p. 1060, 22.0-47.5 eV), ``tbl_d3f670af9b86`` (p. 1061, 48.0-70.5 eV);
- Ag, Table 20: ``tbl_a521b89d0a14`` (p. 1071), ``tbl_ba923cb150a2`` (p. 1072),
  ``tbl_d7484476dcbe`` (p. 1073);
- Au, Table 24: ``tbl_991d52c88cc2`` (p. 1083), ``tbl_f4fe6ec8cdaa`` (p. 1084),
  ``tbl_4c92ff38fd6e`` (p. 1085).

The first parts' rows below 6 eV and the last parts of Au were compared with the images
cell by cell; every row of every part satisfies ELF = epsilon_2 / (epsilon_1^2 +
epsilon_2^2) to the printed rounding (at most 0.0015 off, checked 2026-09-27), for DFT
and REELS alike, which ties each epsilon pair to the independently printed ELF. The DFT
and REELS energy columns are identical in every row (checked here). The ELF, w_b and w_s
have no quantity: they stay in ``printed`` with the rest of the row. w_s is printed "-"
above 49.5 eV.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Werner_2009_OpticalConstantsInelastic"
PARTS = {
    "Cu": ("tbl_bbde835d361c", "tbl_77acafd389e3", "tbl_d3f670af9b86"),
    "Ag": ("tbl_a521b89d0a14", "tbl_ba923cb150a2", "tbl_d7484476dcbe"),
    "Au": ("tbl_991d52c88cc2", "tbl_f4fe6ec8cdaa", "tbl_4c92ff38fd6e"),
}
HEADER = [
    "E(eV)",
    "\\varepsilon_1",
    "\\varepsilon_2",
    "ELF",
    "E(eV)",
    "\\varepsilon_1",
    "\\varepsilon_2",
    "ELF",
    "w_b(meV^{-1})",
    "w_s(arb. un.)",
]
# (method, column of epsilon_1); epsilon_2 is the next column.
METHODS = (("DFT", 1), ("REELS", 5))
ROWS = 150


def extract() -> dict[str, list[Point]]:
    """Return epsilon_1 and epsilon_2 of each metal by DFT and by REELS."""
    files: dict[str, list[Point]] = {}
    for metal, parts in PARTS.items():
        stem = f"{metal}/werner_2009__"
        for method, _ in METHODS:
            for part in ("eps1", "eps2"):
                files[f"{stem}{part}_{method.lower()}_{metal}.csv"] = []
        previous = -1.0
        for table_id in parts:
            rows = json_table(DIRECTORY, table_id)
            if rows[0][0] != metal or rows[1][0] != "DFT" or rows[2] != HEADER:
                msg = f"{table_id}: unexpected header {rows[:3]}"
                raise ValueError(msg)
            for cells in rows[3:]:
                energy = number(cells[0])
                if number(cells[4]) != energy or float(energy) <= previous:
                    msg = f"{table_id}: energies out of step at {cells}"
                    raise ValueError(msg)
                previous = float(energy)
                printed = " | ".join(cells)
                for method, column in METHODS:
                    for part, offset in (("eps1", 0), ("eps2", 1)):
                        value = number(cells[column + offset])
                        files[f"{stem}{part}_{method.lower()}_{metal}.csv"].append(
                            Point(energy, energy, value, "", printed)
                        )
        counts = {len(points) for key, points in files.items() if key.startswith(stem)}
        if counts != {ROWS}:
            msg = f"{metal}: {counts} rows, not {ROWS}"
            raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

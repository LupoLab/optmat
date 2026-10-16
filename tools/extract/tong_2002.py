# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tong, Zhao & Lin (2002): the MO-ADK structure parameters (Tables I and II).

Table I (``tbl_3048e91bd1cc``, cells file only; image checked, p. 3) prints, per
diatomic valence orbital, the ionisation energy I_p (eV), the internuclear distance R
(angstrom) and the coefficients C_l of Eq. (1). It has two blocks: the homonuclear
molecules with columns l = 0, 2, 4 (the cells file puts "l=0" in the first header row
and "l=2", "l=4" in the second, under the spanning "C_l"), and a second header row
"I_p (eV) | R | l=0 | l=1 | l=2 | l=3" for the heteronuclear ones. The script reads
the partial waves from the header in force. A blank C_l cell (l = 0 of a pi orbital,
which has no l = 0 wave) is not a value and is skipped. Only the curated species are
taken: H2+ (an ion), the hydrogen row printed "D_2", N2, O2 and CO; F2, S2, NO and SO
are not curated (see the record). R is written as the per-row column
``internuclear_distance`` of each C_l file.

Table II (``tbl_58844fa6da0d``, CSV; image checked) prints, per rare-gas atom, I_p
(eV), this paper's C_l (self-interaction-free density-functional theory) and the ADK
coefficient C_l(ADK), in columns headed "He (1s)" ... "Xe (5p)"; the single partial
wave is the orbital's l (s: 0, p: 1).
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Tong_2002_TheoryMolecularTunneling"
TABLE_I = "tbl_3048e91bd1cc"
TABLE_II = "tbl_58844fa6da0d"
SIGMA = "\N{GREEK SMALL LETTER SIGMA}"
# Table I row label -> curated material.
MOLECULES = {
    f"H_2^+({SIGMA}_g)": "H2+",
    "D_2(\\sigma_g)": "D2",
    "N_2(\\sigma_g)": "N2",
    "O_2(\\pi_g)": "O2",
    f"CO ({SIGMA})": "CO",
}
SKIPPED = {"F_2(\\pi_g)", "S_2(\\pi_g)", "NO (π)", "SO (π)"}
_WAVE = re.compile(r"^l=(\d)$")
_ATOM = re.compile(r"^(He|Ne|Ar|Kr|Xe) \((\d)([sp])\)$")
ORBITAL_L = {"s": "0", "p": "1"}


def _table_i() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_I)
    # The first block's partial waves are split over its two header rows.
    waves = {3: "0", 4: "2", 5: "4"}
    if [grid[0][3], grid[1][4], grid[1][5]] != ["l=0", "l=2", "l=4"]:
        msg = f"unexpected Table I header {grid[:2]!r}"
        raise ValueError(msg)
    for cells in grid[2:]:
        label = cells[0].strip()
        if not label:  # the second block's header row
            waves = {
                column: match.group(1)
                for column, cell in enumerate(cells)
                if (match := _WAVE.match(cell.strip()))
            }
            continue
        if label in SKIPPED:
            continue
        gas = MOLECULES[label]
        printed = " | ".join(cells)
        files[f"{gas}/tong_2002__ip_{gas}.csv"] = [
            Point("I_p", "", number(cells[1]), "", printed)
        ]
        files[f"{gas}/tong_2002__cl_{gas}.csv"] = [
            Point(
                wave,
                wave,
                number(cells[column]),
                "",
                printed,
                {"internuclear_distance": number(cells[2])},
            )
            for column, wave in waves.items()
            if cells[column].strip()
        ]
    return files


def _table_ii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE_II)
    rows = {cells[0].strip(): cells for cells in grid[1:]}
    if set(rows) != {"Ip ~eV!", "Cl", "Cl~ADK!"}:
        msg = f"unexpected Table II rows {sorted(rows)}"
        raise ValueError(msg)
    for column, heading in enumerate(grid[0][1:], start=1):
        match = _ATOM.match(heading.strip())
        if match is None:
            msg = f"unexpected Table II column {heading!r}"
            raise ValueError(msg)
        gas, wave = match.group(1), ORBITAL_L[match.group(3)]
        printed = " | ".join(
            f"{name}: {rows[name][column]}" for name in ("Ip ~eV!", "Cl", "Cl~ADK!")
        )
        printed = f"{heading} | {printed}"
        files[f"{gas}/tong_2002__ip_{gas}.csv"] = [
            Point("I_p", "", number(rows["Ip ~eV!"][column]), "", printed)
        ]
        files[f"{gas}/tong_2002__cl_{gas}.csv"] = [
            Point(wave, wave, number(rows["Cl"][column]), "", printed)
        ]
        files[f"{gas}/tong_2002__cl_adk_{gas}.csv"] = [
            Point(wave, wave, number(rows["Cl~ADK!"][column]), "", printed)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the ionisation energies and C_l of Tables I and II."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

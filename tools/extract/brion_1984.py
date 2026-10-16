# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Brion & Thomson (1984): total photoabsorption and photoionisation df/dE.

The compilation tabulates, per molecule, dipole (e,e) total photoabsorption and dipole
(e,e+ion) total photoionisation oscillator strengths df/dE in 10^-2 eV^-1 against
photon energy (eV), with the partial (ionic-fragment) channels. The "(e,e+ion)" tables
of the curated species are read (cells files only; every table image checked, all
cells agree): H2 ``tbl_55480943b731`` (p. 318), CO ``tbl_72dab91fceaa`` (p. 319), N2
``tbl_a794d910da1a`` (p. 320), O2 ``tbl_00a15c604fe1`` (p. 321), CH4
``tbl_689c845cc210`` (p. 326), N2O ``tbl_c5ae487e5662`` (p. 327), CO2
``tbl_f55f568f42e8`` (p. 328), SF6 ``tbl_643f8c4f78a2`` (p. 331). Columns 0-2 are the
energy, the total photoabsorption and the total photoionisation; the partial channels
(columns 3 on) are not extracted. Handled structurally:

- the title strips ("CO (e,e+ion) J. Phys. B 9 (1976) 675-689.") are extracted as
  separate one-cell tables, the SF6 one reading "SF_5" for SF6, and the CH4 (e,e+ion)
  strip of p. 326 is missing from the extraction (the page image shows it; the CH4
  strip the extraction has, ``tbl_e7bbffed7df0``, belongs to the (e,2e) table of p.
  313): the gas of each table is set here, from the page;
- the three header rows are skipped; rows with no photoionisation value (below the
  ionisation threshold, or beyond the ion data) give an absorption point only.

The photoionisation values are written in the printed unit; the source record turns
them into a cross section with the paper's own factor, sigma(Mb) = 109.75 df/dE
(eV^-1) (p. 303), as the dataset scale 1.0975.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Brion_1984_CompilationDipoleOscillator"
TABLES = {
    "H2": "tbl_55480943b731",
    "CO": "tbl_72dab91fceaa",
    "N2": "tbl_a794d910da1a",
    "O2": "tbl_00a15c604fe1",
    "CH4": "tbl_689c845cc210",
    "N2O": "tbl_c5ae487e5662",
    "CO2": "tbl_f55f568f42e8",
    "SF6": "tbl_643f8c4f78a2",
}
HEADER_ROWS = 3


def extract() -> dict[str, list[Point]]:
    """Return every data file of the source, keyed by its ``<gas>/<file>`` name."""
    files: dict[str, list[Point]] = {}
    for gas, table_id in TABLES.items():
        absorption: list[Point] = []
        ionization: list[Point] = []
        for cells in json_table(DIRECTORY, table_id)[HEADER_ROWS:]:
            energy, total, ionized = cells[0], cells[1], cells[2]
            printed = f"{energy} | {total} | {ionized}"
            energy = number(energy)
            absorption.append(Point(energy, energy, number(total), "", printed))
            if ionized:
                ionization.append(Point(energy, energy, number(ionized), "", printed))
        files[f"{gas}/brion_1984__absorption_{gas}.csv"] = absorption
        files[f"{gas}/brion_1984__ionization_{gas}.csv"] = ionization
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

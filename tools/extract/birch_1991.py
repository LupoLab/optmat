# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Birch (1991): molar refractions and refractivities of gases at 633 nm.

- Table 1 (``tbl_6072c809ac8e``, cells file only): per gas, the runs (temperature,
  pressure, measured (n - 1) x 1e8, real-gas density, molar refraction [R] in 1e-30 m^3
  per molecule) and the mean molar refraction and polarizability alpha-bar
  = 3 [R] / (4 pi) (1e-30 m^3, a polarizability volume) printed once per gas on the
  middle row. Per gas: ``alpha_<gas>`` holds alpha-bar; ``runs_<gas>`` the per-run [R]
  (rows numbered 1, 2, ... within the gas; the run's T, p, n - 1 and density are in
  the printed cells, as the format has no per-row state). Gas labels carry the purity
  ("Nitrogen(99.998%)"); the gas name is taken before the parenthesis. Digit groups
  printed with spaces ("101 524.2", and "10 1188.6" for the first CO2 run, a misplaced
  space in the print, image checked) stay in the printed cells only.
- Table 2 (``tbl_10576f9abc2b``, cells file only): 10^6 (n - 1) at 632.99 nm,
  101 325 Pa and 15 degrees C, "value +- uncertainty"; the first row of each gas is
  Birch's (derived from [R]), the others are reproduced from the reference numbers in
  the last column.

Air (Table 3 and the composition table) is not a curated species. The images agree
with the cells.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Birch_1991_PreciseDeterminationRefractometric"
TABLE_1 = "tbl_6072c809ac8e"
TABLE_2 = "tbl_10576f9abc2b"
GASES = {
    "Nitrogen": "N2",
    "Oxygen": "O2",
    "Argon": "Ar",
    "CarbonDioxide": "CO2",
    "Carbon dioxide": "CO2",
    "Neon": "Ne",
    "Helium": "He",
    "NitrousOxide": "N2O",
    "Nitrous oxide": "N2O",
    "Krypton": "Kr",
}
# Table 2 reference number -> dataset name.
REFERENCES = {
    "21": "kerl_1982",
    "16": "peck_1966",
    "14": "cuthbertson_1910",
    "17": "peck_1964",
    "9": "edlen_1966",
    "20": "simmons_1978",
    "15": "cuthbertson_1932",
    "18": "mansfield_peck_1969",
    "19": "hey_1988",
}
WAVELENGTH = "632.99"


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas = ""
    run = 0
    for cells in json_table(DIRECTORY, TABLE_1)[1:]:
        if cells[0]:
            gas = GASES[cells[0].split("(")[0].strip()]
            run = 0
        run += 1
        printed = " | ".join(cells)
        files.setdefault(f"{gas}/birch_1991__runs_{gas}.csv", []).append(
            Point(str(run), WAVELENGTH, number(cells[5]), "", printed)
        )
        if cells[7]:
            files[f"{gas}/birch_1991__alpha_{gas}.csv"] = [
                Point(WAVELENGTH, WAVELENGTH, number(cells[7]), "", printed)
            ]
    return files


def _table_2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_2)[1:]:
        if cells[0]:
            gas = GASES[cells[0]]
        value, _, sigma = cells[1].partition("\N{PLUS-MINUS SIGN}")
        name = REFERENCES[cells[2]] if cells[2] else "refractivity"
        point = Point(
            WAVELENGTH, WAVELENGTH, number(value), sigma.strip(), " | ".join(cells)
        )
        path = f"{gas}/birch_1991__{name}_{gas}.csv"
        assert path not in files, path
        files[path] = [point]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table 1 (alpha-bar and the runs) and Table 2."""
    return _table_1() | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Orcutt & Cole (1967): dielectric virial coefficients (He-Kr, H2, N2).

- Table II (``tbl_f6b30fa838e4``, cells file only): the final A_eps (cm^3/mol) of this
  work, Newell and Baird's and the NBS (Maryott and Buckley) values, and B_eps
  (cm^6/mol^2, "(corr)": corrected for the pressure dependence of the cells'
  geometric capacitance) of this work, the literature values (marks c, d, e) and the
  calculated fluctuation-induced contribution. One file per gas and column. The
  argon literature cell holds two values, "1.0+-0.7^{c}2.2+-0.2^{d}" (image: two
  lines, Refs. 1 and 5), split by their marks. The quadrupole-induced column (a
  partial contribution) is not extracted.
- Table III (``tbl_3e1a5c368033``): the authors' fits of others' refractive-index
  data; only the CO2 row at the single wavelength 5461 A (Phillips 1920, mark c,
  printed "6" in the library's footnote) is extracted; the rows averaged over
  4471-6678 A have no single wavelength and are left to the source notes.

Table I (per-run values before the capacitance correction; its cells merge the
stacked runs) is not extracted. Images checked.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Orcutt_1967_DielectricConstantsImperfect"
TABLE_II = "tbl_f6b30fa838e4"
TABLE_III = "tbl_3e1a5c368033"
GASES = {"He": "He", "Ne": "Ne", "Ar": "Ar", "Kr": "Kr", "H_{2}": "H2", "N_{2}": "N2"}
LITERATURE = {"c": "johnston_1960", "d": "orcutt_cole_1965", "e": "unidentified_e"}
_PM = "\N{PLUS-MINUS SIGN}"
_MARKED = re.compile(r"(-?[\d.]+)" + _PM + r"([\d.]+?)\^\{([a-e])\}")


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``0.519+-0.001``."""
    value, _, sigma = cell.partition(_PM)
    return number(value), sigma.strip()


def _table_ii() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_II)
    files: dict[str, list[Point]] = {}

    def add(gas: str, name: str, point: Point) -> None:
        path = f"{gas}/orcutt_1967__{name}_{gas}.csv"
        assert path not in files, path
        files[path] = [point]

    for cells in rows[1:7]:  # A (cm^3/mol)
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, name in ((1, "A_this_work"), (2, "A_newell_baird"), (4, "A_nbs")):
            if cells[column] != "...":
                add(gas, name, Point("322.15", "inf", *split(cells[column]), printed))
    for cells in rows[9:15]:  # B (cm^6/mol^2)
        gas = GASES[cells[0].removesuffix(" (corr)")]
        printed = " | ".join(cells)
        add(gas, "B_this_work", Point("322.15", "inf", *split(cells[1]), printed))
        for value, sigma, mark in _MARKED.findall(cells[2]):
            name = f"B_{LITERATURE[mark]}"
            add(gas, name, Point("unknown", "inf", number(value), sigma, printed))
        add(gas, "B_fluctuation", Point("322.15", "inf", number(cells[3]), "", printed))
    return files


def _table_iii() -> dict[str, list[Point]]:
    cells = json_table(DIRECTORY, TABLE_III)[3]
    assert cells[:2] == ["CO2, 34°C", "5461"], cells
    printed = " | ".join(cells)
    return {
        "CO2/orcutt_1967__A_n_phillips_1920_CO2.csv": [
            Point("307.15", "5461", *split(cells[2]), printed)
        ],
        "CO2/orcutt_1967__B_n_phillips_1920_CO2.csv": [
            Point("307.15", "5461", *split(cells[3]), printed)
        ],
    }


def extract() -> dict[str, list[Point]]:
    """Return Table II and the single-wavelength row of Table III."""
    return _table_ii() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sneep and Ubachs (2005): Rayleigh cross sections at 532.2 nm.

Table 2 (``tbl_f4af3ec06557``, cells file only; image checked) gives, at the Nd:YAG
second harmonic 18 788.4 cm^-1, the measured extinction cross section with its
one-sigma uncertainty in parentheses (``4.45 (0.3)``, 10^-27 cm^2), the n-based cross
section, the (n - 1) and King factor behind it, the purity and the number of pressure
ramps. The cells lose the subscripts of the gas names (``N2``) and print some decimal
points as colons (``1:035``); the measured and n-based columns are clean. The measured
and n-based cross sections of Ar, N2, CO, CO2, CH4, N2O and SF6 are extracted. The
(n - 1) and F_k columns are recorded in the source record (they check the formulas).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Sneep_2005_DirectMeasurementRayleigh"
TABLE_2 = "tbl_f4af3ec06557"
GASES = ("Ar", "N2", "CO", "CO2", "CH4", "N2O", "SF6")
WAVENUMBER = "18788.4"
_VALUE_SIGMA = re.compile(r"^\s*([\d.]+)\s*\(([\d.]+)\)\s*$")


def extract() -> dict[str, list[Point]]:
    """Return the measured and n-based cross sections of Table 2 as data files."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_2)[2:]:  # a header and an empty row
        gas = cells[0].strip()
        if gas not in GASES:
            continue
        printed = " | ".join(cells)
        match = _VALUE_SIGMA.match(cells[1])
        if match is None:
            msg = f"unexpected measured cell {cells[1]!r}"
            raise ValueError(msg)
        prefix = f"{gas}/sneep_2005__"
        files[f"{prefix}rayleigh_{gas}.csv"] = [
            Point(WAVENUMBER, WAVENUMBER, match.group(1), match.group(2), printed)
        ]
        files[f"{prefix}n_based_{gas}.csv"] = [
            Point(WAVENUMBER, WAVENUMBER, number(cells[2]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

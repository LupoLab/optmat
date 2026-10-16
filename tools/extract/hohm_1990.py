# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hohm & Kerl (1990): noble-gas polarizabilities at 632.99 nm (Table 1).

Table 1 has no CSV rendering (the gas cell spans each group of rows) and is read
from its JSON cells, checked against the image. Values are ``10^30 alpha / m^3``
(polarizability volumes in angstrom^3), printed with a raised decimal point
(``0·20608``) and, for this work and Burns et al., the standard deviation in the last
digits in parentheses (``0·20608 (15)``), which the script turns into an absolute
sigma in the same unit. The reference marks (daggers) are dropped from the author cell.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hohm_1990_InterferometricMeasurementsDipole"
TABLE_1 = "tbl_f8502d3f8f19"
GASES = {"Helium": "He", "Neon": "Ne", "Argon": "Ar", "Krypton": "Kr", "Xenon": "Xe"}
AUTHORS = {
    "This work": "this_work",
    "Burns et al. [37]": "burns_1986",
    "Kerl [33]": "kerl_1982",
    "Koide et al. [40]": "koide_1982",
    "Buckingham et al. [38]": "buckingham_graham_1974",
    "Koch [41]": "koch_1949",
    "Kumar et al. [42]": "kumar_1985",
    "Bridge et al. [43]": "bridge_1964",
}
WAVELENGTH = "632.99"  # nm, the caption's
_CELL = re.compile(r"^([\d.]+)(?: \((\d+)\))?$")


def _sigma(value: str, digits: str) -> str:
    """Return a last-digits standard deviation as a number in the value's unit."""
    if not digits:
        return ""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and author."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE_1)[1:]:
        gas = GASES.get(cells[0], gas)
        author = AUTHORS[cells[1].rstrip("†‡§").strip()]
        match = _CELL.match(cells[2].replace("·", ".").strip())
        if match is None:
            msg = f"unexpected Table 1 cell {cells[2]!r}"
            raise ValueError(msg)
        value = number(match[1])
        point = Point(
            WAVELENGTH,
            WAVELENGTH,
            value,
            _sigma(value, match[2] or ""),
            " | ".join(cells),
        )
        files[f"{gas}/hohm_1990__{author}_{gas}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rourke (2021): recommended A_eps, A_R and the nitrogen B_R (Table 3).

Table 3 (``tbl_9524b9b87512``, cells file only; image checked) lists the recommended
reference data for helium, neon, argon and nitrogen with their best reference. Only the
rows with a number and a curated quantity are extracted: A_eps (static molar
polarizability, He, Ne, Ar), A_R (molar refractivity at 632.9908 nm, note a; He, Ne,
Ar, N2) and B_R of N2 (second refractivity virial coefficient at 302.919 K, 632.9908
nm). A_mu (molar magnetic susceptibility), C_eps, C_R, D_eps and the density virial
coefficients have no quantity in the vocabulary; the rows that point at supplementary
tables hold no number. One file per row. Values are printed with digit-group spaces
(``0.517 254 08(5)``): the spaces are removed and the parenthesised last-digit standard
uncertainty goes to the sigma column; the unit, the temperature and the note mark are
stripped (the record keeps them).
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Rourke_2021_PerspectiveRefractiveIndex"
TABLE = "tbl_9524b9b87512"
GASES = {"He": "He", "Ne": "Ne", "Ar": "Ar", "N_2": "N2"}
# Property -> (file stem, abscissa: static or the 632.9908 nm of note a).
PROPERTIES = {
    "A_{\\varepsilon}": ("A_eps", "inf"),
    "A_{R}": ("A_R", "632.9908"),
    "B_{R}": ("B_R", "632.9908"),
}
_VALUE = re.compile(
    r"^(?P<number>-?[\d. ]+)\((?P<digits>\d+)\) cm\^[36] mol\^\{-[12]a?\}"
)
_TEMPERATURE = re.compile(r"at (?P<t>[\d.]+) K")


def _split(number: str, digits: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty in the value's unit."""
    value = number.replace(" ", "")
    decimals = len(value.split(".")[1]) if "." in value else 0
    return value, f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file per extracted row of Table 3."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        quantity, gas, cell = cells[1].strip(), cells[2].strip(), cells[5].strip()
        if quantity not in PROPERTIES:
            continue
        match = _VALUE.match(cell)
        if match is None:
            continue  # a pointer to a supplementary table, no number
        stem, abscissa = PROPERTIES[quantity]
        value, sigma = _split(match["number"], match["digits"])
        found = _TEMPERATURE.search(cell)
        row = found["t"] if found else ("static" if abscissa == "inf" else abscissa)
        name = GASES[gas]
        path = f"{name}/rourke_2021__{stem}_{name}.csv"
        files[path] = [Point(row, abscissa, value, sigma, " | ".join(cells))]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

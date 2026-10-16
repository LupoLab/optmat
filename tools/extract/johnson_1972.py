# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Johnson & Christy (1972): n and k of Cu, Ag and Au, 0.64-6.60 eV.

Table I (``tbl_c125cb44d750``, p. 4374) prints, against the photon energy (eV; no
wavelength column), n and k of copper, silver and gold, and one pair of error columns,
Delta n and Delta k, common to the three metals ("the approximate errors in n and k",
caption): 49 rows. The library publishes it as JSON and HTML only (a merged "Copper |
Silver | Gold | Error" header over "n | k ..."), so it is read from the JSON cells.
Every cell was compared with the image (2026-09-28); they agree. The errors are printed
with a "+-" sign ("±0.18"): the script drops the sign and writes the number as each
point's sigma, for every metal's n (Delta n) and k (Delta k). The table has no redundant
column (no epsilon, R or wavelength), so the image is the only check.

Table II (``tbl_d5e8fa5d3264``, p. 4376) prints the optical mass m_0 (electron masses)
and the relaxation time tau (s) of each metal as "value +- error" (tau as "(6.9 \pm
0.7) \times 10^{-15}"): the value and its error are split here, and tau's column
factor is the dataset's scale. The CSV agrees with the image.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Johnson_1972_OpticalConstantsNoble"
TABLE_I = "tbl_c125cb44d750"
TABLE_II = "tbl_d5e8fa5d3264"
HEADER_I = [
    ["eV", "Copper", "", "Silver", "", "Gold", "", "Error", ""],
    ["", "n", "k", "n", "k", "n", "k", "Δn", "Δk"],
]
# The metal's n column in Table I; k follows it.
METALS = {"Cu": 1, "Ag": 3, "Au": 5}
ERRORS = {"n": 7, "k": 8}
PLUS_MINUS = "±"
HEADER_II = ["", "m_{0} (electron masses)", "\\tau (sec)"]
NAMES = {"Copper": "Cu", "Silver": "Ag", "Gold": "Au"}
# "1.49 \pm 0.06" (m_0) and "(6.9 \pm 0.7) \times 10^{-15}" (tau): value and error.
_VALUE_ERROR = {
    "m0": re.compile(r"^([\d.]+) \\pm ([\d.]+)$"),
    "tau": re.compile(r"^\(([\d.]+) \\pm ([\d.]+)\) \\times 10\^\{-15\}$"),
}


def _error(cell: str) -> str:
    """Return a printed error without its "+-" sign."""
    if not cell.startswith(PLUS_MINUS):
        msg = f"an error cell without its sign: {cell!r}"
        raise ValueError(msg)
    return number(cell.removeprefix(PLUS_MINUS))


def _table_i() -> dict[str, list[Point]]:
    """Return n and k of each metal against the photon energy, with the errors."""
    rows = json_table(DIRECTORY, TABLE_I)
    if rows[:2] != HEADER_I:
        msg = f"{TABLE_I}: unexpected header {rows[:2]}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"{metal}/johnson_1972__{part}_{metal}.csv": []
        for metal in METALS
        for part in ("n", "k")
    }
    for cells in rows[2:]:
        energy = number(cells[0])
        printed = " | ".join(cells)
        for metal, first in METALS.items():
            for part, offset in (("n", 0), ("k", 1)):
                files[f"{metal}/johnson_1972__{part}_{metal}.csv"].append(
                    Point(
                        energy,
                        energy,
                        number(cells[first + offset]),
                        _error(cells[ERRORS[part]]),
                        printed,
                    )
                )
    return files


def _table_ii() -> dict[str, list[Point]]:
    """Return each metal's optical mass and relaxation time, with their errors."""
    rows = table(DIRECTORY, TABLE_II)
    if rows[0] != HEADER_II or [row[0] for row in rows[1:]] != list(NAMES):
        msg = f"{TABLE_II}: unexpected layout {rows}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        metal = NAMES[cells[0]]
        printed = " | ".join(cells)
        for part, column in (("m0", 1), ("tau", 2)):
            match = _VALUE_ERROR[part].match(cells[column].strip())
            if match is None:
                msg = f"{TABLE_II}: unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            value, error = number(match.group(1)), number(match.group(2))
            files[f"{metal}/johnson_1972__{part}_{metal}.csv"] = [
                Point("1", "", value, error, printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the data files of Tables I and II."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

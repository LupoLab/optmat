# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Thalman, Zarzana, Tolbert & Volkamer (2014): Rayleigh cross sections.

Table 2 of the paper (``tbl_576e4a7df80e``, CSV; image checked, cells agree) gives, at
370.0, 405.8, 532.2 and 660.0 nm, the refractive-index-based ("Theory") and measured
("Exp") Rayleigh scattering cross sections of N2, Ar, O2 and air in 10^-27 cm^2 (the
caption prints "cm^-3 x 10^27"; the erratum reprints it as "cm^2 x 10^27"). The
measured cells carry their uncertainty in the last digit in parentheses, ``20.2(3)``:
one sigma for the BBCES values (Ar, O2, air), the range of repeat runs for the CRDS
values (N2; note c). A dash is no value. Air is not a registered material and is not
extracted (the record's notes give its values).

The values are written as printed in the 2014 table. The erratum (JQSRT 189, 281
(2017), appended to the library PDF; its reprinted Table 2 is ``tbl_c363d0a7b75b``)
changes six Theory cells; those corrections are ``[[errata]]`` entries of the record,
not edits here. The N2 and O2 measurements are one file per wavelength, because each
wavelength has its own King factor.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Thalman_2014_RayleighScatteringCross"
TABLE_2 = "tbl_576e4a7df80e"
# Column of each gas's Theory cell in the CSV; its Exp cell is the next column.
COLUMNS = {"N2": 1, "Ar": 3, "O2": 5}
# The measured values split by wavelength (their King factors differ), as file suffixes.
SPLIT = {"N2": "crds", "O2": "bbces"}
SUFFIX = {"370.0": "370", "405.8": "405", "532.2": "532", "660.0": "660"}
_DIGITS = re.compile(r"^(?P<value>\d+\.\d+)\((?P<digits>\d+)\)$")


def _measured(cell: str) -> tuple[str, str]:
    """Return a measured value and its last-digit uncertainty in the value's unit."""
    match = _DIGITS.match(cell.strip())
    if match is None:
        msg = f"unreadable measured cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["value"], match["digits"]
    decimals = len(value.split(".")[1])
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the theory and measured cross sections of N2, Ar and O2."""
    rows = table(DIRECTORY, TABLE_2)
    assert [rows[0][c] for c in COLUMNS.values()] == list(COLUMNS), rows[0]
    assert rows[1][1] == "Theory" and rows[1][2].startswith("Exp"), rows[1]
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for gas, column in COLUMNS.items():
            theory = f"{gas}/thalman_2014__theory_{gas}.csv"
            files.setdefault(theory, []).append(
                Point(wavelength, wavelength, number(cells[column]), "", printed)
            )
            if cells[column + 1].strip() == "-":
                continue
            value, sigma = _measured(cells[column + 1])
            if gas in SPLIT:
                name = f"{SPLIT[gas]}_{gas}_{SUFFIX[wavelength]}"
            else:
                name = f"bbces_{gas}"
            files.setdefault(f"{gas}/thalman_2014__{name}.csv", []).append(
                Point(wavelength, wavelength, value, sigma, printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

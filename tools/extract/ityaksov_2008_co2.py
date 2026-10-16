# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ityaksov, Linnartz and Ubachs (2008): CO2 deep-UV extinction.

The only table (``tbl_d53eeed0ae27``; image checked) lists, at 22 vacuum wavelengths
from 197.70 to 270.15 nm, the calculated Rayleigh cross section and the measured
extinction cross section (+- one sigma), both in 10^-25 cm^2. The CSV's first row is
the header, its second the units (``(nm)``, ``(cm2)``), dropped here.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Ityaksov_2008_DeepUVAbsorption"
TABLE = "tbl_d53eeed0ae27"
_VALUE_SIGMA = re.compile(r"^\s*([\d.]+)\s*±\s*([\d.]+)\s*$")
PREFIX = "CO2/ityaksov_2008_co2__"


def extract() -> dict[str, list[Point]]:
    """Return the calculated and measured cross sections as data files."""
    calculated, measured = [], []
    for cells in table(DIRECTORY, TABLE)[2:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        calculated.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        match = _VALUE_SIGMA.match(cells[2])
        if match is None:
            msg = f"unexpected cell {cells[2]!r}"
            raise ValueError(msg)
        measured.append(
            Point(wavelength, wavelength, match.group(1), match.group(2), printed)
        )
    return {
        f"{PREFIX}extinction_CO2.csv": measured,
        f"{PREFIX}calculated_CO2.csv": calculated,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

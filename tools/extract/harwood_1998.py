# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Harwood, Burkholder & Ravishankara (1998): N2O5 cross sections, NO3 yields.

* **Table 2** (``tbl_52089e995460``, CSV; image checked, p. 1312): the N2O5
  cross section, 208-398 nm every 2 nm, in 10^-20 cm^2 molecule^-1, printed in three
  column pairs (wavelength, sigma), read down each pair. The spectrum falls
  monotonically over the whole range, which is checked.
* **Table 3** (``tbl_952e9553cc98``, cells only; image checked, p. 1313): the NO3
  yield of N2O5 photolysis by calibration method and photolysis wavelength, "value +-
  2 sigma" of the precision. Each filled cell is a row, its id the method and
  wavelength; the last row ("value used in BrONO2 analysis") is not a measurement and
  is not written.

Register: TC-974.
"""

import itertools
import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Harwood_1998_PhotodissociationBrONO2N2O5"
TABLE_2 = "tbl_52089e995460"
TABLE_3 = "tbl_952e9553cc98"
METHODS = {
    "HNO_3 photolysis": "HNO3",
    "HNO_3/F_2": "HNO3/F2",
    "O_3 photolysis": "O3",
    "O_3/Cl_2": "O3/Cl2",
}
_YIELD = re.compile(r"^([\d.]+) ± ([\d.]+)$")


def _cross_sections() -> list[Point]:
    """Return Table 2's cross sections, ordered by wavelength."""
    points = []
    for pair in range(3):
        for cells in table(DIRECTORY, TABLE_2)[2:]:
            wavelength, sigma = cells[2 * pair], cells[2 * pair + 1]
            points.append(
                Point(
                    wavelength, number(wavelength), number(sigma), "", " | ".join(cells)
                )
            )
    values = [float(p.value) for p in points]
    if [float(p.x) for p in points] != sorted(float(p.x) for p in points) or any(
        later >= earlier for earlier, later in itertools.pairwise(values)
    ):
        msg = "Table 2 is not read in wavelength order, or does not fall monotonically"
        raise ValueError(msg)
    return points


def _yields() -> list[Point]:
    """Return Table 3's NO3 quantum yields by calibration method and wavelength."""
    grid = json_table(DIRECTORY, TABLE_3)
    wavelengths = [cell.removesuffix(" nm") for cell in grid[1][1:]]
    points = []
    for cells in grid[2:]:
        if cells[0] not in METHODS:
            continue  # the values used in the BrONO2 analysis
        for wavelength, cell in zip(wavelengths, cells[1:], strict=True):
            if not cell:
                continue
            match = _YIELD.match(cell)
            if match is None:
                msg = f"unreadable cell {cell!r}"
                raise ValueError(msg)
            row = f"{METHODS[cells[0]]} {wavelength}"
            points.append(
                Point(
                    row, wavelength, match.group(1), match.group(2), " | ".join(cells)
                )
            )
    return points


def extract() -> dict[str, list[Point]]:
    """Return the N2O5 cross sections and the NO3 yields by method."""
    return {
        "N2O5/harwood_1998__cross_section.csv": _cross_sections(),
        "N2O5/harwood_1998__no3_yield_by_method.csv": _yields(),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Simmons (1978): N2 and CO2 at 22002, 15803 and 944 cm^-1 (the table).

The only table (uncaptioned in the library, no CSV rendering) is read from its
JSON cells and was checked against the image. Its heading row gives the three wave
numbers, each with two columns: the Lorentz-Lorenz function at 300 K (cm^3/mol) and
the refractive index at NTP, each ``value±uncertainty``. The refractive index is
printed as n (``1.000302±0.000003``); the script stores n - 1 by exact decimal
subtraction of 1 (``0.000302``), the uncertainty as printed. Propane is not curated.
"""

from decimal import Decimal

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Simmons_1978_RefractiveIndexLorentz"
TABLE = "tbl_fe193a93ae5a"
GASES = {"Nitrogen": "N2", "Carbon-dioxide": "CO2"}


def extract() -> dict[str, list[Point]]:
    """Return the LL functions and the refractivities per gas."""
    rows = json_table(DIRECTORY, TABLE)
    wavenumbers = [cell.removesuffix("cm-1").strip() for cell in rows[1][1::2]]
    files: dict[str, list[Point]] = {}
    for cells in rows[3:]:
        if cells[0] not in GASES:
            continue
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for index, wavenumber in enumerate(wavenumbers):
            ll, ll_sigma = cells[1 + 2 * index].split("±")
            n, n_sigma = cells[2 + 2 * index].split("±")
            files.setdefault(f"{gas}/simmons_1978__LL_{gas}.csv", []).append(
                Point(wavenumber, wavenumber, number(ll), ll_sigma.strip(), printed)
            )
            refractivity = str(Decimal(number(n)) - 1)
            files.setdefault(f"{gas}/simmons_1978__ntp_{gas}.csv", []).append(
                Point(wavenumber, wavenumber, refractivity, n_sigma.strip(), printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Essen (1953): microwave refractivities of O2, N2, H2, D2 and He.

The only table (``tbl_648ebcc40786``, cells file only; image checked) gives
``(n - 1) 10^6`` of dry gases at 0 C and 760 mm Hg (water vapour at 20 C, 10 mm Hg):
Essen and Froome's 1951 values at 24 000 Mc/s with their +-, the present values (the
column is headed '9000 Mc/s'; the text gives the apparatus frequency as 9200 Mc/s),
the number of measurements, the standard deviation and the limits of error including
systematic errors. The decimal point is printed raised (``266·2``). The present
values are extracted with the limits of error as sigma, Essen and Froome's with
theirs. Water vapour and air are not curated species.

The wavelength column holds the vacuum wavelength c/f in m (the Newell and Baird
record's convention): c/(9200 Mc/s) for the present values and c/(24 000 Mc/s) for
Essen and Froome's; the printed frequency stays in the printed cells.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Essen_1953_RefractiveIndicesWater"
TABLE = "tbl_648ebcc40786"
C = 299_792_458.0  # speed of light, m/s
GASES = {
    "Oxygen": "O2",
    "Nitrogen": "N2",
    "Hydrogen": "H2",
    "Deuterium": "D2",
    "Helium": "He",
}
PRESENT = 9200e6  # Hz, the apparatus frequency (text, p. 3)
ESSEN_FROOME = 24000e6  # Hz


def _decimal(cell: str) -> str:
    """Read the printed raised decimal point as a decimal point."""
    return number(cell.replace("\N{MIDDLE DOT}", ".").replace(" ", ""))


def extract() -> dict[str, list[Point]]:
    """Return the present and Essen and Froome's values as data files."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        limit = _decimal(cells[5].removeprefix("\N{PLUS-MINUS SIGN}"))
        wavelength = f"{C / PRESENT:.6g}"
        files[f"{gas}/essen_1953__microwave_{gas}.csv"] = [
            Point("9200 Mc/s", wavelength, _decimal(cells[2]), limit, printed)
        ]
        if cells[1]:
            value, _, sigma = cells[1].partition("\N{PLUS-MINUS SIGN}")
            files[f"{gas}/essen_1953__essen_froome_1951_{gas}.csv"] = [
                Point(
                    "24000 Mc/s",
                    f"{C / ESSEN_FROOME:.6g}",
                    _decimal(value),
                    _decimal(sigma),
                    printed,
                )
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

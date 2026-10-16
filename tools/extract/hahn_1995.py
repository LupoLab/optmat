# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hahn & Lee (1995): nonresonant CARS susceptibilities (Tables 1-3).

- Table 1 (``tbl_50a8b1a2127a``; the cells agree with the image, whose header the
  library garbles: "10218" is 10^-18, "cm21" is cm^-1) gives per gas the measured
  ratio of the interference-fringe slope to argon's and the effective chi_NR(3) in
  10^-18 cm^3/(erg amagat), each with its uncertainty in parentheses, and the Raman
  shift of the measurement.
- Table 2 (``tbl_4564b15b3033``, JSON cells only) gives the calculated off-resonant
  Raman terms at 2157 and 2334 cm^-1. The library lost every minus sign of this column,
  reading the minus glyph as a "2": the image prints -12.00, -0.35, -0.62, -0.67,
  -2.03 and -0.24 where the cells hold 212.00, 20.35, 20.62, 20.67, 22.03 and 20.24
  (and 4.17, 2.30, 12.5, 0.77 without sign). The rows in ``NEGATIVE`` are restored
  here, after checking that each cell is "2" followed by the printed magnitude.
- Table 3 (``tbl_c1933c977eba``; the cells agree with the image) gives the electronic
  chi_NR(3) of this experiment and of six other determinations; argon's "(9.46)" in
  the first column is the reference value, not a result, and is not extracted.

Acetylene, propane and Freon (CF4) are not curated species and are skipped.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Hahn_1995_MeasurementNonresonantThird"
TABLE_1 = "tbl_50a8b1a2127a"
TABLE_2 = "tbl_4564b15b3033"
TABLE_3 = "tbl_c1933c977eba"
PUMP = "532"  # nm, the injection-seeded frequency-doubled Nd:YAG (p. 1023)
GASES = {
    "Argon": "Ar",
    "Carbon dioxide": "CO2",
    "Methane": "CH4",
    "Nitrogen": "N2",
    "Oxygen": "O2",
    "Carbon monoxide": "CO",
    "Hydrogen": "H2",
}
# Table 2: (gas, Raman band) rows whose printed off-resonant term is negative (image).
NEGATIVE = {("CO2", "1285"), ("CO2", "1388"), ("O2", "1555"), ("CO", "2143")}
# Table 3 column -> (dataset prefix, pump wavelength in nm as printed in the header)
TABLE_3_COLUMNS = {
    1: ("electronic", "532"),
    2: ("irs_rosasco_hurst", "676"),
    3: ("cars_farrow", "532"),
    4: ("eshg_shelton_rice", "694"),
    5: ("cars_used", "532"),
    6: ("cars_lundeen", "532"),
    7: ("cars_rado", "694"),
}
_PAREN = re.compile(r"^(?P<value>[\d.]+) \((?P<sigma>[\d.]+)\)$")


def _paren(cell: str) -> tuple[str, str]:
    """Return the value and the parenthesized uncertainty of a Table 1 cell."""
    match = _PAREN.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(match["value"]), number(match["sigma"])


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and quantity."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_1)[2:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        shift = cells[3].strip()
        ratio, ratio_sigma = _paren(cells[1])
        value, sigma = _paren(cells[2])
        files[f"{gas}/hahn_1995__ratio_{gas}.csv"] = [
            Point(shift, PUMP, ratio, ratio_sigma, printed)
        ]
        files[f"{gas}/hahn_1995__effective_{gas}.csv"] = [
            Point(shift, PUMP, value, sigma, printed)
        ]
    for cells in json_table(DIRECTORY, TABLE_2)[3:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        band = cells[1]
        cell = cells[4] if cells[4] != "-" else cells[5]
        if (gas, band) in NEGATIVE:
            if not cell.startswith("2"):
                msg = f"expected a lost minus in {cell!r}"
                raise ValueError(msg)
            cell = "-" + cell[1:]
        files.setdefault(f"{gas}/hahn_1995__off_resonant_{gas}.csv", []).append(
            Point(band, PUMP, number(cell), "", printed)
        )
    for cells in table(DIRECTORY, TABLE_3)[1:]:
        gas = GASES.get(cells[0])
        if gas is None:
            continue
        printed = " | ".join(cells)
        for column, (prefix, wavelength) in TABLE_3_COLUMNS.items():
            cell = cells[column].strip()
            if cell in {"", "-"} or cell.startswith("("):  # none, or the reference
                continue
            files[f"{gas}/hahn_1995__{prefix}_{gas}.csv"] = [
                Point("chi_NR", wavelength, number(cell), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

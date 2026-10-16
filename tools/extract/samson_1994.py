# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Samson et al. (1994): helium photoionization cross sections, 24.6 eV-8 keV.

- Table 1 (``tbl_95893c11a7e1``): energy (eV) and cross section (Mb) in three column
  pairs, 24.587-280 eV, read column pair by column pair (the print order). Split at
  120 eV: the values to 120 eV are the best fit to the authors' measurements
  (``measured``); 125-280 eV are that fit extrapolated to Watson's and Denne's data
  (``extrapolated``).
- Table 2 (``tbl_168a48a2c629``): the predicted cross sections 280-8000 eV in barns
  (two column pairs), stored with scale 1e-6 (b to Mb); its 280 eV row repeats Table
  1's last value.

The images agree with the CSV cells. Table 3 (oscillator strengths of energy ranges)
is not extracted.
"""

from _library import Point, number, table, write

DIRECTORY = "Samson_1994_PrecisionMeasurementsAbsolute"
TABLE_1 = "tbl_95893c11a7e1"
TABLE_2 = "tbl_168a48a2c629"
PREFIX = "He/samson_1994__"


def _pairs(table_id: str, n_pairs: int) -> list[tuple[str, str]]:
    rows = table(DIRECTORY, table_id)[1:]
    return [
        (cells[2 * k].strip(), cells[2 * k + 1].strip())
        for k in range(n_pairs)
        for cells in rows
        if cells[2 * k].strip()
    ]


def extract() -> dict[str, list[Point]]:
    """Return Table 1 (measured and extrapolated) and Table 2."""
    measured, extrapolated = [], []
    for energy, value in _pairs(TABLE_1, 3):
        point = Point(energy, number(energy), number(value), "", f"{energy} | {value}")
        (measured if float(energy) <= 120 else extrapolated).append(point)
    predicted = [
        Point(energy, number(energy), number(value), "", f"{energy} | {value}")
        for energy, value in _pairs(TABLE_2, 2)
    ]
    return {
        f"{PREFIX}measured.csv": measured,
        f"{PREFIX}extrapolated.csv": extrapolated,
        f"{PREFIX}predicted.csv": predicted,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

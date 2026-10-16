# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Carusotto et al. (1985): Kerr constants of six gases at 514.5 nm.

Both tables checked against their images (no differences in the cells).

* Table I (``tbl_77609c382209``, CSV): the individual measurements, per gas the density
  (amagat), the field E_0 (statV/cm), the ellipticity and K x 10^17 +- its statistical
  error. One file per gas (``table1_<gas>``), rows "1", "2", ... in printed order, with
  the density and the field as extra columns; the ellipticity stays in ``printed``.
* Table II (``tbl_569a1459caec``, JSON cells only): the averaged K x 10^17 of this work
  (the last six rows, ``table2_<gas>``) and, above them, earlier values converted by the
  authors (``table2_ref<n>_<gas>``, with the temperature appended where one gas has two
  rows). The header prints the unit "StatV^-1 cm^2 amagat^-1"; see the record. Argon is
  printed "A".
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Carusotto_1985_MeasurementKerrConstant"
TABLE_1 = "tbl_77609c382209"
TABLE_2 = "tbl_569a1459caec"
GASES = {"He": "He", "A": "Ar", "O_2": "O2", "N_2": "N2", "H_2": "H2", "D_2": "D2"}
WAVELENGTH = "514.5"  # nm, the argon-ion line of this work
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
REFERENCE_WAVELENGTHS = {"1": "632.8", "2": "632.8", "6": "632.8", "13": "650"}
TWO_ROWS = {("N_2", "6"), ("H_2", "2")}  # one gas, two temperatures: two datasets


def _split(cell: str) -> tuple[str, str]:
    """Split "0.129 +- 0.011" into value and error, as printed."""
    value, _, error = cell.partition(PLUS_MINUS)
    return number(value), number(error) if error.strip() else ""


def _table1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_1)[1:]:
        gas = GASES[cells[0]]
        name = f"{gas}/carusotto_1985__table1_{gas}.csv"
        value, error = _split(cells[4])
        row = str(len(files.get(name, [])) + 1)
        extra = {"density": number(cells[1]), "field": number(cells[2])}
        files.setdefault(name, []).append(
            Point(row, WAVELENGTH, value, error, " | ".join(cells), extra)
        )
    return files


def _table2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, TABLE_2)
    if not grid[0][1].startswith("Kerr constant"):
        msg = f"unexpected Table II header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        gas = GASES[cells[0]]
        value = number(cells[1])
        error = number(cells[2].replace(PLUS_MINUS, "")) if cells[2] else ""
        printed = " | ".join(cells)
        reference = cells[5].strip("() ")
        if reference == "present work":
            files[f"{gas}/carusotto_1985__table2_{gas}.csv"] = [
                Point("present work", WAVELENGTH, value, error, printed)
            ]
            continue
        stem = f"table2_ref{reference}_{gas}"
        if (cells[0], reference) in TWO_ROWS:
            stem += "_" + cells[4].split(".")[0]
        files[f"{gas}/carusotto_1985__{stem}.csv"] = [
            Point(
                f"ref {reference}",
                REFERENCE_WAVELENGTHS[reference],
                value,
                error,
                printed,
            )
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I and Table II files."""
    return _table1() | _table2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

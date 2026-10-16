# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Li et al. (1989): saturation intensities and THG chi(3) of the rare gases.

Phys. Rev. A 39, 5751-5761 (1 June 1989); journal page = PDF page + 5750. A scan with
a text layer. Two tables, each a row per quantity and a column per gas (Xe, Kr, Ar,
Ne, He), every cell printed as ``x[y]`` = x * 10^y (the captions' notation):

* **Table I** (``tbl_8884f312bffe``, p. 5752): the saturation intensities I_s (W
  cm^-2) at 1064 nm and 36 ps, from the authors' earlier ion measurements (Ref. 12,
  L'Huillier et al., J. Phys. B 16, 1363 (1983); 'The saturation intensities for He
  and Ne had been overestimated by a factor of 2'), +-50 %. One file per gas.
* **Table II** (``tbl_c24350e35014``, p. 5754): Delta k_3 and Delta k_5 (cm^-1, from
  Leonard's refractive indices, Ref. 22), chi_0^(3) (Lehmeier et al.'s THG values,
  Ref. 29, times 4), chi^(3) and chi^(5) (esu, 'extracted from the experimental
  results and normalized to chi_0^(3) in Xe', +-50 %). The library's CSV labels the
  fourth and fifth rows 'chi^(3)_0' and 'chi^(5)_0'; the table image prints chi^(3)
  and chi^(5) without the subscript: the rows are taken by position (the label of
  each is checked against the CSV as it is). Only the chi^(3) row is written, for Kr,
  Ar, Ne and He (Xe's equals chi_0^(3) by the normalisation). Xe's Delta k_3 prints
  '2.[-1]' on the image and '2[-1]' in the CSV; it is not written.

Each value is the cell's mantissa times its power of ten, as text (``1.2e13``); the
cell is kept in ``printed``. Every cell was checked against the table images
(2026-10-09).

:func:`checks` tests the chi_0^(3) row against lehmeier_1985's own esu values
(``gamma_thg_esu_<gas>``, 10^-39 esu) times 4, as the caption says: each printed
mantissa is within one unit of its last figure of the product (Ne's 6.56 and He's
3.656 are printed 6.5 and 3.6, truncated rather than rounded), and that Xe's chi^(3)
equals its chi_0^(3).
"""

import csv
import re
from typing import Final

from _library import DATA, Point, table, write

DIRECTORY: Final = "Li_1989_MultipleHarmonicGeneration"
TABLE_I: Final = "tbl_8884f312bffe"
TABLE_II: Final = "tbl_c24350e35014"
GASES: Final = ("Xe", "Kr", "Ar", "Ne", "He")
"""The columns of both tables, in order."""
TABLE_II_ROWS: Final = (
    r"\Delta k_{3} (cm^{-1})",
    r"\Delta k_{5} (cm^{-1})",
    r"\chi^{(3)}_{0} (esu)",
    r"\chi^{(3)}_{0} (esu)",
    r"\chi^{(5)}_{0} (esu)",
)
"""Table II's row labels as the library's CSV holds them (rows 4 and 5 mislabelled)."""
CHI3_ROW: Final = 3
"""The index of the chi^(3) row (the image's fourth row) among Table II's data rows."""
CHI0_ROW: Final = 2
"""The index of the chi_0^(3) row (Lehmeier's values times 4)."""
CELL: Final = re.compile(r"^(\d+(?:\.\d*)?)\[(-?\d+)\]$")
"""A cell, ``x[y]`` = x * 10^y."""
LEHMEIER_SCALE: Final = 1e-39
"""The scale of lehmeier_1985's ``gamma_thg_esu_<gas>`` files (its datasets' scale)."""


def _cell(cell: str) -> tuple[float, float, str]:
    """Return a cell's mantissa, its power of ten and the value as text."""
    match = CELL.match(cell.strip())
    if match is None:
        msg = f"an unexpected cell {cell!r}"
        raise ValueError(msg)
    mantissa, power = match.groups()
    return float(mantissa), 10.0 ** int(power), f"{mantissa}e{power}"


def _rows(table_id: str, labels: tuple[str, ...]) -> list[list[str]]:
    """Return a table's data rows, checking the header and the row labels."""
    grid = table(DIRECTORY, table_id)
    if tuple(c.strip() for c in grid[0][1:]) != GASES:
        msg = f"{table_id}: the gas columns moved: {grid[0]}"
        raise ValueError(msg)
    rows = grid[1:]
    if tuple(row[0].strip() for row in rows) != labels:
        msg = f"{table_id}: the rows moved: {[row[0] for row in rows]}"
        raise ValueError(msg)
    return rows


def saturation_rows() -> list[list[str]]:
    """Return Table I's one data row."""
    return _rows(TABLE_I, (r"I_{s} (W cm^{-2})",))


def table_ii_rows() -> list[list[str]]:
    """Return Table II's five data rows."""
    return _rows(TABLE_II, TABLE_II_ROWS)


def _lehmeier_esu(gas: str) -> tuple[float, str]:
    """Return lehmeier_1985's gamma(3) of a gas in esu, and its file's value string."""
    name = f"lehmeier_1985__gamma_thg_esu_{gas}.csv"
    path = DATA / "values" / "third_order" / gas / name
    with path.open(newline="", encoding="utf-8") as handle:
        (row,) = list(csv.DictReader(handle))
    return float(row["value"]) * LEHMEIER_SCALE, row["value"]


def checks() -> dict[str, object]:
    """Return the chi_0^(3) row against lehmeier_1985 times 4.

    Raises unless every chi_0^(3) mantissa is within one unit of its last printed
    figure of 4 x lehmeier_1985's esu value, and Xe's chi^(3) equals its chi_0^(3).
    """
    rows = table_ii_rows()
    out: dict[str, object] = {}
    failures = []
    for j, gas in enumerate(GASES, start=1):
        mantissa, power, _ = _cell(rows[CHI0_ROW][j])
        lehmeier, printed = _lehmeier_esu(gas)
        product = 4 * lehmeier / power
        # One unit of the last printed figure: the paper truncates as well as rounds.
        unit = 10.0 ** -len(rows[CHI0_ROW][j].split("[")[0].partition(".")[2])
        if abs(product - mantissa) > unit + 1e-9:
            failures.append(f"{gas}: chi_0 {mantissa} against 4 x {printed} e-39")
        out[f"{gas}: chi_0^(3) mantissa, 4 x lehmeier_1985"] = (
            mantissa,
            round(product, 3),
        )
    if rows[CHI3_ROW][1] != rows[CHI0_ROW][1]:
        failures.append("Xe: chi^(3) is not chi_0^(3)")
    out["Xe: chi^(3) = chi_0^(3) (the normalisation)"] = rows[CHI3_ROW][1]
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return one file per gas of I_s, and of chi^(3) for Kr, Ar, Ne and He."""
    checks()
    files: dict[str, list[Point]] = {}
    (saturation,) = saturation_rows()
    for j, gas in enumerate(GASES, start=1):
        files[f"{gas}/li_1989__isat_{gas}.csv"] = [
            Point(
                gas,
                "1064",
                _cell(saturation[j])[2],
                "",
                " | ".join(c.strip() for c in saturation),
            )
        ]
    chi3 = table_ii_rows()[CHI3_ROW]
    for j, gas in enumerate(GASES, start=1):
        if gas == "Xe":
            continue
        files[f"{gas}/li_1989__chi3_thg_{gas}.csv"] = [
            Point(
                gas,
                "1.064",
                _cell(chi3[j])[2],
                "",
                " | ".join(c.strip() for c in chi3),
            )
        ]
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

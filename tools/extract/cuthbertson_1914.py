# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract C. and M. Cuthbertson (1914): ozone and nitrous oxide (Tables IX-XI, XIV).

The library stores the paper under a title with a junk prefix
(``Cuthbertson_1914_InfoGettingLost``). All values are (mu - 1) at the 'theoretical'
density D (see the source record). Decimal points are printed as raised dots ('·'),
converted here.

- Table IX (``tbl_6f597b2ac7a8``, JSON only) and Table X (``tbl_f49a5176a7d2``, JSON
  only): ozone refractivity x 10^6 at the green mercury line, one per experiment, by
  two methods of estimating the ozone content. The image shows the remark 'Not
  reliable, percentage of ozone too small' against experiment 7 of Table IX and 'Not
  very trustworthy' against experiment 3 of Table X, the two values printed in square
  brackets; the library attaches both remarks to experiment 1. The bracketed values
  are written to separate files.
- Table XI (``tbl_4bf7db766c5d``, JSON only, image checked): refractivities x 10^8 of
  ozone-oxygen mixtures at eight wavelengths; only the last two rows are extracted:
  oxygen (reproduced from the authors' 1909 paper) and ozone (derived from the mixture
  means and the oxygen row).
- The unlabelled table ``tbl_6f5db40b38bd`` (JSON only): seven determinations of N2O's
  refractivity x 10^7 at the green line in three series. The library merges the values
  of each series into one cell ('50925102', '509750985099', '50875091'); the image
  shows one four-digit value per line (5092, 5102 / 5097, 5098, 5099 / 5087, 5091), so
  the cells are split into groups of four digits.
- Table XIV (``tbl_d05d1c4d670a``, JSON only; its caption lacks the gas, which the text
  shows to be N2O): observed and calculated refractivity x 10^8 at eight wavelengths.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Cuthbertson_1914_InfoGettingLost"
TABLE_IX = "tbl_6f597b2ac7a8"
TABLE_X = "tbl_f49a5176a7d2"
TABLE_XI = "tbl_4bf7db766c5d"
N2O_SERIES = "tbl_6f5db40b38bd"
TABLE_XIV = "tbl_d05d1c4d670a"

GREEN = "5461"  # the green mercury line, lambda = 5461 (air) in the text


def _number(cell: str) -> str:
    return number(cell.replace("\N{MIDDLE DOT}", ".").strip("[]"))


def _ozone_experiments(table_id: str, name: str) -> dict[str, list[Point]]:
    kept: list[Point] = []
    bracketed: list[Point] = []
    for cells in json_table(DIRECTORY, table_id)[1:]:
        if not cells[0]:  # the 'Mean' row
            continue
        printed = " | ".join(cells[:3])
        point = Point(f"experiment {cells[0]}", GREEN, _number(cells[2]), "", printed)
        (bracketed if cells[2].startswith("[") else kept).append(point)
    return {
        f"O3/cuthbertson_1914__{name}.csv": kept,
        f"O3/cuthbertson_1914__{name}_bracketed.csv": bracketed,
    }


def _table_xi() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_XI)
    wavelengths = rows[1][1:9]
    files: dict[str, list[Point]] = {}
    for cells in rows:
        if cells[0].startswith("Refractivities of O_2"):
            gas, name = "O2", "oxygen_1909"
        elif cells[0].startswith("Refractivities of O_3"):
            gas, name = "O3", "ozone_dispersion"
        else:
            continue
        printed = " | ".join(cells[:9])
        files[f"{gas}/cuthbertson_1914__{name}.csv"] = [
            Point(w, w, _number(v), "", printed)
            for w, v in zip(wavelengths, cells[1:9], strict=True)
        ]
    return files


def _n2o_series() -> dict[str, list[Point]]:
    points: list[Point] = []
    for cells in json_table(DIRECTORY, N2O_SERIES)[1:4]:
        merged = cells[1]
        values = [merged[i : i + 4] for i in range(0, len(merged), 4)]
        printed = " | ".join(cells)
        points.extend(
            Point(f"series {cells[0]} value {k}", GREEN, number(v), "", printed)
            for k, v in enumerate(values, start=1)
        )
    return {"N2O/cuthbertson_1914__green_line.csv": points}


def _table_xiv() -> dict[str, list[Point]]:
    observed: list[Point] = []
    calculated: list[Point] = []
    for cells in json_table(DIRECTORY, TABLE_XIV)[2:]:
        wavelength = _number(cells[0])
        printed = " | ".join(cells)
        observed.append(Point(wavelength, wavelength, _number(cells[1]), "", printed))
        calculated.append(Point(wavelength, wavelength, _number(cells[2]), "", printed))
    return {
        "N2O/cuthbertson_1914__dispersion.csv": observed,
        "N2O/cuthbertson_1914__dispersion_calculated.csv": calculated,
    }


def extract() -> dict[str, list[Point]]:
    """Return the ozone, oxygen and nitrous-oxide data files."""
    return (
        _ozone_experiments(TABLE_IX, "ozone_method1")
        | _ozone_experiments(TABLE_X, "ozone_method2")
        | _table_xi()
        | _n2o_series()
        | _table_xiv()
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

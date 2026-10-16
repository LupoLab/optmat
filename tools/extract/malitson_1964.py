# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Malitson (1964): the index of BaF2 at 25 degrees C, 0.2652-10.346 um.

- **Table II** (``tbl_85b36c3e0390``, cells only): 46 wavelengths (um, the lines' air
  wavelengths), the spectral source, the measured and the computed index at 25 degrees
  C, and the residual x 1e5. The measured and computed columns are extracted; the
  residual is their difference to the last digit (all 46 rows, checked) and is not.
- **Table I** (``tbl_a7ee026029c3``, cells only): the 1944 MIT experimental sample at
  15, 35 and 55 degrees C (nine visible lines), and -dn/dt "x 10^-6 per degree C"
  averaged over 15-55 degrees C. Two wavelengths carry the footnote mark "a"
  ("Intensity weighted mean of doublet", pdf p. 1): stripped, and named in the notes.
  The last row (Abbe numbers) is left to the notes.
- **Table IV** (``tbl_3302616c5282``): n and -dn/dlambda from the formula at regular
  wavelengths, in two column pairs (the right pair ends one row early); the n columns
  are extracted (85 points).

All cells agree with the images. Table III (``tbl_c6dc19d9a616``; the library missed its
caption, printed above it on pdf p. 2) is the record's formula; this script checks the
record's strings against its cells. Rows are the wavelengths as printed (um).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Malitson_1964_RefractivePropertiesBarium"
TABLE_1 = "tbl_a7ee026029c3"
TABLE_2 = "tbl_85b36c3e0390"
TABLE_3 = "tbl_c6dc19d9a616"
TABLE_4 = "tbl_3302616c5282"
PREFIX = "BaF2/malitson_1964__"

TABLE_3_PRINTED = [
    ["A_{1}=0.643356", "\\lambda_{1}=0.057789", "\\lambda_{1}^{2}=0.0033396"],
    ["A_{2}=0.506762", "\\lambda_{2}=0.10968", "\\lambda_{2}^{2}=0.012030"],
    ["A_{3}=3.8261", "\\lambda_{3}=46.3864", "\\lambda_{3}^{2}=2151.70"],
]
TABLE_1_COLUMNS = {1: "mit_n_15C", 2: "mit_n_35C", 3: "mit_n_55C", 4: "mit_dndt_15_55C"}
DOUBLET = "^a"


def _check_table_3() -> None:
    """Check the record's transcription of Table III against the library cells."""
    found = table(DIRECTORY, TABLE_3)
    if found != TABLE_3_PRINTED:
        msg = f"Table III reads {found}, not {TABLE_3_PRINTED}"
        raise ValueError(msg)


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    if rows[1][2:4] != ["Measured", "Computed"] or len(rows) != 48:
        msg = f"unexpected Table II layout {rows[1]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, name in ((2, "measured_25C"), (3, "computed_25C")):
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _table_1() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[1][1:4] != ["15°C", "35°C", "55°C"] or rows[-1][0] != "Abbe No.":
        msg = f"unexpected Table I layout {rows[1]!r}"
        raise ValueError(msg)
    body = rows[2:-1]
    if len(body) != 9:
        msg = f"Table I has {len(body)} wavelengths, not 9"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        # The footnote mark "a" (a doublet's intensity-weighted mean) is not a digit.
        wavelength = number(cells[0].removesuffix(DOUBLET))
        printed = " | ".join(cells)
        for column, name in TABLE_1_COLUMNS.items():
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _table_4() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_4)
    if len(rows) != 44 or rows[0][1] != "n" or rows[0][4] != "n":
        msg = f"unexpected Table IV layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    points: list[Point] = []
    # The left column pair (0.26-2.00 um), then the right (2.20-10.40 um), whose last
    # row is empty.
    for first in (0, 3):
        for cells in rows[1:]:
            if not cells[first]:
                continue
            wavelength = number(cells[first])
            printed = " | ".join(cells[first : first + 3])
            points.append(
                Point(wavelength, wavelength, number(cells[first + 1]), "", printed)
            )
    return {f"{PREFIX}table4_n_25C.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return Tables I, II and the n of Table IV."""
    _check_table_3()
    return {**_table_2(), **_table_1(), **_table_4()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

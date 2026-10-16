# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Malitson (1963): the index of CaF2 at 24 degrees C, 0.23-9.7 um.

- **Table II** (``tbl_519685ce8f61``): 46 wavelengths (um, the lines' air wavelengths),
  the spectral source, the index computed from the Table III formula at 24 degrees C,
  the residuals observed minus computed (x 1e-5) of the synthetic and of the natural
  prism, and -dn/dt x 1e6 for the mean temperature 19 degrees C. The cells agree with
  the image, signs included (the "+" of a positive residual is dropped as printed
  numbers are).
- **Table I** (``tbl_ae8cd73c9fb7``, cells only): the averaged index of two synthetic
  prisms measured in 1944 at 15, 35 and 55 degrees C (9 visible lines), and -dn/dt x
  1e6 for 15-35 and 35-55 degrees C. Two wavelengths carry the footnote mark "a"
  ("Intensity weighted mean of doublet"): a marker, not a digit; the script strips it
  and the record's notes say which rows. The last row (Abbe numbers) is not per
  wavelength and is left to the notes. The cells agree with the image.
- **Table IV** (``tbl_9273a9ecad2b``): n and -dn/dlambda computed from the formula at
  regular wavelengths, 0.22-9.80 um, printed in two column pairs; the n columns are
  extracted (84 points), -dn/dlambda has no quantity. The cells agree with the image.

Table III (the dispersion constants, cells only) is the record's formula; this script
checks the record's strings against the cells. Rows are the wavelengths as printed (um).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Malitson_1963_RedeterminationSomeOptical"
TABLE_1 = "tbl_ae8cd73c9fb7"
TABLE_2 = "tbl_519685ce8f61"
TABLE_3 = "tbl_b43d14ecd64e"
TABLE_4 = "tbl_9273a9ecad2b"
PREFIX = "CaF2/malitson_1963__"

# Table III as printed (image and cells): lambda_j, lambda_j^2 and A_j.
TABLE_3_PRINTED = (
    (
        "\\lambda_{1} = 0.050263605",
        "\\lambda_{1}^{2} = 0.002526430",
        "A_{1} = 0.5675888",
    ),
    ("\\lambda_{2} = 0.1003909", "\\lambda_{2}^{2} = 0.01007833", "A_{2} = 0.4710914"),
    ("\\lambda_{3} = 34.649040", "\\lambda_{3}^{2} = 1200.5560", "A_{3} = 3.8484723"),
)

TABLE_2_COLUMNS = {
    2: "computed_24C",
    3: "residual_synthetic_24C",
    4: "residual_natural_24C",
    5: "dndt_19C",
}
TABLE_1_COLUMNS = {
    1: "stockbarger_n_15C",
    2: "stockbarger_n_35C",
    3: "stockbarger_n_55C",
    4: "dndt_15_35C",
    5: "dndt_35_55C",
}
DOUBLET = "^a"


def _check_table_3() -> None:
    """Check the record's transcription of Table III against the library cells."""
    rows = json_table(DIRECTORY, TABLE_3)[1:]
    found = tuple(tuple(r) for r in rows)
    if found != TABLE_3_PRINTED:
        msg = f"Table III reads {found}, not {TABLE_3_PRINTED}"
        raise ValueError(msg)


def _table_2() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_2)
    if rows[0][2] != "Computed index" or len(rows) != 47:
        msg = f"unexpected Table II layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, name in TABLE_2_COLUMNS.items():
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _table_1() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[2] != ["", "15°C", "35°C", "55°C", "15-35°C", "35-55°C"]:
        msg = f"unexpected Table I header {rows[2]!r}"
        raise ValueError(msg)
    body = rows[3:]
    if len(body) != 10 or body[-1][0] != "Abbe No.":
        msg = "Table I does not have nine wavelengths before the Abbe numbers"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body[:-1]:
        # The footnote mark "a" (a doublet's intensity-weighted mean) is not a digit.
        wavelength = number(cells[0].removesuffix(DOUBLET))
        printed = " | ".join(cells)
        for column, name in TABLE_1_COLUMNS.items():
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _table_4() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_4)
    if len(rows) != 43 or rows[0][1] != "n" or rows[0][4] != "n":
        msg = f"unexpected Table IV layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    points: list[Point] = []
    # The two column pairs in reading order: the left pair (0.22-1.40 um), then the
    # right (1.60-9.80 um).
    for first in (0, 3):
        for cells in rows[1:]:
            wavelength = number(cells[first])
            printed = " | ".join(cells[first : first + 3])
            points.append(
                Point(wavelength, wavelength, number(cells[first + 1]), "", printed)
            )
    return {f"{PREFIX}table4_n_24C.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return Tables I, II and the n of Table IV."""
    _check_table_3()
    return {**_table_2(), **_table_1(), **_table_4()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

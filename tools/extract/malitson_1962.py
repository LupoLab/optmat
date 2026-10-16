# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Malitson (1962): the ordinary index of synthetic sapphire at 24 degrees C.

- **Table II** (``tbl_9f97a99b9fc8``): 46 wavelengths (um, the lines' air wavelengths),
  the spectral source, the observed and the computed ordinary index n_0 at 24 degrees
  C, and the residual (n_0 - n_c) x 1e5. The observed and computed columns are
  extracted; the residual column is not (the difference of the two, except at two rows
  the record's notes discuss). The cells agree with the image.
- **Table III** (``tbl_a823708dead2``): the formula's index at regular wavelengths,
  0.270-5.60 um, in two column pairs (the right pair ends one row early); 61 points.
  The cells agree with the image.

Table I (``tbl_96d9dc95c187``) is the record's formula; this script checks the record's
strings against its cells. Rows are the wavelengths as printed (um).
"""

from _library import Point, number, table, write

DIRECTORY = "Malitson_1962_RefractionDispersionSynthetic"
TABLE_1 = "tbl_96d9dc95c187"
TABLE_2 = "tbl_9f97a99b9fc8"
TABLE_3 = "tbl_a823708dead2"
PREFIX = "Al2O3/malitson_1962__"

TABLE_1_PRINTED = [
    ["\\lambda_1 = 0.06144821", "\\lambda_1^2 = 0.00377588", "A_1 = 1.023798"],
    ["\\lambda_2 = 0.1106997", "\\lambda_2^2 = 0.0122544", "A_2 = 1.058264"],
    ["\\lambda_3 = 17.92656", "\\lambda_3^2 = 321.3616", "A_3 = 5.280792"],
]


def _check_table_1() -> None:
    """Check the record's transcription of Table I against the library cells."""
    found = table(DIRECTORY, TABLE_1)
    if found != TABLE_1_PRINTED:
        msg = f"Table I reads {found}, not {TABLE_1_PRINTED}"
        raise ValueError(msg)


def _table_2() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_2)
    if rows[0][2:4] != ["Observed index", "Computed index"] or len(rows) != 47:
        msg = f"unexpected Table II layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, name in ((2, "observed_o_24C"), (3, "computed_o_24C")):
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def _table_3() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_3)
    if len(rows) != 32 or rows[0][1] != "n" or rows[0][3] != "n":
        msg = f"unexpected Table III layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    points: list[Point] = []
    # The left column pair (0.270-0.870 um), then the right (0.890-5.60 um), whose last
    # row is empty.
    for first in (0, 2):
        for cells in rows[1:]:
            if not cells[first]:
                continue
            wavelength = number(cells[first])
            printed = " | ".join(cells[first : first + 2])
            points.append(
                Point(wavelength, wavelength, number(cells[first + 1]), "", printed)
            )
    return {f"{PREFIX}table3_n_o_24C.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return Table II's observed and computed n_0 and Table III's n."""
    _check_table_1()
    return {**_table_2(), **_table_3()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Malitson (1965): fused-silica index at 20 degrees C (Tables I and II).

Table I (``tbl_b692ffe8e259``, cells only: the library has no CSV) lists 60
wavelengths (um), the spectral source, the index computed from Eq. (1), and four
columns of residuals in units of 1e-6: the arithmetical-mean table of the three prisms
(C-D-G.E.), and Corning code 7940, Dynasil and General Electric type 151 (one prism
each, Gaertner spectrometer). The last row is the average of the absolute residuals,
not a wavelength, and is not extracted.

Structural fix (checked on the table image): the cells give the Corning column's
header as "residual \times 10^9"; the image prints "residual x 10^6" over every residual
column, and the text (p. 1207) gives the residuals "about 12 x 10^-6" in each case. The
extraction takes the factor from the image, 1e-6, for all four columns (the record's
``scale``), and does not read the header factor.

Table II (``tbl_46b93cc0e194``, cells) gives the residuals (units of 1e-6, measured
minus Eq. (1)) of the twelve prisms, four per brand, at nine of the 18 wavelengths used
(Wild spectrometer), then the average residual over all 18 wavelengths and the
nu-value; those two rows are not per wavelength and are left to the record's notes. The
cells agree with the image.

The Table II residuals carry the uncertainty the text states for the Wild spectrometer
(pp. 1207-1208): "+-5 x 10^-6 for 14 visible lines and +-20 x 10^-6 for 4 near infrared
wavelengths". The 18 wavelengths "from 0.4047 to 1.08 u" are Table I's lines in that
range, fourteen up to 0.706519 um and four from 0.852111 um, so the printed rows at
0.8944, 1.014 and 1.083 um get 20 (units of 1e-6) in the sigma column, the others 5.

Rows are the wavelengths as printed (um).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Malitson_1965_InterspecimenComparisonRefractive"
TABLE_1 = "tbl_b692ffe8e259"
TABLE_2 = "tbl_46b93cc0e194"
PREFIX = "SiO2/malitson_1965__"

# Table I column -> dataset id.
TABLE_1_COLUMNS = {
    2: "computed",
    3: "residual_mean",
    4: "residual_corning_7940",
    5: "residual_dynasil",
    6: "residual_ge_151",
}

# Table II column -> dataset id: brand and specimen number.
BRANDS = ("corning_7940", "dynasil", "ge_151")
TABLE_2_COLUMNS = {
    1 + 4 * b + s: f"specimen_{brand}_{s + 1}"
    for b, brand in enumerate(BRANDS)
    for s in range(4)
}


# Table I's fourteen visible lines of the 18 end at 0.706519 um, the four near-infrared
# ones start at 0.852111 um.
NEAR_INFRARED_UM = 0.8


def _table_1() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_1)
    header, body = rows[0], rows[1:]
    if not header[4].startswith("Corning residual"):
        msg = f"unexpected Table I header {header!r}"
        raise ValueError(msg)
    if not body[-1][0].startswith("Average of absolute values"):
        msg = f"unexpected last Table I row {body[-1]!r}"
        raise ValueError(msg)
    body = body[:-1]
    if len(body) != 60:
        msg = f"Table I has {len(body)} wavelengths, not 60"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for column, dataset in TABLE_1_COLUMNS.items():
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{dataset}.csv", []).append(point)
    return files


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    if rows[1][0] != "\\mu" or rows[1][1:] != ["1", "2", "3", "4"] * 3:
        msg = f"unexpected Table II specimen row {rows[1]!r}"
        raise ValueError(msg)
    body = [cells for cells in rows[2:] if cells[0][0].isdigit()]
    if len(body) != 9 or not rows[-2][0].startswith("Av."):
        msg = "Table II does not have nine wavelength rows before the averages"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        wavelength = number(cells[0])
        # The stated uncertainty, in the residuals' unit of 1e-6 (see the docstring).
        sigma = "20" if float(wavelength) > NEAR_INFRARED_UM else "5"
        printed = " | ".join(cells)
        for column, dataset in TABLE_2_COLUMNS.items():
            value = number(cells[column])
            point = Point(wavelength, wavelength, value, sigma, printed)
            files.setdefault(f"{PREFIX}{dataset}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the computed index and every residual column of Tables I and II."""
    return {**_table_1(), **_table_2()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

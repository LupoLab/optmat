# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Leviton & Frey: absolute index of Corning 7980 fused silica, 30-300 K.

All tables were checked cell by cell against the PDF's own text layer (``pdftotext
-layout`` of ``original/source_01``) and against the table images.

- **Table 4** (``tbl_e21403d142ed``): the absolute index n at 16 wavelengths
  ("0.4 microns" ... "2.6 microns") and 13 temperatures (30-300 K), computed from the
  temperature-dependent Sellmeier fit of Table 3 (its caption). One file per
  temperature; rows are the wavelengths in um as printed. The cells agree with the
  image.
- **Table 6** (``tbl_1e2b6a0f092c``): dn/dT (1/K) on the same grid. Structural fix: the
  CSV joins the 13 temperature columns into one space-separated cell, in the header
  and in every row; the image shows 13 columns. The script splits them and checks
  that each row has 13 numbers.
- **Table 2** (``tbl_164512d3b768``, cells only): PMIC minus CHARMS, Delta n (ppm), and
  the spread in dn from PMIC's fits (ppm), at 160-230 K and 295 K (at 632.8 nm, PMIC's
  HeNe wavelength: text p. 2). Structural fix: the cells merge the two rows of the
  first eight temperatures into one cell, and join the last Delta n (-30, at 230 K)
  with the first spread (66, at 160 K) into "-3066"; the image and the text layer
  show "-30" and "66". The script splits that token and checks that 16 numbers
  result. The spread is the file's sigma column.
- **Table 7** (``tbl_835faf4136ad``): Waxler & Cleek's indices, adjusted to vacuum,
  minus CHARMS, at 10 wavelengths and 6 temperatures. One file per temperature.
- **Table 8** (``tbl_bb97e696ab87``): Matsuoka et al.'s indices minus CHARMS, at 3
  wavelengths and 5 temperatures. One file per temperature.

Table 1 (the uncertainty of n), Table 3 (the Sellmeier coefficients; a formula in the
record) and Table 5 (dn/dlambda) are not data files; the record's notes and formula
hold them. Table 3's CSV merges the S2 and S3 columns and the column labels into the
first row's cells; the record's coefficients were transcribed from the image and the
text layer, and this script checks them against the cells (so a change of the library
table is noticed).
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Leviton_2008_TemperatureDependentAbsolute"
TABLE_2 = "tbl_164512d3b768"
TABLE_3 = "tbl_071707d808a6"
TABLE_4 = "tbl_e21403d142ed"
TABLE_6 = "tbl_1e2b6a0f092c"
TABLE_7 = "tbl_835faf4136ad"
TABLE_8 = "tbl_bb97e696ab87"
PREFIX = "SiO2/leviton_2008__"

TEMPERATURES = ("30", "40", "50", "60", "80", "100", "120", "160", "200", "240", "275")
TEMPERATURES += ("295", "300")

# Table 3 as printed (image and text layer), rows T^0 ... T^4, columns S1, S2, S3,
# lambda1, lambda2, lambda3: the strings of the record's formula.
TABLE_3_PRINTED = (
    ("1.10127E+00", "1.78752E-05", "7.93552E-01", "-8.90600E-02", "2.97562E-01"),
    ("-4.94251E-05", "4.76391E-05", "-1.27815E-03", "9.08730E-06", "-8.59578E-04"),
    ("5.27414E-07", "-4.49019E-07", "1.84595E-05", "-6.53638E-08", "6.59069E-06"),
    ("-1.59700E-09", "1.44546E-09", "-9.20275E-08", "7.77072E-11", "-1.09482E-08"),
    ("1.75949E-12", "-1.57223E-12", "1.48829E-10", "6.84605E-14", "7.85145E-13"),
)
LAMBDA_3_PRINTED = ("9.34454E+00", "-7.09788E-03", "1.01968E-04", "-5.07660E-07")
LAMBDA_3_PRINTED += ("8.21348E-10",)
_E_NUMBER = re.compile(r"-?\d\.\d+E[+-]\d+")


def _check_table_3() -> None:
    """Check the record's transcription of Table 3 against the library cells."""
    rows = table(DIRECTORY, TABLE_3)[1:]
    if len(rows) != 5:
        msg = f"Table 3 has {len(rows)} coefficient rows, not 5"
        raise ValueError(msg)
    expected_rows = zip(TABLE_3_PRINTED, LAMBDA_3_PRINTED, strict=True)
    for cells, (printed, last) in zip(rows, expected_rows, strict=True):
        found = [x for cell in cells[1:] for x in _E_NUMBER.findall(cell)]
        if tuple(found) != (*printed, last):
            msg = f"Table 3 row {cells[0]!r} reads {found}, not {(*printed, last)}"
            raise ValueError(msg)


def _wavelength(cell: str) -> str:
    """Return the wavelength in um of a row label such as "0.4 microns"."""
    text = cell.removesuffix(" microns")
    return number(text)


def _grid(rows: list[list[str]], quantity: str) -> dict[str, list[Point]]:
    """Return one file per temperature column of a wavelength x temperature grid."""
    files: dict[str, list[Point]] = {}
    for cells in rows:
        if len(cells) != 1 + len(TEMPERATURES):
            msg = f"row {cells!r} does not have {len(TEMPERATURES)} temperatures"
            raise ValueError(msg)
        wavelength = _wavelength(cells[0])
        printed = " | ".join(cells)
        for temperature, cell in zip(TEMPERATURES, cells[1:], strict=True):
            point = Point(wavelength, wavelength, number(cell), "", printed)
            files.setdefault(f"{PREFIX}{quantity}_{temperature}K.csv", []).append(point)
    return files


def _table_4() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_4)
    header = [cell.removesuffix(" K") for cell in rows[0][1:]]
    if tuple(header) != TEMPERATURES or len(rows) != 17:
        msg = f"unexpected Table 4 layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    return _grid(rows[1:], "n")


def _table_6() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_6)
    # Structural fix: the 13 temperature columns are one space-joined cell (see the
    # docstring); the image shows 13 columns.
    header = rows[0][1].replace(" K", "").split()
    if tuple(header) != TEMPERATURES or len(rows) != 17:
        msg = f"unexpected Table 6 layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    split = [[cells[0], *cells[1].split()] for cells in rows[1:]]
    return _grid(split, "dndt")


def _table_2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_2)
    temperatures = [*rows[0][1].split(), rows[0][2]]
    tokens = rows[1][1].split()
    # Structural fix: "-3066" is the last Delta n (-30) joined to the first spread (66);
    # the image and the text layer print them in separate rows.
    if tokens[7] != "-3066":
        msg = f"unexpected Table 2 cell {rows[1][1]!r}"
        raise ValueError(msg)
    tokens[7:8] = ["-30", "66"]
    if len(temperatures) != 9 or len(tokens) != 16:
        msg = f"unexpected Table 2 layout {rows!r}"
        raise ValueError(msg)
    delta = [*tokens[:8], rows[1][2]]
    spread = [*tokens[8:], rows[2][2]]
    files: dict[str, list[Point]] = {}
    for temperature, value, sigma in zip(temperatures, delta, spread, strict=True):
        printed = f"{temperature} K | {value} | {sigma}"
        # PMIC's wavelength, 0.6328 um (text p. 2); the table prints none.
        point = Point("0.6328", "0.6328", number(value), number(sigma), printed)
        files[f"{PREFIX}pmic_minus_charms_{temperature}K.csv"] = [point]
    return files


def _differences(table_id: str, name: str, count: int) -> dict[str, list[Point]]:
    """Return one file per temperature of Table 7 or 8 (rows are wavelengths, um)."""
    rows = table(DIRECTORY, table_id)
    temperatures = [cell.removesuffix(" K") for cell in rows[0][1:]]
    if not rows[0][0].startswith("wavelength [microns]") or len(rows) != count + 1:
        msg = f"unexpected layout of {table_id}: {rows[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for temperature, cell in zip(temperatures, cells[1:], strict=True):
            point = Point(wavelength, wavelength, number(cell), "", printed)
            files.setdefault(f"{PREFIX}{name}_{temperature}K.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 4 and 6 per temperature, and the differences of Tables 2, 7, 8."""
    _check_table_3()
    return {
        **_table_4(),
        **_table_6(),
        **_table_2(),
        **_differences(TABLE_7, "waxler_cleek_minus_charms", 10),
        **_differences(TABLE_8, "matsuoka_minus_charms", 3),
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

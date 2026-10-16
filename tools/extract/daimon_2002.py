# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Daimon & Masumura (2002): CaF2 index relative to nitrogen, 138-2326 nm.

Table 1 is printed over two pages (``tbl_8ae247914098``, pdf p. 4, 59 wavelengths, and
``tbl_b68ea97710c0``, pdf p. 5, "Continued", 10 wavelengths and the standard-deviation
row); the library publishes both as cells only (no CSV). The columns: the vacuum
wavelength (nm); the relative index in nitrogen (1013.25 hPa) measured at 20 and 25
degrees C, each with its residual from the Table 2 Sellmeier (x 1e-7); the uncertainty
of the index (x 1e-7), "Except lambda" and, below 283 nm only, "Include lambda" (with
the uncertainty of the Pt-Ne line wavelengths); and the relative dn/dt in nitrogen
(x 1e-6 per degree C), measured, residual from Eq. (7) and uncertainty. All 69 rows of
both parts were compared with the images: the cells agree, signs included.

Extracted, one file each (rows are the wavelengths in nm as printed):

- ``n_20C`` and ``n_25C``: the measured relative index. The sigma column is the
  table's uncertainty, "Include lambda" where printed (below 283 nm) and "Except
  lambda" otherwise. The table prints it in units of 1e-7; the sigma column carries
  the printed digits with that factor written out ("9e-7"), since the index itself has
  no column factor.
- ``dndt``: the measured relative dn/dt, in units of 1e-6 per degree C (the record's
  scale), with the table's uncertainty (same units) as sigma.

The residual columns and the standard-deviation row are not data files: they are the
formulas' own checks, reproduced in the record (``check``). Table 2 (the Sellmeier
constants, cells only) and Table 3 (the dn/dt constants) are formulas in the record;
this script checks the record's strings against the cells, so that a change of the
library tables is noticed.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Daimon_2002_HighAccuracyMeasurements"
TABLE_1 = ("tbl_8ae247914098", "tbl_b68ea97710c0")
TABLE_2 = "tbl_1b907e57ed10"
TABLE_3 = "tbl_e48e1a7b603c"
PREFIX = "CaF2/daimon_2002__"

# Table 2 as printed (image and cells): A_i and lambda_i^2 at 20 and 25 degrees C, the
# strings of the record's formulas, with their printed powers of ten.
TABLE_2_PRINTED = (
    (
        ("4.43749998 \\times 10^{-1}", "1.78027854 \\times 10^{-3}"),
        ("4.44930066 \\times 10^{-1}", "7.88536061 \\times 10^{-3}"),
        ("1.50133991 \\times 10^{-1}", "1.24119491 \\times 10^{-2}"),
        ("8.85319946", "2.75228175 \\times 10^3"),
    ),
    (
        ("4.37387571 \\times 10^{-1}", "1.73799328 \\times 10^{-3}"),
        ("4.49211397 \\times 10^{-1}", "7.82718648 \\times 10^{-3}"),
        ("1.52068715 \\times 10^{-1}", "1.24086125 \\times 10^{-2}"),
        ("1.30020420 \\times 10^1", "4.03976541 \\times 10^3"),
    ),
)
TABLE_3_PRINTED = (
    ("D_{0}", "-2.8422 \\times 10^{-5}"),
    ("E_{0}", "4.6988 \\times 10^{-7}"),
    ("\\lambda_{0}^{2}", "1.3459 \\times 10^{-2}"),
)
HEADER = [
    "",
    "Measured",
    "Residual (\\times 10^{-7})",
    "Measured",
    "Residual (\\times 10^{-7})",
]


def _check_formulas() -> None:
    """Check the record's transcription of Tables 2 and 3 against the library cells."""
    rows = json_table(DIRECTORY, TABLE_2)
    body = rows[2:]
    if [r[0] for r in body] != ["1", "2", "3", "4"]:
        msg = f"unexpected Table 2 layout {rows!r}"
        raise ValueError(msg)
    for column, printed in zip((1, 3), TABLE_2_PRINTED, strict=True):
        found = tuple((r[column], r[column + 1]) for r in body)
        if found != printed:
            msg = f"Table 2 column {column} reads {found}, not {printed}"
            raise ValueError(msg)
    found3 = tuple(tuple(r) for r in table(DIRECTORY, TABLE_3))
    if found3 != TABLE_3_PRINTED:
        msg = f"Table 3 reads {found3}, not {TABLE_3_PRINTED}"
        raise ValueError(msg)


def _rows() -> list[list[str]]:
    """Return Table 1's 69 wavelength rows, both parts, in printed order."""
    body: list[list[str]] = []
    for part in TABLE_1:
        rows = json_table(DIRECTORY, part)
        if rows[2][:5] != HEADER or rows[2][5] not in {"Except λ", "Except \\lambda"}:
            msg = f"unexpected Table 1 header in {part}: {rows[2]!r}"
            raise ValueError(msg)
        data = rows[3:]
        if data[-1][0] == "Standard deviation":
            data = data[:-1]
        body += data
    if len(body) != 69:
        msg = f"Table 1 has {len(body)} wavelengths, not 69"
        raise ValueError(msg)
    return body


def extract() -> dict[str, list[Point]]:
    """Return the measured relative index at 20 and 25 degrees C and dn/dt."""
    _check_formulas()
    files: dict[str, list[Point]] = {}
    for cells in _rows():
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        # The uncertainty with the wavelength's where printed (below 283 nm), else the
        # one without; printed in units of 1e-7 (see the docstring).
        uncertainty = number(cells[6] or cells[5])
        sigma = f"{uncertainty}e-7"
        for column, name in ((1, "n_20C"), (3, "n_25C")):
            point = Point(wavelength, wavelength, number(cells[column]), sigma, printed)
            files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
        dndt, dndt_sigma = number(cells[7]), number(cells[9])
        point = Point(wavelength, wavelength, dndt, dndt_sigma, printed)
        files.setdefault(f"{PREFIX}dndt.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

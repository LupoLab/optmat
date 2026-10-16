# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Frey, Leviton & Madison (2006): absolute index of Si and Ge, 30-295 K.

The library copy is a scan (NASA NTRS) with an OCR text layer. Every cell of the
extracted tables was compared with the table images and with the PDF's text layer
(``pdftotext -layout`` of ``original/source_01``, its OCR confusions S/8, O/0, l/1
allowed for): no misread cell. Rows the text layer splits or drops (Si Table 2 at 1.1
um; Ge Table 7 at 1.8, 2.5, 3.0, 4.0, 4.5, 5.0, 5.5 um; Ge Table 9 at 2.5 and 5.0 um)
were read on the images.

The library attaches the captions to the wrong tables (each caption printed below a
figure and a table is read as the next item's): ``tbl_032f91d73970`` carries Figure 1's
caption but is Table 2, and ``tbl_7450bffe6825`` (unlabelled in the library) is Table 4
(dn/dT), whose caption the library gives to ``tbl_e700ce199afe`` (Table 5). The images
and the text layer fix which is which; the table ids below are by content.

- **Table 2** (``tbl_032f91d73970``): Si absolute index at 13 wavelengths ("1.1
  microns" ... "5.5 microns") and 12 temperatures (30-295 K). One file per
  temperature; rows are the wavelengths in um as printed.
- **Table 4** (``tbl_7450bffe6825``): Si dn/dT (1/K) on the same grid.
- **Table 7** (``tbl_f6ae92c9ec46``): Ge absolute index at 9 wavelengths (1.8-5.5 um)
  and the same 12 temperatures.
- **Table 9** (``tbl_b208df1e1ae7``): Ge dn/dT (1/K) on the Table 7 grid.

Tables 1 and 6 (the uncertainty of n) and Tables 3 and 8 (dn/dlambda) are not data
files; the record's notes hold them. Tables 5 and 10 (the Sellmeier coefficients) are
the record's formulas, transcribed from the images; this script checks them against
the library cells (``tbl_e700ce199afe``, ``tbl_3cbb049ce020``; JSON only), so a change
of the library table is noticed.
"""

import tomllib
from typing import Final

from _json_table import json_table
from _library import DATA, Point, number, table, write

DIRECTORY: Final = "Frey_2006_TemperatureDependentRefractive"
TEMPERATURES: Final = (
    "30",
    "40",
    "50",
    "60",
    "70",
    "80",
    "90",
    "100",
    "150",
    "200",
    "250",
    "295",
)
# (library table, material, quantity prefix, number of wavelength rows)
GRIDS: Final = (
    ("tbl_032f91d73970", "Si", "n", 13),
    ("tbl_7450bffe6825", "Si", "dndt", 13),
    ("tbl_f6ae92c9ec46", "Ge", "n", 9),
    ("tbl_b208df1e1ae7", "Ge", "dndt", 9),
)
# (library table, the record's formula id): the Sellmeier coefficient tables.
COEFFICIENTS: Final = (
    ("tbl_e700ce199afe", "sellmeier_Si"),
    ("tbl_3cbb049ce020", "sellmeier_Ge"),
)
TERMS: Final = ("S1", "S2", "S3", "L1", "L2", "L3")


def _wavelength(cell: str) -> str:
    """Return the wavelength in um of a row label such as "1.1 microns"."""
    return number(cell.removesuffix(" microns"))


def _grid(table_id: str, material: str, quantity: str, count: int) -> dict:
    """Return one file per temperature column of a wavelength x temperature table."""
    rows = table(DIRECTORY, table_id)
    header = [cell.removesuffix(" K") for cell in rows[0][1:]]
    if tuple(header) != TEMPERATURES or len(rows) != count + 1:
        msg = f"unexpected layout of {table_id}: {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        wavelength = _wavelength(cells[0])
        printed = " | ".join(cells)
        for temperature, cell in zip(TEMPERATURES, cells[1:], strict=True):
            point = Point(wavelength, wavelength, number(cell), "", printed)
            name = f"{material}/frey_2006__{quantity}_{material}_{temperature}K.csv"
            files.setdefault(name, []).append(point)
    return files


def _check_coefficients() -> None:
    """Check the record's Sellmeier coefficients against the library cells."""
    with (DATA / "sources" / "frey_2006.toml").open("rb") as handle:
        formulas = {f["id"]: f for f in tomllib.load(handle)["formulas"]}
    for table_id, formula_id in COEFFICIENTS:
        rows = json_table(DIRECTORY, table_id)
        body = [cells for cells in rows if cells[0].endswith("term")]
        if len(body) != 5:
            msg = f"{table_id}: {len(body)} coefficient rows, not 5"
            raise ValueError(msg)
        terms = formulas[formula_id]["terms"]
        for column, term in enumerate(TERMS, start=1):
            found = [number(cells[column]) for cells in body]
            if found != terms[term]:
                msg = f"{formula_id} {term}: library {found}, record {terms[term]}"
                raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Tables 2, 4, 7 and 9 per temperature."""
    _check_coefficients()
    files: dict[str, list[Point]] = {}
    for table_id, material, quantity, count in GRIDS:
        files.update(_grid(table_id, material, quantity, count))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

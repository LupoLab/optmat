# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gorshelev et al. (2014): O3 cross sections at selected wavelengths, ~293 K.

Four comparison tables (CSV; each image checked, every cell agrees), each a column per
source, one file per column:

- Table 4 (``tbl_102450bc135d``), the UV, "cm^2 molecule^-1 x 10^-20": Orphal (2003),
  this work, Hearn (1961), BMD (1998), Burrows et al. (1999a), Bogumil et al. (2003) and
  Petersen et al. (2012) "using Hearn/BMD value" (two values per cell, ``964/950``: two
  files, ``petersen_2012_hearn`` and ``petersen_2012_bmd``);
- Table 5 (``tbl_8eb4cf894101``), the minimum near 380 nm, "x 10^-23", each value
  followed by the typical FWHM in nm in parentheses (``4.9 (0.02)``: the value is kept,
  the FWHM stays in ``printed``);
- Table 6 (``tbl_a995b8119c7a``), the visible, "x 10^-21";
- Table 7 (``tbl_6fa5f0d2dae8``), the near infrared, "x 10^-22".

The wavelength column of Tables 4 and 6 is headed "vacuum (air)" and prints both,
``253.73 (253.65)``, or the vacuum value alone: the vacuum value is the abscissa. Table
4's last wavelength cell is ``e 302.24 (302.15)``: the image prints the stray "e" too
(a typesetting remnant), dropped here. Orphal's mean values carry a relative spread,
``1141\N{PLUS-MINUS SIGN}0.9%``: the value is kept, the spread stays in ``printed``.
"-" and empty cells are no value. Table 4's row 248.32 nm is printed a factor of ten
low in every column (an erratum in the source record, not a structural fix).
"""

import re

from _library import Point, table, write

DIRECTORY = "Gorshelev_2014_HighSpectralResolution"
PREFIX = "O3/gorshelev_2014__"

# Each table: its id, the number of header rows, and the dataset of each value column.
TABLES = {
    "tbl_102450bc135d": (
        3,
        (
            "orphal_2003_uv",
            "this_work_uv",
            "hearn_1961_uv",
            "bmd_uv",
            "burrows_1999_uv",
            "bogumil_2003_uv",
            ("petersen_2012_hearn", "petersen_2012_bmd"),
        ),
    ),
    "tbl_8eb4cf894101": (
        1,
        (
            "this_work_380",
            "bmd_380",
            "burrows_1999_380",
            "bogumil_2003_380",
            "burkholder_talukdar_1994_380",
            "axson_2011_380",
        ),
    ),
    "tbl_a995b8119c7a": (
        1,
        (
            "orphal_2003_vis",
            "this_work_vis",
            "bmd_vis",
            "burrows_1999_vis",
            "bogumil_2003_vis",
            "anderson_mauersberger_1992_vis",
            "burkholder_talukdar_1994_vis",
            "el_helou_2005_vis",
        ),
    ),
    "tbl_6fa5f0d2dae8": (
        2,
        (
            "this_work_nir",
            "bmd_nir",
            "burrows_1999_nir",
            "bogumil_2003_nir",
            "anderson_1993_nir",
            "burkholder_talukdar_1994_nir",
            "el_helou_2005_nir",
        ),
    ),
}
_WAVELENGTH = re.compile(r"^(?:e )?([\d.]+)(?: \(([\d.]+)\))?$")
_VALUE = re.compile(r"^([\d.]+)(?:\s*\N{PLUS-MINUS SIGN}\s*[\d.]+\s*%| \([\d.]+\))?$")
_NONE = {"", "-", "\N{EN DASH}"}


def _value(cell: str) -> str:
    match = _VALUE.match(cell)
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return match.group(1)


def extract() -> dict[str, list[Point]]:
    """Return one file per source column of Tables 4-7."""
    files: dict[str, list[Point]] = {}
    for table_id, (header_rows, columns) in TABLES.items():
        for raw in table(DIRECTORY, table_id)[header_rows:]:
            cells = [c.strip() for c in raw]
            match = _WAVELENGTH.match(cells[0])
            if match is None:
                msg = f"unexpected wavelength {cells[0]!r}"
                raise ValueError(msg)
            wavelength = match.group(1)
            printed = " | ".join(cells)
            for column, cell in zip(columns, cells[1:], strict=True):
                if cell in _NONE:
                    continue
                if isinstance(column, tuple):  # Petersen: "Hearn/BMD"
                    names, values = column, cell.split("/")
                else:
                    names, values = (column,), [cell]
                for name, text in zip(names, values, strict=True):
                    point = Point(wavelength, wavelength, _value(text), "", printed)
                    files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

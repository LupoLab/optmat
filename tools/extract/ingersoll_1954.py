# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ingersoll and Liebenberg (1954): Verdet constants of He, Ne, H2, D2, N2, O2.

- Table I (``tbl_794626e36b6c``; image checked): Verdet constants (microminutes per
  oersted-centimeter, gas at 76 cm and 0 C) read from smooth curves at 13 filter
  wavelengths, 3635-9875 A. Helium has no value at 9875 A (printed ``...``), skipped.
  n-Butane is not curated.
- Table III (``tbl_b0f0c3f0522a``; image checked): the same smooth-curve values at
  4360, 5460, 5780, 5893 A, with de Mallemann's values in parentheses (``(M 1.80)``)
  and Bizette and Tsai's for helium (``(BT 0.45)``). Both the authors' values and the
  parenthesised ones are extracted; the first row prints ``4360 A``.

Table II (dispersion relative to 4000 A) is derived from Table I; Table IV (pressure
effect, individual determinations) is not extracted (see the source record).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Ingersoll_1954_FaradayEffectGases"
TABLE_I = "tbl_794626e36b6c"
TABLE_III = "tbl_b0f0c3f0522a"
GASES = {1: "He", 2: "Ne", 3: "H2", 4: "D2", 5: "N2", 6: "O2"}
_WITH_OTHER = re.compile(r"^\s*([\d.]+)\s*(?:\((M|BT)\s*([\d.]+)\))?\s*$")
OTHERS = {"M": "mallemann", "BT": "bizette_tsai"}


def _table_i() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        wavelength = number(cells[0])
        for column, gas in GASES.items():
            if cells[column].strip() == "...":
                continue
            files.setdefault(f"{gas}/ingersoll_1954__verdet_{gas}.csv", []).append(
                Point(
                    wavelength, wavelength, number(cells[column]), "", " | ".join(cells)
                )
            )
    return files


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_III)[1:]:
        wavelength = number(cells[0].removesuffix(" A"))
        printed = " | ".join(cells)
        for column, gas in GASES.items():
            match = _WITH_OTHER.match(cells[column])
            if match is None:
                msg = f"unexpected Table III cell {cells[column]!r}"
                raise ValueError(msg)
            prefix = f"{gas}/ingersoll_1954__"
            files.setdefault(f"{prefix}verdet_lines_{gas}.csv", []).append(
                Point(wavelength, wavelength, match.group(1), "", printed)
            )
            if match.group(2):
                other = OTHERS[match.group(2)]
                files.setdefault(f"{prefix}{other}_{gas}.csv", []).append(
                    Point(wavelength, wavelength, match.group(3), "", printed)
                )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I and III as data files."""
    return _table_i() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract W.-C. Liu (1997): helium dynamic polarizability, 0-58 eV.

- Table I (``tbl_c824036c5cc9``): the dynamic polarizability alpha(omega) of ground-
  state helium in atomic units at photon energies omega = 0 to 0.75 a.u. (below the
  first excitation threshold).
- Table II (``tbl_936a0db2e146``): above the ionisation threshold, the complex
  polarizability printed as ``-0.758+1.911i`` (split here into its parts: only the
  imaginary part is extracted, as the format has no quantity for a negative real part
  above threshold), the photoionisation cross section (present work) and Samson et
  al.'s measured cross section, at photon energies in eV.
- Table IV (``tbl_5f861154394b``): the q = 1 row, |chi_1| in esu (a polarizability
  volume in cm^3) at "lambda ~ 46 um (low frequency limit)", 1064, 532 and 248 nm. All
  four are printed with ``x 10^-25``; the script checks that exponent.

The table images agree with the extracted cells. The higher-order susceptibilities
(Tables III and IV, q >= 3) are nonlinear and not extracted here.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Liu_1997_HighOrderNonlinear"
TABLE_I = "tbl_c824036c5cc9"
TABLE_II = "tbl_936a0db2e146"
TABLE_IV = "tbl_5f861154394b"
PREFIX = "He/liu_1997__"
# Table IV columns -> wavelength as printed (um for the first, nm for the others).
TABLE_IV_WAVELENGTHS = {1: "46000", 2: "1064", 3: "532", 4: "248"}
_COMPLEX = re.compile(r"^([-+]?[\d.]+)([-+][\d.]+)i$")
_POWER = re.compile(r"^\{([\d.]+)\}\s*\\times\s*\{10\}\^\{-\{(\d+)\}\}$")


def _table_i() -> list[Point]:
    rows = table(DIRECTORY, TABLE_I)[2:]  # two header rows: names and symbols
    return [
        Point(number(c[0]), number(c[0]), number(c[1]), "", " | ".join(c)) for c in rows
    ]


def _table_ii() -> dict[str, list[Point]]:
    imaginary, cross_section, samson = [], [], []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        energy, printed = number(cells[0]), " | ".join(cells)
        parts = _COMPLEX.match(cells[1].strip().replace("\N{MINUS SIGN}", "-"))
        if parts is None:
            msg = f"cannot read the complex polarizability {cells[1]!r}"
            raise ValueError(msg)
        imaginary.append(Point(energy, energy, number(parts.group(2)), "", printed))
        cross_section.append(Point(energy, energy, number(cells[2]), "", printed))
        samson.append(Point(energy, energy, number(cells[3]), "", printed))
    return {
        f"{PREFIX}alpha_imaginary.csv": imaginary,
        f"{PREFIX}cross_section.csv": cross_section,
        f"{PREFIX}samson_cross_section.csv": samson,
    }


def _table_iv() -> list[Point]:
    _, first, *_ = table(DIRECTORY, TABLE_IV)
    if first[0].strip() != "1":
        msg = f"Table IV's first row is not q = 1: {first!r}"
        raise ValueError(msg)
    points = []
    for column, wavelength in TABLE_IV_WAVELENGTHS.items():
        match = _POWER.match(first[column].strip())
        if match is None or match.group(2) != "25":
            msg = f"unexpected chi_1 cell {first[column]!r}"
            raise ValueError(msg)
        points.append(
            Point(wavelength, wavelength, number(match.group(1)), "", first[column])
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return Tables I and II (imaginary part, cross sections) and IV (q = 1)."""
    files = {f"{PREFIX}alpha.csv": _table_i()}
    files |= _table_ii()
    files[f"{PREFIX}chi1.csv"] = _table_iv()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

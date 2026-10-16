# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Franta, Necas & Zajickova (2013): the atomic densities of Table 1.

Table 1 (p. 436) gives, per material, the factor U, the density parameter N_a (eV^2),
the atomic density and the electron density (1/m^3) and the mass density (kg/m^3). Only
the atomic density is extracted (``atom_density``); the rest are in the record's notes.

Structural fix: the library's cells lost the superscripts of the powers of ten, so
4.994 (times sign) 10^28 reads as 4.994, the times sign, then "1028". The table image
shows every density as a mantissa times 10 to a two-digit power; the script restores
"4.994e28" (the mantissa, then the digits after "10" as the exponent) and keeps the
raw cells in ``printed``.
"""

import re

from _library import Point, table, write

DIRECTORY = "Franta_2013_ApplicationThomasReiche"
TABLE = "tbl_7c9db7f2834e"
MATERIALS = {
    "c-Si": "Si",
    "c-C (diamond)": "C_diamond",
    "c-C (graphite)": "C_graphite",
    "a-SiO2": "SiO2",
    "a-Si:H": "a-Si_H",
    "a-C:H (DLC)": "a-C_H",
}
_POWER = re.compile(r"^(\d+\.\d+)\N{MULTIPLICATION SIGN}10(\d{2})$")


def _power(cell: str) -> str:
    """Return mantissa, times sign, "10dd" (a lost superscript) as "m.mmmedd"."""
    match = _POWER.match(cell.strip())
    if match is None:
        msg = f"not a mantissa times a power of ten: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return one atomic-density dataset per material."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        material = MATERIALS[cells[0]]
        files[f"{material}/franta_2013__atom_density_{material}.csv"] = [
            Point("1", "", _power(cells[3]), "", " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

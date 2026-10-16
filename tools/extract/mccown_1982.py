# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract McCown, Ediger & Eden (1982): two-photon ionisation of Xe at 193 nm.

Table I (p. 3322, library table ``tbl_0641e78d7d33``, JSON cells only) prints the
coefficient alpha of the two-photon ionisation of xenon and propane at 193 nm in
cm4 W-1, defined by the rate equation dn/dt = alpha [Xe] I^2 / hbar omega (Eq. 9).
Each cell holds the column's name and its number, ``Present work4 \times 10^{-32}``;
the script splits them and checks the power of ten against the record's scale
(:data:`XENON`). Xenon's row: the present measurement and three values the authors
reproduce (Bischel et al., Ref. 11; Hodges, Lee and Moseley, Ref. 12, a lower limit;
McGuire, Ref. 13, a calculation for circular polarisation), one file each. Propane is
not a curated material and is not read. The cells were compared with the table image
and the page (2026-10-08): no misread cell.

**Checks** (:func:`checks`): sigma_2 = alpha hbar omega, the generalised cross section
on the photon flux that alpha's definition implies (the rate alpha I^2 / hbar omega
equals sigma_2 Phi^2 with Phi = I / hbar omega), at the ArF line centre the paper prints
(two photons 103 466 cm-1, p. 3321), must round to the 4e-50 cm4 s that Santra and
Greene quote for this measurement (santra_2004, p. 7 of the preprint).
"""

import math
import re

from _json_table import json_table
from _library import Point, number, write

from optmat.constants import HBAR, C

DIRECTORY = "McCown_1982_TwoPhotonIonization"
TABLE_I = "tbl_0641e78d7d33"
WAVELENGTH_NM = "193"

# Table I's xenon row, by column: (the column's printed name, the row id, the file,
# the power of ten the record's scale holds).
XENON = (
    ("Present work", "present", "Xe/mccown_1982__alpha_Xe.csv", -32),
    ("Bischel et al.^{a}", "bischel", "Xe/mccown_1982__bischel_Xe.csv", -34),
    ("Hodges, Lee, and Moseley^{b}", "hodges", "Xe/mccown_1982__hodges_Xe.csv", -34),
    ("McGuire^{c}", "mcguire", "Xe/mccown_1982__mcguire_Xe.csv", -32),
)
_CELL = re.compile(r"^(?P<name>.*?)(?P<number>[\d.]+ \\times 10\^\{(?P<power>-\d+)\})$")

# The ArF line centre: twice its energy is printed as 103466 +- 70 cm-1 (p. 3321).
TWO_PHOTON_CM1 = 103466.0
# Santra & Greene's quotation of this measurement, 'about 4 x 10^-50 cm^4 s' (p. 7).
SANTRA_SIGMA2 = 4e-50


def extract() -> dict[str, list[Point]]:
    """Return the xenon coefficients of Table I, one file per column."""
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0][0] != "Xenon" or grid[2][0] != "Propane":
        msg = f"Table I's layout has changed: {grid[0][0]!r}, {grid[2][0]!r}"
        raise ValueError(msg)
    files = {}
    for column, (name, row, file, power) in enumerate(XENON):
        cell = grid[1][column]
        match = _CELL.match(cell)
        if match is None or match["name"] != name or int(match["power"]) != power:
            msg = f"Table I, xenon column {column}: unexpected cell {cell!r}"
            raise ValueError(msg)
        value = number(match["number"])
        files[file] = [Point(row, WAVELENGTH_NM, value, "", cell)]
    return files


def sigma2() -> float:
    """Return alpha hbar omega at the ArF line centre, in cm4 s."""
    (point,) = extract()["Xe/mccown_1982__alpha_Xe.csv"]
    alpha_cm4_per_w = float(point.value) * 1e-32
    photon_j = HBAR * 2 * math.pi * C * (TWO_PHOTON_CM1 / 2) * 100
    return alpha_cm4_per_w * photon_j


def checks() -> dict[str, str]:
    """Raise unless alpha hbar omega rounds to Santra and Greene's 4e-50 cm4 s."""
    value = sigma2()
    if f"{value:.0e}" != f"{SANTRA_SIGMA2:.0e}":
        msg = f"alpha hbar omega = {value:.3e} cm4 s, not about {SANTRA_SIGMA2:.0e}"
        raise ValueError(msg)
    return {"alpha hbar omega": f"{value:.3e} cm4 s (santra_2004: about 4e-50)"}


def main() -> None:
    """Write the data files and print the check."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

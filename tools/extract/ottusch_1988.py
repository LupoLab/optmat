# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ottusch and Rockwell (1988): Raman gains of H2, D2 and CH4 at 532 nm.

Table I (``tbl_6fef6b554123``, p. 2077; CSV checked against the image) lists the
spectroscopic parameters of the gain calculation, Eq. (4): the Stokes wavelength, the
population difference, the FWHM linewidth law (the record's formulas) and the cross
section dsigma/dOmega (cm^2/sr, ``0.81 \times 10^{-30}``: the factor 10^-30 is the
dataset scale). Only the cross sections are written here.

Table II (``tbl_963e1d09eaae``, p. 2079; JSON cells, image checked) gives, per gas and
density (amagat), the calculated gain, this work's measured gain ``value ± error`` and a
previous measurement, in cm/GW at a 532 nm pump. STRUCTURAL FIX: the cells read the
methane row as ``CH_{3} \nu_{1}``; the image prints CH4 nu1.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Ottusch_1988_MeasurementRamanGain"
TABLE_I = "tbl_6fef6b554123"
TABLE_II = "tbl_963e1d09eaae"
PUMP = "532"  # nm, the pump of every value
TABLE_I_GASES = {r"\text{H}_2 Q(1)": "H2", r"\text{D}_2 Q(2)": "D2"}
TABLE_I_GASES |= {r"\text{CH}_4 \nu_1": "CH4"}
# The image prints CH4 nu1 where the cells read CH_{3} (see the module docstring).
TABLE_II_GASES = {"H_{2} Q(1)": "H2", "D_{2} Q(2)": "D2", r"CH_{3} \nu_{1}": "CH4"}
_FACTOR = re.compile(r"^([\d.]+) \\times 10\^\{-30\}$")
_MARK = re.compile(r"\^\{[ab]\}$")


def _table_i() -> dict[str, list[Point]]:
    """Return the cross sections of Table I."""
    files = {}
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        gas = TABLE_I_GASES[cells[0]]
        match = _FACTOR.match(cells[4].strip())
        if match is None:
            msg = f"unexpected cross section {cells[4]!r}"
            raise ValueError(msg)
        point = Point(PUMP, PUMP, number(match.group(1)), "", " | ".join(cells))
        files[f"{gas}/ottusch_1988__dsigma_{gas}.csv"] = [point]
    return files


def _table_ii() -> dict[str, list[Point]]:
    """Return the calculated, measured and previous gains of Table II."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_II)[2:5]:
        gas = TABLE_II_GASES[cells[0]]
        density, printed = cells[1], " | ".join(cells)
        row = f"{density} amagat"
        stem = f"{gas}/ottusch_1988__"
        files[f"{stem}gain_calc_{gas}.csv"] = [
            Point(row, density, number(cells[2]), "", printed)
        ]
        value, sigma = cells[3].split("±")
        files[f"{stem}gain_{gas}.csv"] = [
            Point(row, density, number(value), number(sigma), printed)
        ]
        previous = _MARK.sub("", cells[4].strip())
        if previous:
            parts = previous.split(r"\pm")
            sigma = number(parts[1]) if len(parts) > 1 else ""
            files[f"{stem}gain_previous_{gas}.csv"] = [
                Point(row, density, number(parts[0]), sigma, printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables I and II as data files."""
    return _table_i() | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

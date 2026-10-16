# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hanna, Pointer and Pratt (1986): Table I, the calculated SRS parameters.

Table I (``tbl_555a4008bd8a``, p. 333; no CSV, the JSON cells match the image) gives,
for H2 Q(1), D2 Q(2) and the CH4 Q branch, the Raman shift, and at the pumps 1.06 and
0.53 um the Stokes wavelength, the cross section dsigma/dOmega (m^2/sr), the
steady-state gain g_R (m/W) and, per molecule, the dephasing time T2 (ps), 'both
calculated at 30 atm pressure'. The cross sections are the numbers as printed,
``1.2x10-36`` written ``1.2e-36`` (a change of notation only); the D2 ones carry an
asterisk, 'Scaled from the measured value at 694nm reported by Devir [22]', kept in the
printed column. T2 is printed once per molecule (a merged cell) and is written once,
at 30 atm.

Rows: each molecule's first row (the 1.06 um pump) names the gas, its second the
transition and its third the shift and the 0.53 um pump.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hanna_1986_StimulatedRamanScattering"
TABLE_I = "tbl_555a4008bd8a"
# Library row of the 1.06 um line of each molecule -> (material, transition label).
BLOCKS = {
    1: ("H2", "H_{2}", "Q(1)"),
    4: ("D2", "D_{2}", "Q(2)"),
    7: ("CH4", "CH_{4}", "Q"),
}
PUMPS = ("1.06", "0.53")  # um, the pump wavelengths as printed
_POWER = re.compile(r"^([\d.]+)x10-(\d+)\*?$")


def _power(cell: str) -> str:
    """Return ``1.2x10-36`` (a trailing asterisk dropped) as ``1.2e-36``."""
    match = _POWER.match(cell.strip())
    if not match:
        msg = f"not a power of ten: {cell!r}"
        raise ValueError(msg)
    return number(f"{match.group(1)}e-{match.group(2)}")


def extract() -> dict[str, list[Point]]:
    """Return Table I as data files."""
    grid = json_table(DIRECTORY, TABLE_I)
    files: dict[str, list[Point]] = {}
    for row, (gas, printed_gas, transition) in BLOCKS.items():
        first, label, third = grid[row], grid[row + 1], grid[row + 2]
        if (first[0], label[0]) != (printed_gas, transition):
            msg = f"row {row} is {first[0]!r}/{label[0]!r}"
            raise ValueError(msg)
        shift = re.fullmatch(r"Shift (\d+)cm-1", third[0])
        if shift is None:
            msg = f"row {row + 2} has no shift: {third[0]!r}"
            raise ValueError(msg)
        stem = f"{gas}/hanna_1986__"
        files[f"{stem}shift_{gas}.csv"] = [
            Point(transition, transition, shift.group(1), "", third[0])
        ]
        for pump, cells in zip(PUMPS, (first, third), strict=True):
            if not cells[1].startswith(f"{pump}μm"):
                msg = f"pump cell {cells[1]!r} is not {pump} um"
                raise ValueError(msg)
            printed = " | ".join(cells[1:5])
            files.setdefault(f"{stem}dsigma_{gas}.csv", []).append(
                Point(f"{pump} um pump", pump, _power(cells[3]), "", printed)
            )
            files.setdefault(f"{stem}gain_{gas}.csv", []).append(
                Point(f"{pump} um pump", pump, _power(cells[4]), "", printed)
            )
        files[f"{stem}T2_{gas}.csv"] = [
            Point("30 atm", "30", number(first[5]), "", f"{first[0]} | {first[5]}")
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

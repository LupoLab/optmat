# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rosasco, Lempert, Hurst & Fein (1983), Chem. Phys. Lett. 97, 435.

The unlabelled table on p. 438 (``tbl_ee3e296d7c2c``, JSON cells only; its caption,
"Table 1. Line-broadening and mixing-coefficients at 298 K", is lost in the extraction)
was checked against the page rendered at 400 dpi. It gives, per J = 0-18, for N2 and
CO: the line-broadening coefficient gamma_j (cm-1 atm-1) with its error in parentheses
in units of the last decimal ("0.065(15)" is 0.065 +- 0.015), and the line-mixing
coefficient Y_j (atm-1). The extraction misreads the heads: "\aleph_2" is N2 and
"y_f"/"Y_f" are gamma_j and Y_j (image). Only gamma_j is written here (one file per
gas); Y_j has no quantity and is in the record's formulas.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rosasco_1983_LineInterferenceEffects"
TABLE = "tbl_ee3e296d7c2c"
GASES = {1: "N2", 3: "CO"}
_CELL = re.compile(r"^([\d.]+)\((\d+)\)$")


def _last_decimal(value: str, digits: str) -> str:
    """Return an uncertainty printed in units of the value's last decimal."""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) * 10**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the N2 and CO line-broadening coefficients."""
    grid = json_table(DIRECTORY, TABLE)
    # The heads as extracted (N2 read as \aleph_2, gamma_j as y_f; image checked).
    expected = [
        ["J", "\\aleph_2", "", "CO", ""],
        [
            "",
            "y_f (cm^{-1} atm^{-1})",
            "Y_f (atm^{-1})",
            "y_f (cm^{-1} atm^{-1})",
            "Y_f (atm^{-1})",
        ],
    ]
    if grid[:2] != expected:
        msg = f"unexpected table heads {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        label = f"Q({int(cells[0])})"
        for column, gas in GASES.items():
            match = _CELL.match(cells[column])
            if match is None:
                msg = f"not a 'value(n)' cell: {cells[column]!r}"
                raise ValueError(msg)
            value = number(match.group(1))
            files.setdefault(f"{gas}/rosasco_1983__table1_{gas}.csv", []).append(
                Point(
                    label,
                    label,
                    value,
                    _last_decimal(value, match.group(2)),
                    " | ".join(cells),
                )
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

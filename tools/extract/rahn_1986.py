# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rahn & Palmer (1986), J. Opt. Soc. Am. B 3, 1164: N2 Q-branch widths.

Table 1 (``tbl_7d46e56fb640``, JSON cells only; checked against the table image and the
page, p. 1167): the collisional self-broadening coefficients 2 gamma_j (FWHM, cm-1
atm-1) of the N2 Q(J) lines at 295, 500, 730, 1000 and 1500 K, each with its 1-sigma
error in parentheses in units of the last decimal ("0.0851 (46)" is 0.0851 +- 0.0046).
One file per temperature, one row per J, written as the Q(J) line. The 295 K column is
Rosasco et al.'s (note a, "From Ref. 6"). A dash ("-") marks a line not measured at
that temperature and a blank cell a J beyond the measured range; neither is a value.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rahn_1986_StudiesNitrogenSelf"
TABLE_1 = "tbl_7d46e56fb640"
COLUMNS = {1: "295K", 2: "500K", 3: "730K", 4: "1000K", 5: "1500K"}
HEADER = {
    1: "295\\text{K}^{a}",
    2: "500\\text{K}",
    3: "730\\text{K}",
    4: "1000\\text{K}",
    5: "1500\\text{K}",
}
_CELL = re.compile(r"^([\d.]+)\s*\((\d+)\)$")


def _last_decimal(value: str, digits: str) -> str:
    """Return an uncertainty printed in units of the value's last decimal."""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) * 10**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file of Q(J) coefficients per temperature."""
    grid = json_table(DIRECTORY, TABLE_1)
    if grid[0][:2] != ["J", "2\\gamma_{j} (cm^{-1} atm^{-1})"]:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    for column, text in HEADER.items():
        if grid[1][column] != text:
            msg = f"unexpected Table 1 column {column}: {grid[1][column]!r}"
            raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        label = f"Q({int(cells[0])})"
        for column, name in COLUMNS.items():
            cell = cells[column]
            if cell in {"", "-"}:  # not measured at this temperature
                continue
            match = _CELL.match(cell)
            if match is None:
                msg = f"not a 'value (n)' cell: {cell!r}"
                raise ValueError(msg)
            value = number(match.group(1))
            files.setdefault(f"N2/rahn_1986__table1_{name}.csv", []).append(
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

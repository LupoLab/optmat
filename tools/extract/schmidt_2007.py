# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Schmidt et al. (2007): the theoretical A_eps they used (Table I).

Table I (``tbl_db1498b084c8``, cells file only) lists the inputs of the uncertainty
budget; only the theoretical molar polarizability A_eps = 0.517 254 19 (10) cm^3/mol
(row 2, Lach et al. 2004) is a dispersion datum and is extracted, as a reproduction.
The digit groups are printed with spaces and the uncertainty after a space; both are
removed. The image agrees with the cells. The measured A_eps is text-only (source
record).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Schmidt_2007_PolarizabilityHeliumGas"
TABLE_I = "tbl_db1498b084c8"
_VALUE = re.compile(r"^(-?[\d.]+)\((\d+)\)$")


def extract() -> dict[str, list[Point]]:
    """Return the theoretical A_eps of Table I."""
    cells = json_table(DIRECTORY, TABLE_I)[2]
    assert cells[0].startswith("A_{\\varepsilon}"), cells
    match = _VALUE.match(cells[1].replace(" ", ""))
    assert match is not None, cells
    value, digits = match.groups()
    decimals = len(value.split(".")[1])
    sigma = f"{int(digits) * 10.0**-decimals:.{decimals}f}"
    point = Point("static", "inf", number(value), sigma, " | ".join(cells))
    return {"He/schmidt_2007__lach_2004.csv": [point]}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

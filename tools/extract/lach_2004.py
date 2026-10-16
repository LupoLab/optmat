# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lach, Jeziorski & Szalewicz (2004): the static polarizability of 4He.

Table IV (``tbl_81cac1358e95``) lists the contributions to the polarizability "in
atomic units x 10^-6". Only the "Total" row, 1 383 760.79(23), is extracted (the
dataset's scale is 1e-6); its cell is intact in the CSV. The other rows lost their
minus signs and parentheses in the CSV (unmapped glyphs); they are read from the table
image into the source record's notes, not extracted. Tables I-III (convergence of the
Bethe logarithm and its field derivative) hold no optical data.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Lach_2004_RadiativeCorrectionsPolarizability"
TABLE_IV = "tbl_81cac1358e95"
_VALUE = re.compile(r"^([\d.]+)\((\d+)\)$")


def extract() -> dict[str, list[Point]]:
    """Return the total static polarizability."""
    (total,) = [c for c in table(DIRECTORY, TABLE_IV) if c[0].strip() == "Total"]
    match = _VALUE.match(total[1].replace(" ", ""))
    assert match is not None, total
    value, digits = match.groups()
    decimals = len(value.split(".")[1])
    sigma = f"{int(digits) * 10.0**-decimals:.{decimals}f}"
    point = Point("static", "inf", number(value), sigma, " | ".join(total))
    return {"He/lach_2004__alpha_static.csv": [point]}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

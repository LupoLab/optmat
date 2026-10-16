# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hölsch et al. (2019): the ionisation energy of ortho-H2.

Table II (``tbl_7f08590ec4cc``, image checked, p. 4): the energy intervals (cm^-1)
that link X (v = 0, N = 1) of ortho-H2 to X+ (v+ = 0, N+ = 1) of ortho-H2+ and give
the dissociation energy. Only row (4), E_I^ortho(H2) = (1) + (2) + (3), is a
quantity of the vocabulary; the other intervals (and D0 in row (9)) are quoted in the
record's notes. The cell prints the value with digit-group spaces and its uncertainty
in units of the last digit, "124 357.238 003(11)": the spaces are dropped and the
parenthesised digits scaled to the last decimal.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Holsch_2019_BenchmarkingTheoryImproved"
TABLE_2 = "tbl_7f08590ec4cc"
_VALUE = re.compile(r"^([\d .]+)\((\d+)\)$")


def value_and_sigma(cell: str) -> tuple[str, str]:
    """Return the value and the absolute uncertainty of a "124 357.238 003(11)" cell."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"not a 'value(uncertainty)' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1).replace(" ", ""))
    decimals = len(value.partition(".")[2])
    return value, f"{int(match.group(2)) * 10**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the ionisation energy of ortho-H2 (Table II, row (4))."""
    rows = {cells[0].strip(): cells for cells in table(DIRECTORY, TABLE_2)[1:]}
    cells = rows["(4)"]
    if not cells[1].startswith("E_{\\text{I}}^{\\text{ortho}}"):
        msg = f"Table II row (4) is not E_I^ortho: {cells!r}"
        raise ValueError(msg)
    value, sigma = value_and_sigma(cells[2])
    return {
        "H2/holsch_2019__ie_ortho_H2.csv": [
            Point("(4)", "ortho-H2", value, sigma, " | ".join(cells))
        ]
    }


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

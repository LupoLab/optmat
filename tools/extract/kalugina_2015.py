# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kalugina & Thakkar (2015): CCSD Cauchy moments of ozone (Table 4).

Table 4 (``tbl_fa5727ebe567``, checked against the image) gives S(0), S(-2), S(-4) and
S(-6) in atomic units from CCSD with aug-cc-pVnZ basis sets, n = 2-6. Only the
largest basis, n = 6, is extracted (the basis-set limit is reached for S(-2k),
k = 1-3, p. 2941); the orders k are read from the column headers 'S(k)'.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Kalugina_2015_DipolePolarisabilityDipole"
TABLE_4 = "tbl_fa5727ebe567"
_ORDER = re.compile(r"S\(\s*([-\N{MINUS SIGN}]?\s*\d+)\s*\)")


def extract() -> dict[str, list[Point]]:
    """Return the n = 6 CCSD Cauchy moments."""
    header, *rows = table(DIRECTORY, TABLE_4)
    (largest,) = (cells for cells in rows if cells[0].strip() == "6")
    points = []
    for column in range(1, len(header)):
        match = _ORDER.search(header[column])
        if match is None:
            msg = f"no order in header {header[column]!r}"
            raise ValueError(msg)
        k = number(match.group(1).replace(" ", ""))
        points.append(Point(k, k, number(largest[column]), "", " | ".join(largest)))
    return {"O3/kalugina_2015__ccsd_cauchy.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

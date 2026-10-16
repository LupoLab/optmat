# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mizrahi & Shelton (1985, PRL): chi(3) ratios of Ar, N2, O2 to He (Table I).

Table I (``tbl_1008a085d5aa``, CSV; image checked) gives chi(3)(omega)/chi(3)_He(omega)
for ESHG at 19 lines (air wavelength and vacuum wavenumber printed; the wavenumber
cells carry a thin space, "14 282", removed here), with 'the uncertainty of the last
digits ... in parenthesis': "29.67(25)" is 29.67 +- 0.25. Four Ar cells carry the mark
"a" ('obtained by first measuring against N2 used as a transfer standard'): the mark is
stripped and the record names those rows; all O2 values carry the mark in the column
heading. O2 has no value at 457.9 nm ("..."). One file per gas,
``ratio_<gas>_over_He``; row id: the printed air wavelength; abscissa: the wavenumber.
The printed Ar value at 600.0 nm, 20.99(16), is corrected by an erratum in the record
(to 30.99; the image shows 20.99).
"""

import re
from decimal import Decimal

from _library import Point, number, table, write

DIRECTORY = "Mizrahi_1985_DispersionNonlinearSusceptibilities"
TABLE = "tbl_1008a085d5aa"
COLUMNS = {2: "Ar", 3: "N2", 4: "O2"}
_PARENTHESES = re.compile(r"^(\d+\.(\d+))\((\d+)\)$")


def parenthesised(cell: str) -> tuple[str, str]:
    """Return the value and sigma of a cell printed as "29.67(25)"."""
    match = _PARENTHESES.match(cell.strip())
    if match is None:
        msg = f"not a value(uncertainty) cell: {cell!r}"
        raise ValueError(msg)
    value, decimals, digits = match.groups()
    return value, f"{Decimal(digits).scaleb(-len(decimals)):f}"


def extract() -> dict[str, list[Point]]:
    """Return one ratio file per gas."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        label = number(cells[0])
        wavenumber = number(cells[1].replace(" ", ""))
        for column, gas in COLUMNS.items():
            cell = cells[column]
            if cell == "...":
                continue
            value, sigma = parenthesised(cell.removesuffix("a"))
            files.setdefault(
                f"{gas}/mizrahi_1985_arn2o2__ratio_{gas}_over_He.csv", []
            ).append(Point(label, wavenumber, value, sigma, " | ".join(cells)))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

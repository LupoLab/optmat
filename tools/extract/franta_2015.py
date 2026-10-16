# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Franta, Necas & Ohlidal (2015, Appl. Opt.): the hafnia band gap of Table 2.

Table 2 (p. 9117) lists the UDM's interband parameters as "value ± sigma"; the band gap
E_g is the only one with a quantity (``band_gap``). The other parameters are the
record's formula coefficients. The cells agree with the table image.
"""

from _library import Point, number, table, write

DIRECTORY = "Franta_2015_UniversalDispersionModel"
TABLE = "tbl_ad5b7faa65e7"


def extract() -> dict[str, list[Point]]:
    """Return the band gap with its standard deviation."""
    rows = {cells[0]: cells for cells in table(DIRECTORY, TABLE)[1:]}
    cells = rows[r"E_{\text{g}}"]
    value, sigma = (part.strip() for part in cells[1].split("\N{PLUS-MINUS SIGN}"))
    point = Point("1", "", number(value), number(sigma), " | ".join(cells))
    return {"HfO2/franta_2015__band_gap.csv": [point]}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

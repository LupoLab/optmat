# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Farrow, Lucht & Rahn (1987): the N2 reference chi_nr of Table 2.

Table 2 (``tbl_4c51c969c816``; the cells agree with the table image) lists nonresonant
susceptibilities in 10^-18 (cm^3/erg)/amagat of propane, n-butane, water vapour (not
curated species) and N2. Only N2's "Electronic chi_nr(3)" cell, "8.5^c pm 2" with a
LaTeX pm (note c: Ref. 9, Farrow and Rahn 1985), is extracted; its footnote letter is
stripped.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Farrow_1987_MeasurementsNonresonantThird"
TABLE = "tbl_4c51c969c816"
PUMP = "532"  # nm, the frequency-doubled single-mode Nd:YAG pump (p. 1241)
_CELL = re.compile(r"^(?P<value>[\d.]+)\^[a-z] \\pm (?P<sigma>[\d.]+)$")


def extract() -> dict[str, list[Point]]:
    """Return the N2 file."""
    (cells,) = [row for row in table(DIRECTORY, TABLE)[1:] if row[0] == "N_2"]
    match = _CELL.match(cells[2].strip())
    if match is None:
        msg = f"unexpected cell {cells[2]!r}"
        raise ValueError(msg)
    point = Point(
        "chi_nr",
        PUMP,
        number(match["value"]),
        number(match["sigma"]),
        " | ".join(cells),
    )
    return {"N2/farrow_1987__chi_nr_N2.csv": [point]}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

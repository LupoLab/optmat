# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ward & Elliott (1984): third-harmonic chi of CH4 and SF6 (Table I).

Table I (``tbl_67cee6b847af``, published only as JSON cells; image checked, the cells
agree) gives chi x 10^39 (esu/molecule) for THG at 694 nm (column a, this work), dcSHG
(two columns: footnotes c, e, f at 694 nm, and b) and the Kerr effect (d, g at 633 nm).
Uncertainties are printed in parentheses in units of the last digit, ``330(10)``.
Only CH4 and SF6 are curated; the fluoromethanes are skipped. The per-cell footnotes
decide the source; column b's footnote is in its header.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Ward_1984_OpticalThirdHarmonic"
TABLE = "tbl_67cee6b847af"
GASES = {"CH_4": "CH4", "SF_6": "SF6"}
# (column, footnote) -> (dataset prefix, wavelength in nm as printed in the footnotes)
SOURCES = {
    (1, ""): ("thg", "694"),
    (2, "c"): ("dcshg_miller_ward", "694"),
    (2, "f"): ("dcshg_ward_miller", "694"),
    (3, ""): ("dcshg_shelton_buckingham", "541"),  # printed 541: see the errata
    (4, "d"): ("kerr_buckingham_orr", "633"),
    (4, "g"): ("kerr_buckingham_dunmur", "633"),
}
_CELL = re.compile(r"^(?P<value>[\d.]+)\((?P<sigma>\d+)\)(?:\^(?P<mark>[a-z]))?$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column."""
    files = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES.get(cells[0])
        if gas is None:  # a fluoromethane: not a curated species
            continue
        printed = " | ".join(cells)
        for column in (1, 2, 3, 4):
            match = _CELL.match(cells[column].replace(" ", ""))
            if match is None:
                msg = f"unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            name, wavelength = SOURCES[column, match["mark"] or ""]
            files[f"{gas}/ward_1984__{name}_{gas}.csv"] = [
                Point(
                    "chi", wavelength, number(match["value"]), match["sigma"], printed
                )
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

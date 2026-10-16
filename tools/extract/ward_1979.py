# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ward & Miller (1979): chi_parallel(3) and the comparison data (Table III).

Table III (``tbl_d1cb8fd65ca1``, published only as JSON cells) gives per molecule
R(3) and chi_parallel(3)(-2w; 0, w, w) in 10^-39 esu/molecule, then "other
chi_parallel(3) data": dcSHG (note a), three-wave mixing (TWM, note b), Kerr, THG
(note c) and theory. Only the curated species are extracted: H2, N2, O2, SF6, CO2 and
CO (NO, H2S, NH3, H2O, (CH3)2O and CH3OH are not curated). The cells were compared with
the table image: they agree; each value and its "+-" are printed on two lines in the
image and merged into one cell here, with the footnote letter between them
(``47^d±5``), which is stripped. R(3) (column 1) is a ratio of tensor components, not a
dataset.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Ward_1979_MeasurementsNonlinearOptical"
TABLE = "tbl_d1cb8fd65ca1"
RUBY = "694.3"  # nm, the Q-switched ruby laser (p. 827)
HE_NE = "632.8"  # nm, Buckingham and Dunmur's laser (buckingham_1968, p. 1777)
GASES = {
    "H_2": "H2",
    "N_2": "N2",
    "O_2": "O2",
    "SF_6": "SF6",
    "CO_2": "CO2",
    "CO": "CO",
}
# column -> dataset prefix
COLUMNS = {
    2: "chi3",
    3: "dcshg_mayer",
    4: "twm_rado",
    5: "kerr",
    6: "thg_ward_new",
    7: "theory",
}
# Kerr footnotes -> dataset prefix (d: Buckingham and Orr; f: Buckingham et al. 1970;
# h: Buckingham and Dunmur); theory footnotes (e: McLean and Yoshimine; g: Bartlett
# and Purvis).
KERR = {
    "d": "kerr_buckingham_orr",
    "f": "kerr_buckingham_1970",
    "h": "kerr_buckingham_dunmur",
}
THEORY = {"e": "theory_mclean_yoshimine", "g": "theory_bartlett_purvis"}
_CELL = re.compile(r"^(?P<value>[\d.]+)(?:\^(?P<mark>[a-z]))?(?:±(?P<sigma>[\d.]+))?$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES.get(cells[0])
        if gas is None:  # not a curated species
            continue
        printed = " | ".join(cells)
        for column, prefix in COLUMNS.items():
            if not cells[column]:
                continue
            match = _CELL.match(cells[column].replace(" ", ""))
            if match is None:
                msg = f"unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            name = prefix
            if prefix == "kerr":
                name = KERR[match["mark"]]
            elif prefix == "theory":
                name = THEORY[match["mark"]]
            wavelength = RUBY
            if name == "kerr_buckingham_dunmur":
                wavelength = HE_NE
            elif name.startswith("theory"):
                wavelength = "inf"  # chi_zzzz(0; 0, 0, 0), notes e and g
            files[f"{gas}/ward_1979__{name}_{gas}.csv"] = [
                Point(
                    "chi",
                    wavelength,
                    number(match["value"]),
                    match["sigma"] or "",
                    printed,
                )
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

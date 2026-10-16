# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Lundeen, Hou & Nibler (1983): nonresonant chi(3) of gases (Table I).

Table I (``tbl_4a6abc0d694e``, published only as JSON cells; the cells agree with the
table image) lists, in 10^-39 cm^6/(erg molecule), the calculated rotational and
vibrational Raman contributions to the four-wave-mixing (CARS) chi, then the
electronic FWM chi of this work (532 + 683 nm), of Rado (694 + 975 nm, note d) and of
Hauchecorne et al. (694 + 975 nm, note e; printed with a second "c"), and the dcSHG
(694 nm), Kerr (633 nm) and THG (694 nm) values of others. Structure:

- helium has two sub-rows, "expt." and "theory^i" (Sitz and Yaris);
- values in parentheses are the reference values each FWM column is scaled to (N2's
  105 +- 14 in columns c and d, Ar's 118 in column e) or theoretical helium values;
  the parentheses are stripped and the scale references of other authors' columns are
  not extracted (they are this work's N2, or Rado's Ar);
- "..." is no value; helium's dcSHG "1.3 to 11.4" is a range and is not extracted;
- the atoms' rotational and vibrational contributions are 0 by definition and are
  not extracted.

Only curated species are extracted (NO, CS2, H2O, H2S, NH3 and SiF4 are skipped).
Table II (hydrocarbons and halocarbons) repeats CH4's row and is not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Lundeen_1983_NonresonantThirdOrder"
TABLE = "tbl_4a6abc0d694e"
GASES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "H_2": "H2",
    "D_2": "D2",
    "N_2": "N2",
    "O_2": "O2",
    "CO": "CO",
    "CO_2": "CO2",
    "CH_4": "CH4",
    "SF_6": "SF6",
}
ATOMS = {"He", "Ne", "Ar", "Kr", "Xe"}
# column -> (dataset prefix, pump wavelength in nm as printed in the header)
COLUMNS = {
    1: ("rot", "532"),
    2: ("vib", "532"),
    3: ("fwm_532", "532"),
    4: ("fwm_rado", "694"),
    5: ("fwm_hauchecorne", "694"),
    6: ("dcshg", "694"),
    7: ("kerr", "633"),
    8: ("thg", "694"),
}
# Scale references of other authors' columns: this work's N2 (columns 4) and Rado's
# Ar (column 5), printed in parentheses.
REFERENCES = {(4, "N2"), (5, "Ar")}
_CELL = re.compile(r"^\(?(?P<value>-?[\d.]+)(?:±(?P<sigma>[\d.]+))?\)?$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column."""
    files = {}
    gas = ""
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        label = cells[0]
        if label == "He":  # the helium header row; its values are on the sub-rows
            gas = "He"
            continue
        theory = label.startswith("theory")
        if label not in {"expt."} and not theory:
            gas = GASES.get(label, "")
        if not gas:  # not a curated species
            continue
        printed = " | ".join(cells)
        for column, (prefix, wavelength) in COLUMNS.items():
            cell = cells[column].replace(" ", "")
            if cell in {"", "..."} or " to " in cells[column]:  # none, or a range
                continue
            if prefix in {"rot", "vib"} and gas in ATOMS:
                continue
            if (column, gas) in REFERENCES:
                continue
            match = _CELL.match(cell)
            if match is None:
                msg = f"unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            name = f"{prefix}_theory" if theory else prefix
            files[f"{gas}/lundeen_1983__{name}_{gas}.csv"] = [
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

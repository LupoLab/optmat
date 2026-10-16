# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hohm & Thakkar (2012): literature static polarizabilities of He-Xe.

Table 1 (continued) (``tbl_9ebc9ba30435``, cells file only; image checked) lists, per
atom, the literature polarizabilities alpha_lit (a.u.) with their error estimates and
references: "TL-n" is Table n of Thakkar and Lupinetti's (2006) compilation, "12" is
Gaiser and Fellmuth (2010); mark ``^c`` flags experimental values. A noble-gas atom
has two rows, the second continuing the first with an empty atom cell (the label is
carried down). Only He, Ne, Ar, Kr and Xe are extracted, into one file per gas and kind
(experiment, the ``^c`` values; theory, the others). The model columns alpha_calc
(Eq. 4) and a(Z) (Eq. 5) are not extracted. The "group" cell "18" sits on the At row in
the cells file, one row above He in the image; it is not used.

The same rows give the first and second ionisation potentials I_1 and I_2 (eV, columns
2 and 3; "from ref 3 ... except as indicated", no He-Xe cell is marked): I_1 of the
atom (``ip1_<gas>``) and I_2, the ionisation energy of the singly charged ion
(``ip2_<gas>+``, material He+ ... Xe+), from the first row of each atom.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hohm_2011_NewRelationshipsConnecting"
TABLE = "tbl_9ebc9ba30435"
GASES = {"He": "He", "Ne": "Ne", "Ar": "Ar", "Kr": "Kr", "Xe": "Xe"}
_ENTRY = re.compile(
    r"^(?P<value>[\d.]+) \\pm "
    r"(?P<error>[\d.]+(?: \\times 10\^\{-\d+\})?)(?P<mark>\^c)?$"
)
_TIMES = re.compile(r"^(?P<mantissa>[\d.]+) \\times 10\^\{(?P<exponent>-\d+)\}$")


def _error(text: str) -> str:
    r"""Return a printed error estimate as a number (``1 \times 10^{-6}`` -> 1e-6)."""
    match = _TIMES.match(text)
    return f"{match['mantissa']}e{match['exponent']}" if match else number(text)


def extract() -> dict[str, list[Point]]:
    """Return the polarizabilities per gas and kind, and I_1 and I_2 per gas."""
    files: dict[str, list[Point]] = {}
    atom = ""
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        atom = cells[1].strip().removesuffix("^d") or atom
        if atom not in GASES:
            continue
        match = _ENTRY.match(cells[5].strip())
        if match is None:
            msg = f"unreadable cell {cells[5]!r}"
            raise ValueError(msg)
        if cells[1].strip():  # the atom's first row carries I_1 and I_2
            printed = " | ".join(cells)
            files[f"{atom}/hohm_2012__ip1_{atom}.csv"] = [
                Point("I_1", "", number(cells[2]), "", printed)
            ]
            files[f"{atom}+/hohm_2012__ip2_{atom}+.csv"] = [
                Point("I_2", "", number(cells[3]), "", printed)
            ]
        kind = "experiment" if match["mark"] else "theory"
        gas = GASES[atom]
        path = f"{gas}/hohm_2012__{kind}_{gas}.csv"
        files[path] = [
            Point(
                "static",
                "inf",
                match["value"],
                _error(match["error"]),
                " | ".join(cells),
            )
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

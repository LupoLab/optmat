# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shiner et al. (1993): the ionisation potentials of H2 and D2.

Table IV (``tbl_f2f1c838c202``, cells only; image checked, p. 3): per species the
ionisation potential (cm^-1) with its footnote (a: McCormack et al. and Table III; c:
Jungen et al. 1990, Table III and a private communication; d: HD, not curated; e:
Jungen et al. 1992, which uses Table III), the ab initio value (b: Kolos, private
communication) and the difference. H2 has two rows (a and c); the species cell spans
both. One file per species and footnote; the ab initio column in its own file per
species. HD and the difference column are not extracted (record notes).

The cells print digit-group spaces and the uncertainty in units of the last digit,
"124 417.507(12)^{a}": the spaces are dropped, the parenthesised digits scaled to the
last decimal and the footnote mark kept in ``printed``.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Shiner_1993_H2D2HD"
TABLE_4 = "tbl_f2f1c838c202"
SPECIES = {"H_{2}": "H2", "D_{2}": "D2", "HD": "HD"}
MARKS = {"a": "mccormack", "c": "jungen_1990", "e": "jungen_1992"}
_ENTRY = re.compile(r"^([\d .]+?)(?:\((\d+)\))?\^\{([a-e])\}$")


def _entry(cell: str) -> tuple[str, str, str]:
    """Return the value, the absolute uncertainty ("" if none) and the footnote mark."""
    match = _ENTRY.match(cell.strip())
    if match is None:
        msg = f"unreadable Table IV entry {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1).replace(" ", ""))
    sigma = ""
    if match.group(2) is not None:
        decimals = len(value.partition(".")[2])
        sigma = f"{int(match.group(2)) * 10**-decimals:.{decimals}f}"
    return value, sigma, match.group(3)


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 ionisation potentials of Table IV and the ab initio ones."""
    files: dict[str, list[Point]] = {}
    species = ""
    for cells in json_table(DIRECTORY, TABLE_4)[1:]:
        species = SPECIES.get(cells[0].strip(), species)
        if species == "HD":
            continue  # not a curated species (record notes)
        printed = " | ".join(cells)
        value, sigma, mark = _entry(cells[1])
        path = f"{species}/shiner_1993__{MARKS[mark]}_{species}.csv"
        files[path] = [Point("1", species, value, sigma, printed)]
        if cells[2].strip():
            value, sigma, mark = _entry(cells[2])
            if mark != "b" or sigma:
                msg = f"unexpected ab initio entry {cells[2]!r}"
                raise ValueError(msg)
            path = f"{species}/shiner_1993__ab_initio_{species}.csv"
            files[path] = [Point("1", species, value, "", printed)]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

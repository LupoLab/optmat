# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Manakov, Ovsiannikov & Rapoport (1986): He and Xe rows of Table 5.3.

Table 5.3 (``tbl_4be3eb1d2dd5``, cells file; checked cell by cell against the table
image and the rendered page, PDF page 63 = journal p. 381) gives, in atomic units, the
scalar polarizability alpha^S(omega) and the hyperpolarizability components
gamma_0^l(omega) (linear polarization) and gamma_0^c(omega) (circular) of ground-state
atoms at omega_N = 0.043 a.u. and omega_R = 0.0656 a.u., quoted from the authors'
earlier paper (the review's Ref. [78]). Only the noble gases are taken, He and Xe; the
H and alkali rows have no material and are described in the record's notes.

The two frequencies are read from the merged header cells; the column of each
quantity from the second header row. One file per gas and quantity (and polarization
for gamma), rows ``omega_N`` and ``omega_R``, the x the printed frequency in hartree.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Manakov_1986_AtomsLaserField"
TABLE_5_3 = "tbl_4be3eb1d2dd5"
GASES = ("He", "Xe")
_FREQUENCY = re.compile(r"\\omega_\{\\text\{([NR])\}\}\s*=\s*([\d.]+)")
# The column heads of each frequency block, in print order (second header row).
_HEADS = {
    "alpha": r"\alpha^{\text{S}}(\omega)",
    "gamma_linear": r"\gamma_{0}^{\ell}(\omega)",
    "gamma_circular": r"\gamma_{0}^{\text{c}}(\omega)",
}
HARTREE_EV = 27.211386
# The text gives the photon energies as 1.17 and 1.78 eV (p. 376).
_TEXT_EV = {"N": 1.17, "R": 1.78}


def _frequencies(header: list[str]) -> dict[int, tuple[str, str]]:
    """Return each block's first column: its laser (N or R) and printed frequency."""
    blocks = {}
    for column, cell in enumerate(header):
        match = _FREQUENCY.search(cell)
        if match:
            laser, omega = match.groups()
            if abs(float(omega) * HARTREE_EV - _TEXT_EV[laser]) > 0.01:
                text = _TEXT_EV[laser]
                msg = f"omega_{laser} = {omega} a.u. is not the text's {text} eV"
                raise ValueError(msg)
            blocks[column] = (laser, omega)
    if len(blocks) != 2:
        msg = f"expected two frequency blocks, found {blocks}"
        raise ValueError(msg)
    return blocks


def extract() -> dict[str, list[Point]]:
    """Return alpha^S, gamma_0^l and gamma_0^c of He and Xe at both frequencies."""
    header, heads, *rows = json_table(DIRECTORY, TABLE_5_3)
    blocks = _frequencies(header)
    files: dict[str, list[Point]] = {}
    for gas in GASES:
        (cells,) = (cells for cells in rows if cells[0] == gas)
        for offset, (name, head) in enumerate(_HEADS.items()):
            points = []
            for first, (laser, omega) in blocks.items():
                column = first + offset
                if heads[column] != head:
                    msg = f"column {column} is {heads[column]!r}, not {head!r}"
                    raise ValueError(msg)
                printed = f"omega_{laser} = {omega} a.u. | {head} | {cells[column]}"
                points.append(
                    Point(f"omega_{laser}", omega, number(cells[column]), "", printed)
                )
            files[f"{gas}/manakov_1986__{name}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

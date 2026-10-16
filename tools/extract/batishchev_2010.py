# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Batishchev, Tolstikhin & Morishita (2010): the static-field width of Ar(3p).

Table I (``tbl_c150ef90d0bb``, JSON cells only; image checked, p. 023416-8) shows the
convergence of the Siegert eigenvalue E = script-E - i Gamma/2 of Ar(3p, m = 0) in the
single-active-electron potential with the number of adiabatic channels N_ch, for F = 0,
0.05, 0.1, 0.2, 0.5 and 1.0 au: columns F, N_ch, -script-E, Gamma, with F printed only
on the first row of each block. "All digits in the last line for each value of F are
converged" (caption), so the last row of each block is the result; the others are the
convergence study (record notes). The width Gamma of each F > 0 block's last row is
written; F = 0 (Gamma = 0, the bound state) is not a rate and is left in the notes.

Structural fixes (the image shows four columns):
- The cells merge the header row with the first data row ("F Nch0 5", "-E0.600804",
  "0") and, in each later block, F with N_ch ("0.05 5"); they are split on the space.
- Gamma prints "a(b)" for a x 10^b ("0.270883(-5)", caption); it is written
  "0.270883e-5", ``printed`` keeping the cell. The digit-group spaces of the image
  ("0.270 883") are absent from the cells.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Batishchev_2010_AtomicSiegertStates"
TABLE_I = "tbl_c150ef90d0bb"
_POWER = re.compile(r"^([\d.]+)\((-\d+)\)$")
_FIRST = re.compile(r"^F Nch(\S+) (\d+)$")


def width(cell: str) -> str:
    """Return a printed width 'a(b)' as 'aeb', or a plain number as printed."""
    if match := _POWER.match(cell.strip()):
        return f"{match.group(1)}e{match.group(2)}"
    return number(cell)


def blocks() -> dict[str, list[tuple[str, str, str]]]:
    """Return each field's rows (N_ch, -script-E, Gamma) in printed order."""
    grid = json_table(DIRECTORY, TABLE_I)
    found: dict[str, list[tuple[str, str, str]]] = {}
    field = ""
    for index, (label, printed_energy, gamma) in enumerate(grid):
        first = label.strip()
        energy = printed_energy
        if not first:
            continue  # the empty second row of the merged header
        if index == 0:
            match = _FIRST.match(first)
            if match is None or not energy.startswith("-E"):
                msg = f"unexpected Table I first row {grid[0]!r}"
                raise ValueError(msg)
            field, channels, energy = match.group(1), match.group(2), energy[2:]
        elif " " in first:
            field, channels = first.split(" ")
        else:
            channels = first
        found.setdefault(field, []).append((channels, energy.strip(), gamma.strip()))
    if list(found) != ["0", "0.05", "0.1", "0.2", "0.5", "1.0"]:
        msg = f"unexpected Table I fields {list(found)}"
        raise ValueError(msg)
    return found


def extract() -> dict[str, list[Point]]:
    """Return Ar(3p)'s converged widths against field (F > 0)."""
    points = []
    for field, rows in blocks().items():
        if field == "0":
            continue
        channels, energy, gamma = rows[-1]
        printed = f"F = {field} | N_ch = {channels} | -E = {energy} | Gamma = {gamma}"
        points.append(Point(field, number(field), width(gamma), "", printed))
    return {"Ar/batishchev_2010__table1_Ar.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

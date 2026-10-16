# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Santra, Dunford & Young (2006): Kr 4p channel rates at F = 0.10 au.

Table I (``tbl_f7044803390d``, CSV; image checked, p. 043403-3) prints, for the three
metastable 4p states of krypton at the static field E = 0.10 au (the caption), the
projection m, the Stark shift, the tunnel ionisation rate gamma (au) and the squared
overlaps |c(j=1/2)|^2 and |c(j=3/2)|^2 with the field-free 4p_j states. The rows carry
no j; it is read from the overlaps (0.976 for j = 1/2; 0.957 and 0.970 for j = 3/2),
as the text names them (gamma_{3/2,+-1/2} > 8 x gamma_{3/2,+-3/2} and 2.4 x
gamma_{1/2,+-1/2}, p. 043403-4); the script checks that each row's larger overlap is
the j it is given.

One file per channel, ``tableI_j<2j>_m<2|m|>_Kr``, one row: the row id and ``x`` are
the field of the caption, the value is gamma. The cells print per-cell powers of ten
(``4.89 \times 10^{-3}``), turned here into ``4.89e-3``; ``printed`` keeps the row.
"""

import re

from _library import Point, table, write

DIRECTORY = "Santra_2006_SpinOrbitEffect"
TABLE_I = "tbl_f7044803390d"
FIELD = "0.10"  # the caption: "at an electric field strength of E=0.10 a.u."
HEADER = [
    "m",
    r"\Delta E^{(Stark)} (a.u.)",
    r"\gamma (a.u.)",
    "|c(j=1/2)|^2",
    "|c(j=3/2)|^2",
]
# Row order of Table I -> (j, |m|) as file-name tags, and the m the row prints.
CHANNELS = (("12", "12", "+1/2"), ("32", "12", "+1/2"), ("32", "32", "+3/2"))
_POWER = re.compile(r"^\s*([-+]?\d+\.\d+)\s*\\times\s*10\^\{([-+]?\d+)\}\s*$")


def _number(cell: str) -> str:
    r"""Return ``a \times 10^{k}`` as ``ae<k>``, the printed digits kept."""
    match = _POWER.match(cell)
    if match is None:
        msg = f"not a power-of-ten cell: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return the rate of each 4p channel at F = 0.10 au."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER or len(grid) != 1 + len(CHANNELS):
        msg = f"unexpected Table I: {grid!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells, (j, m, printed_m) in zip(grid[1:], CHANNELS, strict=True):
        mostly_half = float(cells[3]) > float(cells[4])  # the j = 1/2 overlap leads
        if cells[0] != printed_m or (j == "12") != mostly_half:
            msg = f"Table I row {cells!r} is not the channel j = {j}/2, m = {m}/2"
            raise ValueError(msg)
        point = Point(FIELD, FIELD, _number(cells[2]), "", " | ".join(cells))
        files[f"Kr/santra_2006__tableI_j{j}_m{m}_Kr.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Vatter, Schmidt, Elias & Schramm (1996): interaction second virials.

Table 1 (``tbl_9d39962ae383``, p. 75, cells file only; image checked against the cells
and a 400 dpi rendering of the page, 2026-10-01) prints B in cm^3/mol at 296.15,
353.15, 413.15 and 463.15 K: first the pure-gas values the authors took from the
literature (Ne, Ar, O2, N2; not extracted, they are other works' values and are listed
in the record's notes), then "This work": pure CO and the interaction (unlike-pair)
coefficients B_12 of Ne/O2, Ar/O2, Ne/N2, Ar/N2, Ne/CO and Ar/CO, each with its
estimated error (``-1.4±1.0``).

Each pair is filed under the first of its two gases in sorted order, with the other as
the dataset's partner. Two cells carry a table note and are not this work's
measurement: CO at 296.15 K (a, "Taken from [21]", Michels et al. 1952) and Ar/N2 at
296.15 K (b, "Taken from [22]", Ewing & Trusler 1992). They go to their own files, as
reproduced values. Tables 2-4 are not extracted (the record's notes say why).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Vatter_1996_MeasurementsInteractionSecond"
TABLE_1 = "tbl_9d39962ae383"
SOURCE = "vatter_1996"
TEMPERATURES = ("296.15 K", "353.15 K", "413.15 K", "463.15 K")
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
# The rows of "This work:" (the row label as printed: the gas, or the pair).
THIS_WORK = ("CO", "Ne/O2", "Ar/O2", "Ne/N2", "Ar/N2", "Ne/CO", "Ar/CO")
# The table-note cells: (row label, temperature) -> (note mark, file suffix).
NOTED = {
    ("CO", "296.15 K"): ("a)", "michels"),
    ("Ar/N2", "296.15 K"): ("b)", "ewing"),
}


def split(cell: str) -> tuple[str, str, str]:
    """Return the value, the error and the table-note mark of ``-11.3±0.5b)``."""
    value, _, rest = cell.partition(PLUS_MINUS)
    mark = ""
    if rest.endswith(("a)", "b)")):
        rest, mark = rest[:-2], rest[-2:]
    return number(value), number(rest), mark


def stem(label: str) -> tuple[str, str]:
    """Return the material of a row and the dataset id's stem (``B12_N2_Ne``)."""
    if "/" not in label:
        return label, f"table1_B_{label}"
    first, second = sorted(label.split("/"))
    return first, f"table1_B12_{first}_{second}"


def extract() -> dict[str, list[Point]]:
    """Return Table 1's results of this work, one file per gas or pair."""
    grid = json_table(DIRECTORY, TABLE_1)
    if tuple(grid[1][1:]) != TEMPERATURES or grid[7][0] != "This work:":
        msg = f"unexpected Table 1 layout: {grid[1]!r}, {grid[7]!r}"
        raise ValueError(msg)
    if tuple(cells[0] for cells in grid[8:]) != THIS_WORK:
        msg = f"unexpected Table 1 rows {[cells[0] for cells in grid[8:]]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    noted = set()
    for cells in grid[8:]:
        label = cells[0]
        material, name = stem(label)
        for heading, cell in zip(TEMPERATURES, cells[1:], strict=True):
            if not cell:
                continue
            value, sigma, mark = split(cell)
            temperature = heading.removesuffix(" K")
            dataset = name
            if mark:
                expected, suffix = NOTED.get((label, heading), ("", ""))
                if mark != expected:
                    msg = f"Table 1 {label} at {heading}: unexpected note {mark!r}"
                    raise ValueError(msg)
                noted.add((label, heading))
                dataset = f"{name}_{suffix}"
            printed = f"{label} | {heading} | {cell}"
            files.setdefault(f"{material}/{SOURCE}__{dataset}.csv", []).append(
                Point(temperature, temperature, value, sigma, printed)
            )
    if noted != set(NOTED):
        msg = f"Table 1 notes found at {sorted(noted)!r}, expected {sorted(NOTED)!r}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

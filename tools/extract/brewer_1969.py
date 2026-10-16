# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Brewer & Vaughn (1969): excess and interaction second virial coefficients.

Table IV (``tbl_5752880ca58e``, p. 2966) prints, for ten binaries of He, Ne, Ar, H2 and
N2 at -125 to 50 degC in steps of 25 degC, the measured excess second virial
coefficient E and the interaction (cross) second virial coefficient B12, both in
cc/mole, in two column blocks of System, Temp (degC), E, B12. Each block is read row by
row. A pair is filed under the first of its two species in sorted order, with the
other as the partner, whichever way round the table names it (``N2-Ar`` and ``Ar-N2``
are the same system): ``<material>/brewer_1969__E_<material>_<partner>.csv`` and
``..._B12_...``, each sorted by temperature.

One label of the table is misprinted, checked against the table image (2026-10-01):
the 25 degC row ``H2-Ne`` (E = 2.83, B12 = 16.01) is the H2-He system. Ne-H2 has its own
25 degC row (1.14, 14.23), H2-He has no other 25 degC row, its E lies on the H2-He run
(2.38 at 50 degC, 3.50 at 0 degC), and its B12 is E + (B_H2 + B_He)/2 from Table I
(16.015), not E + (B_H2 + B_Ne)/2 (15.92). The ``printed`` cell keeps the label as
printed.

Every row's B12 is checked here against E + (B11 + B22)/2 with the pure-gas values of
Table I (``tbl_8e22570edb7a``, p. 2963), the paper's Eq. (12): all 80 agree within half
a unit of the last printed digit.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Brewer_1969_MeasurementCorrelationSome"
TABLE_I = "tbl_8e22570edb7a"
TABLE_IV = "tbl_5752880ca58e"
HEADER = ["System", "Temp (°C)", "E (cc/mole)", "B12 (cc/mole)"]
LABEL_FIXES = {
    # (row index in the extracted table, column): (extracted, printed on the page)
    (18, 0): ("H2-Ne", "H2-He"),
}
GASES = {"He": "He", "H_2": "H2", "Ar": "Ar", "N_2": "N2", "Ne": "Ne"}
TEMPERATURES = 8
TOLERANCE = 0.0051  # half a unit of the printed hundredths, with float slack


def pure_coefficients() -> dict[tuple[str, str], float]:
    """Return Table I's literature B(T) of each gas, keyed by (temperature, gas)."""
    grid = json_table(DIRECTORY, TABLE_I)
    gases = [GASES[cell] for cell in grid[0][1:]]
    return {
        (number(cells[0]), gas): float(number(cell))
        for cells in grid[1:]
        for gas, cell in zip(gases, cells[1:], strict=True)
    }


def extract() -> dict[str, list[Point]]:
    """Return Table IV's E and B12, one file of each per pair."""
    grid = json_table(DIRECTORY, TABLE_IV)
    if grid[0] != HEADER * 2:
        msg = f"unexpected Table IV header {grid[0]!r}"
        raise ValueError(msg)
    labels = [row[:] for row in grid]
    for (row, column), (extracted, printed) in LABEL_FIXES.items():
        if labels[row][column] != extracted:
            msg = f"Table IV cell {(row, column)} is {labels[row][column]!r}"
            raise ValueError(msg)
        labels[row][column] = printed
    pure = pure_coefficients()
    files: dict[str, list[Point]] = {}
    for index, cells in enumerate(grid[1:], start=1):
        for start in (0, 4):
            block = cells[start : start + 4]
            pair = sorted(labels[index][start].split("-"))
            temperature, excess, cross = (number(cell) for cell in block[1:])
            mean = (pure[temperature, pair[0]] + pure[temperature, pair[1]]) / 2
            if abs(float(cross) - float(excess) - mean) > TOLERANCE:
                msg = f"Table IV {block!r}: B12 is not E + (B11 + B22)/2"
                raise ValueError(msg)
            printed = " | ".join(block)
            for name, value in (("E", excess), ("B12", cross)):
                key = f"{pair[0]}/brewer_1969__{name}_{pair[0]}_{pair[1]}.csv"
                files.setdefault(key, []).append(
                    Point(temperature, temperature, value, "", printed)
                )
    for key, points in files.items():
        if len({p.row for p in points}) != TEMPERATURES:
            msg = f"{key}: {len(points)} rows, not one at each of {TEMPERATURES}"
            raise ValueError(msg)
        points.sort(key=lambda p: float(p.x))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

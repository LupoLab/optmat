# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Martin, Trengove, Harris & Dunlop (1982): excess second virial coefficients.

The numbers print a raised decimal point (``20·4``), kept in the ``printed`` column and
written as "." in the numbers. The library attaches the captions of Tables 2 and 3 to
the wrong tables (Table 2's caption to Table 1's notes, the label "Table 3" to Table 2's
grid); the page image (p. 1526) shows the layout read here. All three grids agree with
their images cell for cell (checked 2026-10-01).

- Table 1 (``tbl_8885f566057e``, CSV): the pure gases' second virial coefficients
  B/(cm^3 mol^-1) at 290, 300 and 320 K the authors took from Dymond and Smith's
  compilation, in two column groups. One file per gas (row = the temperature in K).
- Table 2 (``tbl_11bc84c101b5``, cells file only): the individual runs for Ar-N2, two
  per temperature (the temperature printed on the first of each pair): P_1/kPa,
  Delta P/kPa (both with a subscript guard digit, ``190·043_{7}``, kept only in
  ``printed``), E and B_12. Rows ``<T>_1`` and ``<T>_2``.
- Table 3 (``tbl_856ba5037a07``, cells file only): E and B_12 for twenty pairs at 290,
  300 and 320 K. A row label ``+X`` continues the group of the last named first species
  (``Ar+He``, then ``+H_2`` for Ar+H2). One file per pair and quantity, filed under the
  first species in sorted order with the other as the partner. Pairs with CF4 are not
  optmat gases and are skipped (the record's notes).
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Martin_1982_ExcessSecondVirial"
TABLE_1 = "tbl_8885f566057e"
TABLE_2 = "tbl_11bc84c101b5"
TABLE_3 = "tbl_856ba5037a07"
DOT = "\N{MIDDLE DOT}"
TEMPERATURES = ("290", "300", "320")
GASES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "H_{2}": "H2",
    "H_2": "H2",
    "N_{2}": "N2",
    "N_2": "N2",
    "O_{2}": "O2",
    "O_2": "O2",
    "CH_{4}": "CH4",
    "CH_4": "CH4",
    "CO_{2}": "CO2",
    "CO_2": "CO2",
    "CF_{4}": "CF4",
    "CF_4": "CF4",
}
SKIPPED = frozenset({"CF4"})  # not an optmat gas


def _value(cell: str) -> str:
    """Return the number of a cell printed with a raised decimal point."""
    return number(cell.replace(DOT, "."))


def _pure() -> dict[str, list[Point]]:
    """Return Table 1: the pure gases' B from Dymond and Smith, one file per gas."""
    grid = table(DIRECTORY, TABLE_1)
    if grid[0] != ["Gas", "290 K", "300 K", "320 K"] * 2:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        for start in (0, 4):
            if not cells[start]:
                continue
            gas = GASES[cells[start]]
            if gas in SKIPPED:
                continue
            group = cells[start : start + 4]
            printed = " | ".join(group)
            files[f"{gas}/martin_1982__B_{gas}.csv"] = [
                Point(t, t, _value(cell), "", printed)
                for t, cell in zip(TEMPERATURES, group[1:], strict=True)
            ]
    return files


def _runs() -> dict[str, list[Point]]:
    """Return Table 2: the Ar-N2 runs' E and B_12."""
    grid = json_table(DIRECTORY, TABLE_2)
    if not (grid[0][0] == "T/K" and "E" in grid[0][3] and "B_{12}" in grid[0][4]):
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    excess: list[Point] = []
    cross: list[Point] = []
    temperature, run = "", 0
    for cells in grid[1:]:
        if cells[0]:
            temperature, run = cells[0], 0
        run += 1
        row, printed = f"{temperature}_{run}", " | ".join(cells)
        excess.append(Point(row, temperature, _value(cells[3]), "", printed))
        cross.append(Point(row, temperature, _value(cells[4]), "", printed))
    return {
        "Ar/martin_1982__E_runs_Ar_N2.csv": excess,
        "Ar/martin_1982__B12_runs_Ar_N2.csv": cross,
    }


def _pairs() -> dict[str, list[Point]]:
    """Return Table 3: E and B_12 per pair, one file per pair and quantity."""
    grid = json_table(DIRECTORY, TABLE_3)
    if grid[0][1:6:2] != [f"{t} K" for t in TEMPERATURES] or grid[1][1:3] != [
        r"\mathcal{E}",
        "B_{12}",
    ]:
        msg = f"unexpected Table 3 header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    first = ""
    for cells in grid[2:]:
        named, _, second = cells[0].partition("+")
        first = GASES[named] if named else first
        partner = GASES[second]
        if {first, partner} & SKIPPED:
            continue
        material, partner = sorted((first, partner))
        printed = " | ".join(cells)
        for name, offset in (("E", 1), ("B12", 2)):
            files[f"{material}/martin_1982__{name}_{material}_{partner}.csv"] = [
                Point(t, t, _value(cells[offset + 2 * i]), "", printed)
                for i, t in enumerate(TEMPERATURES)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1, 2 and 3."""
    return _pure() | _runs() | _pairs()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

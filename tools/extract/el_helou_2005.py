# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract El Helou et al. (2005): O3 cross sections in the Chappuis and Wulf bands.

All tables are published as cells only (``.json``), each cell "value (uncertainty)" in
10^-23 cm^2 molecule^-1 (images checked, pp. 4-9). Abscissae are FTS wavenumbers in
cm^-1 (vacuum): the five He-Ne lines of the Chappuis band, and eight Wulf band maxima.

* **Table IV** (``tbl_69e3ef0c7fa6``, p. 6) gathers every Chappuis-band row: Anderson
  (reproduced), Orphal (reproduced), the LASIM cell at 293, 223 and 144 K and the Reims
  cell (FTS) at 291 and 225 K. It is the source of those files, because it labels the
  Reims rows' temperatures correctly.
* **Table I** (``tbl_f0168054b0c5``), **Table III** (``tbl_482199636ba1``) and **Table
  II**'s FTS row (``tbl_2cd50501f353``) repeat rows of Table IV: checked equal here, not
  written. Table III labels both its rows "T=225 K"; its first row is the 291 K row of
  Table IV and of its own d sigma / sigma line "(291-225 K)" (image checked).
* **Table II**'s first row (Reims cell + monochromator, 295 K) is Brion et al.'s
  measurement (brion_1998): written as a reproduced file.
* **Table VI** (``tbl_b54e332b1d7d``): the 175 K row (its 293 K row is Table I's,
  checked).
* **Table VII** (``tbl_f1e1ad15b4c2``, p. 8): the Wulf bands, LASIM cell at 294, 222 and
  150 K, and Anderson et al. (1993) at 293 K (reproduced; four bands). **Table VIII**
  (``tbl_8a209ae7d01c``) repeats the 294 and 150 K values: checked equal (its last row
  prints the uncertainties 0.35 where Table VII prints 0.36; the values agree).

Register: TC-972.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Helou_2005_AbsoluteCrossSections"
TABLE_I = "tbl_f0168054b0c5"
TABLE_II = "tbl_2cd50501f353"
TABLE_III = "tbl_482199636ba1"
TABLE_IV = "tbl_69e3ef0c7fa6"
TABLE_VI = "tbl_b54e332b1d7d"
TABLE_VII = "tbl_f1e1ad15b4c2"
TABLE_VIII = "tbl_8a209ae7d01c"

_CELL = re.compile(r"^([\d.]+) \(([\d.]+)\)$")
_WAVENUMBER = re.compile(r"^(\d+\.\d+)")

# Table IV's row headers (as extracted) and the dataset each row is written to.
CHAPPUIS = {
    "Anderson^a \\ \\sigma (\\pm 1 SD)": "anderson_1992",
    "Orphal^b": "orphal_2003",
    "Our cell^c (293 K)": "lasim_293",
    "Reims cell (291 K)": "reims_291",
    "Our cell^c (223 K)": "lasim_223",
    "Reims cell (225 K)": "reims_225",
    "Our cell^c (144 K)": "lasim_144",
}
# Table VII's columns, in order after the wavenumber.
WULF = ("anderson_1993_wulf", "lasim_wulf_294", "lasim_wulf_222", "lasim_wulf_150")


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the bracketed uncertainty of a cell, as printed."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def _rows(table_id: str, first: int) -> dict[str, list[str]]:
    """Return a He-Ne table's rows below its header, keyed by the row header."""
    return {cells[0]: cells[1:] for cells in json_table(DIRECTORY, table_id)[first:]}


def _same(label: str, cells: list[str], reference: list[str]) -> None:
    """Fail unless a repeated row prints the same values as the row it repeats."""
    if [_split(c)[0] for c in cells] != [_split(c)[0] for c in reference]:
        msg = f"{label}: {cells} differs from {reference}"
        raise ValueError(msg)


def _chappuis(lines: list[str], cells: list[str], label: str) -> list[Point]:
    """Return one Chappuis-band row as points at the five He-Ne wavenumbers."""
    printed = " | ".join([label, *cells])
    points = []
    for line, cell in zip(lines, cells, strict=True):
        value, sigma = _split(cell)
        points.append(Point(line, line, value, sigma, printed))
    return points


def extract() -> dict[str, list[Point]]:
    """Return the Chappuis- and Wulf-band cross sections as data files."""
    grid = json_table(DIRECTORY, TABLE_IV)
    lines = grid[1][1:]
    table_iv = {cells[0]: cells[1:] for cells in grid[2:]}
    files: dict[str, list[Point]] = {}
    for label, dataset in CHAPPUIS.items():
        files[f"O3/el_helou_2005__{dataset}.csv"] = _chappuis(
            lines, table_iv[label], label
        )

    table_i = _rows(TABLE_I, 2)
    _same("Table I 293 K", table_i["T=293 K"], table_iv["Our cell^c (293 K)"])
    _same("Table I 223 K", table_i["T=223 K"], table_iv["Our cell^c (223 K)"])
    _same("Table I 144 K", table_i["T=144 K"], table_iv["Our cell^c (144 K)"])
    # Table III prints "T=225 K" on both rows; the first is the 291 K row.
    table_iii = [cells[1:] for cells in json_table(DIRECTORY, TABLE_III)[2:4]]
    _same("Table III first row", table_iii[0], table_iv["Reims cell (291 K)"])
    _same("Table III second row", table_iii[1], table_iv["Reims cell (225 K)"])
    table_ii = json_table(DIRECTORY, TABLE_II)[2:4]
    _same("Table II FTS row", table_ii[1][1:], table_iv["Reims cell (291 K)"])
    files["O3/el_helou_2005__brion_1998_mono.csv"] = _chappuis(
        lines, table_ii[0][1:], table_ii[0][0]
    )
    table_vi = _rows(TABLE_VI, 2)
    _same("Table VI 293 K", table_vi["T=293 K"], table_iv["Our cell^c (293 K)"])
    files["O3/el_helou_2005__lasim_175.csv"] = _chappuis(
        lines, table_vi["T=175 K"], "T=175 K"
    )

    wulf: dict[str, list[str]] = {}
    for cells in json_table(DIRECTORY, TABLE_VII)[2:]:
        match = _WAVENUMBER.match(cells[0])
        if match is None:
            msg = f"unreadable Wulf row {cells[0]!r}"
            raise ValueError(msg)
        line = match.group(1)
        wulf[line] = [cells[2], cells[4]]  # 294 K and 150 K, as Table VIII
        printed = " | ".join(cells)
        for dataset, cell in zip(WULF, cells[1:], strict=True):
            if not cell:
                continue  # Anderson et al. measured four of the eight bands
            value, sigma = _split(cell)
            files.setdefault(f"O3/el_helou_2005__{dataset}.csv", []).append(
                Point(line, line, value, sigma, printed)
            )
    for cells in json_table(DIRECTORY, TABLE_VIII)[2:]:
        line = cells[0].split()[0]
        _same(f"Table VIII {line}", cells[1:3], wulf[line])
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

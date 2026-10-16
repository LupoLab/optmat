# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Komasa et al. (2011), J. Chem. Theory Comput. 7, 3105: H2 and D2 levels.

Only the tables the vocabulary can hold are taken (the record's notes say why Tables
1-3 and 5 are not): the vibrational excitation energies (Tables 4 and 6) and the
quadrupole (Raman-allowed) transition energies (Tables 7 and 8). Every cell prints a
value with its uncertainty on the last digits, "4161.1661(9)"; the sigma column is that
uncertainty in the value's units (``_concise``). Thin spaces inside numbers are dropped.
All tables were checked against their images.

* Table 4 (``tbl_e9a74a8c3fc7``, JSON, p. 3111): the v = 0 -> 1 excitation of H2 and
  D2, theory and experiment (the J = 0 -> 1 column is notes: not a Raman transition).
* Table 6 (``tbl_9ca11361dde2``, CSV, p. 3112): H2 E(v, J = 0) - E(0, 0), v = 1-14,
  theory (the differences from Dabrowski and Wolniewicz stay in ``printed``).
* Table 7 (``tbl_e291b29bef99``, JSON, p. 3113): H2 S0, Q1, S1, Q2, S2 lines, theory
  and experiment. **Structural fix:** the library collapsed the S0(J) block (J = 0-5)
  into one row whose cells hold whole columns run together, and split the column of
  Ref. 53 across two cells ("experim354.373587.032..." and "ent5350(40)11(17)...": the
  first three decimals of each number, then Ref. 53's superscript and each number's
  last two digits with its uncertainty). The script reads the numbers from those
  strings in order and checks their count; the image shows the six rows.
* Table 8 (``tbl_7ba5682f91e7``, JSON, p. 3113): D2 S0, S1, O1, Q1 lines, theory and
  experiment.

Rows are the transitions ("S0(J)", "Q1(J)", ...) for Tables 7-8 and v for Tables 4, 6.
One file per molecule and origin (this work, or each experiment's reference).
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Komasa_2011_QuantumElectrodynamicsEffects"
TABLE_4 = "tbl_e9a74a8c3fc7"
TABLE_6 = "tbl_9ca11361dde2"
TABLE_7 = "tbl_e291b29bef99"
TABLE_8 = "tbl_7ba5682f91e7"
_CONCISE = re.compile(r"^(?P<value>\d+\.\d+)\((?P<digits>\d+)\)")
_NOTE = re.compile(r"\^\{?[a-d]\}?$")


def _concise(cell: str) -> tuple[str, str]:
    """Split "4161.1661(9)" into ("4161.1661", "0.0009"): sigma on the last digits."""
    text = _NOTE.sub("", cell.replace(" ", "").replace("\N{THIN SPACE}", ""))
    match = _CONCISE.fullmatch(text) or _CONCISE.match(text)
    if match is None or match.end() != len(text):
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["value"], match["digits"]
    decimals = len(value.partition(".")[2])
    return value, f"{int(digits) / 10**decimals:.{decimals}f}"


def _all_concise(text: str) -> list[tuple[str, str]]:
    """Return every "v(d)" of a run-together string, in order."""
    found = re.findall(r"(\d+\.\d+)\((\d+)\)", text.replace(" ", ""))
    return [_concise(f"{value}({digits})") for value, digits in found]


def _add(
    files: dict[str, list[Point]], name: str, row: str, cell: str, printed: str
) -> None:
    if not cell:
        return
    value, sigma = _concise(cell)
    files.setdefault(name, []).append(Point(row, row, value, sigma, printed))


def _table_4() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_4)
    if grid[1][2] != "\\nu = 0 \\rightarrow 1" or grid[5][1] != "D_{2}":
        msg = f"unexpected Table 4 layout {grid!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for gas, first in (("H2", 2), ("D2", 7)):
        for label, stem in (("theory", "table4_theory"), ("experiment", "table4_expt")):
            cells = next(row for row in grid[first : first + 3] if row[0] == label)
            _add(
                files,
                f"{gas}/komasa_2011__{stem}_{gas}.csv",
                "1",
                cells[2],
                " | ".join(cells),
            )
    return files


def _table_6() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_6)
    if grid[0][:2] != ["\\nu", "theory"]:
        msg = f"unexpected Table 6 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        _add(
            files,
            "H2/komasa_2011__table6_theory_H2.csv",
            cells[0],
            cells[1],
            " | ".join(cells),
        )
    return files


def _s0_block(block: list[str]) -> dict[str, list[Point]]:
    """Read Table 7's collapsed S0(J) row (see the docstring)."""
    rows = [f"S0({j})" for j in range(6)]
    if block[0] != "J012345" or not block[1].startswith("theory"):
        msg = f"unexpected S0 block {block!r}"
        raise ValueError(msg)
    theory = _all_concise(block[1].removeprefix("theory"))
    # Ref. 53: the leading digits in cell 2 (three decimals each), the rest in cell 3.
    heads = re.findall(
        r"\d[\d ]*\.\d{3}", re.sub(r"^[a-z]+", "", block[2]).removesuffix("Q1(J)")
    )
    tails = re.findall(r"(\d+)\((\d+)\)", block[3].removeprefix("ent53"))
    ref53 = [
        _concise(f"{h.replace(' ', '')}{t}({d})")
        for h, (t, d) in zip(heads, tails, strict=True)
    ]
    ref29 = _all_concise(block[5].removeprefix("experiment29").removesuffix("S1(J)"))
    if not len(theory) == len(ref53) == len(ref29) == len(rows):
        msg = f"S0 block: {len(theory)}, {len(ref53)}, {len(ref29)} values, not 6"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for row, (tv, ts), (av, asg), (bv, bs) in zip(
        rows, theory, ref53, ref29, strict=True
    ):
        printed = f"{row} | {tv}({ts}) | {av}({asg}) | {bv}({bs})"
        for name, value, sigma in (
            ("table7_theory", tv, ts),
            ("table7_jennings_brault", av, asg),
            ("table7_salumbides", bv, bs),
        ):
            files.setdefault(f"H2/komasa_2011__{name}_H2.csv", []).append(
                Point(row, row, value, sigma, printed)
            )
    return files


def _pairs(
    grid: list[list[str]],
    rows: slice,
    branches: tuple[str, str],
    stem: str,
    experiments: tuple[str, str],
) -> dict[str, list[Point]]:
    """Read a block of rows J | theory | expt | delta | theory | expt | delta.

    ``stem`` is the table ("table7"); ``experiments`` names each half's reference.
    """
    gas = "H2" if stem == "table7" else "D2"
    files: dict[str, list[Point]] = {}
    for cells in grid[rows]:
        printed = " | ".join(cells)
        for (theory, expt), branch, source in zip(
            ((1, 2), (4, 5)), branches, experiments, strict=True
        ):
            row = f"{branch}({cells[0]})"
            prefix = f"{gas}/komasa_2011__{stem}"
            _add(files, f"{prefix}_theory_{gas}.csv", row, cells[theory], printed)
            _add(files, f"{prefix}_{source}_{gas}.csv", row, cells[expt], printed)
    return files


def _merge(*parts: dict[str, list[Point]]) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for part in parts:
        for name, points in part.items():
            files.setdefault(name, []).extend(points)
    return files


def _table_7() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_7)
    if grid[2] != ["J", "theory", "experiment54", "δ", "theory", "experiment54", "δ"]:
        msg = f"unexpected Table 7 layout {grid[2]!r}"
        raise ValueError(msg)
    bragg = ("bragg", "bragg")  # Ref. 54 for Q1, S1, Q2 and S2
    return _merge(
        _s0_block(grid[1]),
        _pairs(grid, slice(3, 8), ("Q1", "S1"), "table7", bragg),
        _pairs(grid, slice(9, 13), ("Q2", "S2"), "table7", bragg),
    )


def _table_8() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_8)
    if grid[0][1] != "S_0(J)" or grid[9][1] != "O_1(J)":
        msg = f"unexpected Table 8 layout {grid[0]!r}, {grid[9]!r}"
        raise ValueError(msg)
    mckellar = "mckellar_oka"  # Ref. 56; Ref. 55 (Jennings et al.) for S0
    return _merge(
        _pairs(grid, slice(2, 9), ("S0", "S1"), "table8", ("jennings", mckellar)),
        _pairs(grid, slice(11, 15), ("O1", "Q1"), "table8", (mckellar, mckellar)),
    )


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Tables 4 and 6-8."""
    return _table_4() | _table_6() | _merge(_table_7(), _table_8())


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

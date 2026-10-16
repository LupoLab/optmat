# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Weber (1979): rotational constants from pure rotational Raman spectra.

Chapter 3 (pp. 71-121) of *Raman Spectroscopy of Gases and Liquids*, A. Weber, ed.,
Topics in Current Physics 11 (Springer, 1979); the library's copy is the whole book
(``Weber_1979_RamanSpectroscopyGases``, PDF page n is book page n - 11 here). Both
tables were checked against their images.

* Table 3.2 (``tbl_fefbb3d9eb51``, p. 88; the JSON cells match the image): B_0 and D_0
  of diatomic molecules, each row with its reference. Taken: H2 and D2 (Stoicheff,
  pre-laser, in square brackets, no uncertainty), 14N2 (Butcher et al.'s two- and
  three-constant fits, Bendtsen) and 16O2 (Butcher et al., Edwards et al.). One file per
  molecule, constant and reference; Butcher et al.'s two 14N2 fits are the two rows of
  their files. The isotopologues (HD, HT, DT, T2, 14N15N, 15N2, 16O17O, 16O18O, 17O2,
  18O2) and the other molecules are in the record's notes.
* Table 3.3 (``tbl_c0efcef3f5a7``, p. 93): CO2's B_0 and D_0 in the ground state
  (00^0 0) and for the c levels of (01^1 0).

STRUCTURAL FIX (Table 3.3). The library's cells join the two CO2 rows (the ground and
the excited state) into one ('0.3902634±0.00000440.3906398±0.0000082'); the two lines
are transcribed from the image below and checked to be exactly the library cell when
joined.

The cells are parsed as printed: 'value±sigma' or '(value ±sigma)x10^{-k}' (with a
multiplication sign; the power is the dataset's scale), '[value]' for Stoicheff's
pre-laser values (no uncertainty), and a trailing note letter (^a, ' a') dropped from
the number and kept in the printed column.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Weber_1979_RamanSpectroscopyGases"
TABLE_3_2 = "tbl_fefbb3d9eb51"
TABLE_3_3 = "tbl_c0efcef3f5a7"
STEM = "weber_1979__"
TIMES = "\N{MULTIPLICATION SIGN}"
_CELL = re.compile(
    r"^\s*(?:\[(?P<bracket>[\d.]+)\]"
    r"|\(?(?P<value>[\d.]+)\s*±\s*(?P<sigma>[\d.]+)\)?\s*(?:\u00d710\^\{-(?P<power>\d+)\})?)"
    r"\s*(?:\^\{?[a-e]\}?|\s[a-e])?\s*$"
)


def _parse(cell: str, power: str = "") -> tuple[str, str]:
    """Return the value and sigma of a cell, checking its power of ten."""
    match = _CELL.match(cell)
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    if (match.group("power") or "") != power:
        msg = f"cell {cell!r} is not x 10^-{power}"
        raise ValueError(msg)
    if match.group("bracket"):
        return number(match.group("bracket")), ""
    return number(match.group("value")), number(match.group("sigma"))


# Table 3.2: (library row, molecule cell, reference) -> file tag, row label, the D_0
# power of ten.
TABLE_3_2_ROWS = (
    (1, "H_2", "[3.136]", "H2", "H2", "H2", ""),
    (4, "D_2", "[3.136]", "D2", "D2", "D2", ""),
    (
        13,
        "^{14}N_2",
        "[3.17a]",
        "N2",
        "N2_butcher",
        "two-constant fit, Eq. (3.10)",
        "6",
    ),
    (14, "", "[3.17a]", "N2", "N2_butcher", "three-constant fit, Eq. (3.9)", "6"),
    (15, "", "[3.89]", "N2", "N2_bendtsen", "14N2", "6"),
    (20, "^{16}O_2", "[3.17a]", "O2", "O2_butcher", "16O2", "6"),
    (21, "", "[3.43]", "O2", "O2_edwards", "16O2", "6"),
)


def _table_3_2() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_3_2)
    if grid[0] != ["Molecule", "", "B_0 [cm^{-1}]", "D_0 [cm^{-1}]", "Ref."]:
        msg = f"unexpected Table 3.2 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for row, molecule, ref, material, tag, label, power in TABLE_3_2_ROWS:
        cells = grid[row]
        if (cells[1], cells[4]) != (molecule, ref):
            msg = f"Table 3.2 row {row} is {cells!r}"
            raise ValueError(msg)
        printed = " | ".join(cell for cell in cells if cell)
        for constant, column, scale in (("B0", 2, ""), ("D0", 3, power)):
            value, sigma = _parse(cells[column], scale)
            files.setdefault(
                f"{material}/{STEM}table3_2_{constant}_{tag}.csv", []
            ).append(Point(label, label, value, sigma, printed))
    return files


# Table 3.3, CO2: the two lines of each cell as the image shows them (the ground state,
# note a, then the c levels of (01^1 0), note b).
TABLE_3_3_LINES = (
    (
        "0.3902634±0.0000044",
        f"(1.270±0.035){TIMES}10^{{-7}}",
        "[3.141]",
        "ground",
        "(00^0 0)",
    ),
    (
        "0.3906398±0.0000082",
        f"(1.419±0.048){TIMES}10^{{-7}}",
        "[3.141]",
        "0110c",
        "(01^1 0) c",
    ),
)


def _table_3_3() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_3_3)
    cells = grid[1]
    if cells[:2] != ["Carbon Dioxide", "CO_{2}(00^{0}0)^{a}(01^{1}0)^{b}"]:
        msg = f"unexpected Table 3.3 row {cells!r}"
        raise ValueError(msg)
    for column in (2, 3, 4):
        joined = "".join(line[column - 2] for line in TABLE_3_3_LINES)
        if joined != cells[column]:
            msg = f"Table 3.3 column {column}: {cells[column]!r} != {joined!r}"
            raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for b0, d0, ref, tag, label in TABLE_3_3_LINES:
        printed = f"CO_2 {label} | {b0} | {d0} | {ref}"
        for constant, cell, power in (("B", b0, ""), ("D", d0, "7")):
            value, sigma = _parse(cell, power)
            files[f"CO2/{STEM}table3_3_{constant}_{tag}_CO2.csv"] = [
                Point(label, label, value, sigma, printed)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the files of Tables 3.2 and 3.3."""
    return _table_3_2() | _table_3_3()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

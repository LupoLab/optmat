# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper, Guo, Burton & Brion (1992): Ar, Kr, Xe oscillator strengths.

Continuum (low-resolution dipole (e,e), 1 eV FWHM): Tables I (argon,
``tbl_7266df86dba6`` and continued ``tbl_41f17cecf0c6``), II (krypton,
``tbl_bd52d819c303`` and ``tbl_5de368f90d0b``) and III (xenon, ``tbl_70c8cd79d2c0``
and ``tbl_5469ba982e92``) list photon energy (eV) and df/dE in ``10^-2 eV^-1`` in four
column pairs per page, read column by column (images checked: the cells agree). One
file per gas; the rows are numbered in reading order (``1`` ... ``N``) because the
argon continuation prints a misprinted energy, "208.0" for 308.0 (image checked;
corrected by an erratum), which would collide with the true 208.0 of the first page
as an energy-valued row id. The energy of the normalisation point carries a footnote
mark (``21.2^a``, ``21.0^a``, ``16.8^a``: "Normalized to Ref. [28] at 21.218 eV" or
16.848 eV) that is stripped.

Discrete lines (high-resolution dipole (e,e), 0.048 eV FWHM, and literature values):
Tables IV/V (argon), VI/VII (krypton), VIII/IX (xenon), cells files only (images
checked). The column headers give each transition and its energy in eV (from Moore's
tables); the header row of sub-table (b) of Tables V, VII and IX sits inside the table
and is detected as a row whose cells contain "eV". Each source row (the present work,
other experiments, theories) becomes one file per gas, combining the two tables of the
gas; rows are the transition energies (eV, as printed), the per-cell uncertainties in
parentheses ``0.0662(0.0033)`` go to the sigma column. A theory row appears as a
"Theory" and an "Experiment" row for Geiger [82] (semiempirical calculation and
electron impact), so the section is part of the file name. Handled structurally:

- the "Oscillator strength ratio" column (a_1/a_2 etc.) is derived and skipped;
- the "Total to ionization" column (sub-table (b) of V, VII, IX) goes to a separate
  ``*_total`` file with an index abscissa (it is the whole discrete spectrum up to the
  2P_3/2 threshold, not a line);
- values marked ^b in Tables IV, VI, VIII are the summed a_1 + a_2 (b_1 + b_2, c_1 +
  c_2): the image prints them centred between the two columns; the cells put them
  under the first. They go to a ``*_sum`` file with an index abscissa;
- "(a) dipole length", "(b) dipole velocity", "(a) wave function", "(b)
  semiempirical" rows are variants of the source named on the row above;
- rows holding only a method in parentheses ("(Electron impact)") are skipped; the
  method is in the source record;
- upper limits ("<0.001", "<0.0025") are not numbers and are skipped (listed in the
  record);
- footnote marks on values (``0.395^c``, ``1.640^d``) are stripped; the record says
  what they mean.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1992_AbsoluteOpticalOscillator"
CONTINUUM = {
    "Ar": ("tbl_7266df86dba6", "tbl_41f17cecf0c6"),
    "Kr": ("tbl_bd52d819c303", "tbl_5de368f90d0b"),
    "Xe": ("tbl_70c8cd79d2c0", "tbl_5469ba982e92"),
}
DISCRETE = {
    "Ar": ("tbl_8988f801b4d2", "tbl_699514709b3b"),
    "Kr": ("tbl_85dda44edbca", "tbl_0dd5fa4f52b5"),
    "Xe": ("tbl_ae92208a90e7", "tbl_1dbc4809ef5a"),
}
# Reference number -> file stem (first author and year).
SOURCES = {
    1: "cooper_1962",
    35: "amusya_1990",
    36: "knox_1958",
    37: "dow_knox_1966",
    38: "gruzdev_1967",
    39: "kim_1968",
    40: "aymar_1970",
    41: "lee_lu_1973",
    42: "lee_1974",
    43: "gruzdev_loginov_1975a",
    44: "gruzdev_loginov_1975b",
    45: "albat_1975",
    46: "stewart_1975",
    47: "aymar_coulombe_1978",
    52: "de_jongh_1971",
    53: "westerveld_1979",
    55: "tsurubuchi_1990",
    56: "wilkinson_1965",
    57: "wilkinson_1966",
    58: "griffin_1969",
    59: "chashchina_1967",
    61: "anderson_1965",
    62: "turner_1965",
    63: "morack_1967",
    64: "lawrence_1968",
    65: "irwin_1973",
    66: "wieme_1973",
    67: "matthias_1977",
    68: "chornay_1984",
    70: "lewis_1967",
    72: "copley_1974",
    73: "vallee_1977",
    74: "kramer_1984",
    75: "ferrell_1987",
    76: "mcconkey_1973",
    78: "chamberlain_1965",
    79: "geiger_1970",
    81: "delage_1976",
    82: "geiger_1977",
    83: "geiger_unpublished",
    84: "li_1988",
    85: "takayanagi_1990",
    86: "suzuki_1991",
    87: "natali_unpublished",
    96: "kuyatt_1975",
    97: "wiese_1969",
}
VARIANTS = {
    "(a) dipole length": "length",
    "(b) dipole velocity": "velocity",
    "(a) wave function": "wave_function",
    "(b) semiempirical": "semiempirical",
}
_ENERGY = re.compile(r"\((\d+\.\d+)[^()]*eV")
_REFERENCE = re.compile(r"\[(\d+)\]")
_VALUE = re.compile(
    r"^(?P<value>\d*\.\d+)(?:\((?P<sigma>\d*\.\d+)\))?(?P<mark>\^?[a-d])?$"
)
_MARK = re.compile(r"\^?[a-z]$")


def _continuum(gas: str) -> list[Point]:
    points = []
    for table_id in CONTINUUM[gas]:
        rows = table(DIRECTORY, table_id)[1:]
        for column in range(0, 8, 2):
            for cells in rows:
                energy, strength = cells[column].strip(), cells[column + 1].strip()
                if not energy and not strength:
                    continue  # the short last columns of the continuations
                printed = f"{energy} | {strength}"
                energy = _MARK.sub("", energy)
                index = str(len(points) + 1)
                points.append(
                    Point(index, number(energy), number(strength), "", printed)
                )
    return points


def _source(label: str, parent: str, section: str) -> str:
    """Return the file stem of a source row."""
    if label in VARIANTS:
        return f"{parent}_{VARIANTS[label]}"
    if label.startswith("Present work"):
        return "present"
    found = _REFERENCE.search(label)
    if found is None:
        msg = f"no reference in {label!r}"
        raise ValueError(msg)
    stem = SOURCES[int(found[1])]
    return f"{stem}_theory" if section == "Theory" and stem == "geiger_1977" else stem


def _header(cells: list[str], state: dict[str, dict[int, str]]) -> bool:
    """Read a header or title row into ``state``; return whether the row was one."""
    if not cells[0].strip() and any("eV" in c for c in cells[1:]):
        state["energies"], state["labels"] = {}, {}
        for column, cell in enumerate(cells[1:], start=1):
            found = _ENERGY.search(cell)
            if found:
                state["energies"][column], state["labels"][column] = found[1], cell
        return True
    title = " ".join(cells[1:])
    if "Oscillator strength" in title or "ratio" in title or "Total to" in title:
        for column, cell in enumerate(cells[1:], start=1):
            if "ratio" in cell:
                state["special"][column] = "ratio"
            elif "Total to" in cell:
                state["special"][column] = "total"
        return True
    return False


def _point(
    gas: str, stem: str, cell: str, column: int, cells: list[str], state: dict
) -> tuple[str, Point]:
    """Return the file and point of one value cell of a source row."""
    match = _VALUE.match(cell.replace(" ", ""))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, sigma = match["value"], match["sigma"] or ""
    printed_row = " | ".join(cells)
    prefix = f"{gas}/chan_1992__{stem}"
    if state["special"].get(column) == "total":
        return f"{prefix}_total_{gas}.csv", Point(
            "total", "1", value, sigma, printed_row
        )
    if (match["mark"] or "").endswith("b") and state["first"]:  # summed, note b
        return f"{prefix}_sum_{gas}.csv", Point("sum", "1", value, sigma, printed_row)
    energy = state["energies"][column]
    printed = f"{state['labels'][column]} | {printed_row}"
    return f"{prefix}_{gas}.csv", Point(energy, energy, value, sigma, printed)


def _discrete(gas: str) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id in DISCRETE[gas]:
        state: dict = {"special": {}, "first": table_id == DISCRETE[gas][0]}
        section, parent = "", ""
        for cells in json_table(DIRECTORY, table_id):
            label = cells[0].strip()
            if _header(cells, state):
                continue
            if label in {"Theory", "Experiment"}:
                section = label
                continue
            if not any(c.strip() for c in cells[1:]):
                if not label.startswith("("):  # a parent of (a)/(b) variant rows
                    parent = SOURCES[int(_REFERENCE.findall(label)[0])]
                continue  # or a method line
            stem = _source(label, parent, section)
            for column in range(1, len(cells)):
                cell = cells[column].strip()
                if not cell or cell.startswith("<"):
                    continue  # empty, or an upper limit
                if state["special"].get(column) == "ratio":
                    continue
                path, point = _point(gas, stem, cell, column, cells, state)
                files.setdefault(path, []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the continuum df/dE and the discrete oscillator strengths per gas."""
    files: dict[str, list[Point]] = {}
    for gas in CONTINUUM:
        files[f"{gas}/chan_1992__continuum_{gas}.csv"] = _continuum(gas)
        files |= _discrete(gas)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Olney et al. (1997): dipole sums from dipole (e,e) spectra.

Table 3 (four parts, ``tbl_e403ff88c1f4``, ``tbl_55a125f18af7``, ``tbl_ca4fabc433cb``,
``tbl_569d65b53314``; cells files only, images checked): the dipole sums S(u), u = -1,
-2, -3, -4, -5, -6, -8, -10, in atomic units, of the refined dipole (e,e) spectra
("Expt") and of other work (DOSD, empirical, theory), one row per source with the
reference number in brackets. Only the first two parts hold curated species (He, Ne,
Ar, Kr, Xe, H2, D2, N2, O2, CO, CO2, N2O, CH4); parts 3 and 4 (the larger alkanes,
halides, amines, silanes, phosphines) are read by nobody here. HD, NO, NO2, H2O, H2S
and the hydrocarbons beyond CH4 are skipped. Handled structurally:

- the gas is printed on the first row of its block only (the other rows' first cell
  is empty under the row span): it is carried down;
- the Expt row's S(-2) cell carries the reference of the polarizability the spectrum
  was renormalised to, ``1.383 [79]``: it is stripped (the record says which it is)
  and the S(-2) goes to its own file, ``s-2_<gas>``, because it is a reproduction of
  that polarizability, not an independent sum; the other orders of the Expt row go to
  ``moments_<gas>``;
- powers of ten are printed in parentheses, ``2.404(4)`` = 2.404e4 (the table's own
  convention: S(-10) of O2 lies between Ar's 1800 and Xe's 5.275(4); not an
  uncertainty); written in e-notation, the cell kept in ``printed``;
- empty cells (orders a source does not give) are skipped.

Table 1 (``tbl_96bf5b17dd16``, ``tbl_82467dffab5d``; cells files only, images
checked), column 2 "From dipole (e,e) This work": the static polarizability (au) from
the originally published, valence-shell TRK-normalised spectra, before the refinement
of Table 3, with the reference of the spectrum in brackets (stripped). The CO cell
holds two values, ``12.64 [8]13.18^{aa}``: the published spectrum's and, note aa, the
value after correcting its energy calibration; they go to two files. The other columns
of Table 1 (refractive-index, dielectric, DOSD, theory and additivity values of other
authors, mostly from compilations) are not extracted.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Olney_1997_AbsoluteScaleDetermination"
TABLE_3 = ("tbl_e403ff88c1f4", "tbl_55a125f18af7")
TABLE_1 = ("tbl_96bf5b17dd16", "tbl_82467dffab5d")
GASES = {
    "He": "He",
    "Ne": "Ne",
    "Ar": "Ar",
    "Kr": "Kr",
    "Xe": "Xe",
    "H_2": "H2",
    "D_2": "D2",
    "N_2": "N2",
    "O_2": "O2",
    "CO": "CO",
    "CO_2": "CO2",
    "N_2O": "N2O",
    "CH_4": "CH4",
}
# Table 3's columns 3-10.
ORDERS = ("-1", "-2", "-3", "-4", "-5", "-6", "-8", "-10")
# Reference number of a comparison row -> file stem (first author and year).
SOURCES = {
    "[49]": "zeiss_1977",
    "[50]": "thomas_1977",
    "[52]": "jhanwar_1982",
    "[54]": "kumar_1985",
    "[56]": "kumar_1994",
    "[57]": "dehmer_1975",
    "[58]": "cummings_1975",
    "[59]": "inokuti_1981",
    "[60]": "geertsen_1986",
    "[61]": "packer_1994",
    "[62]": "kosman_1985",
    "[63]": "langhoff_1969",
    "[64]": "starkschall_1971",
    "[65]": "parker_1976",
    "[66]": "pack_1974",
    "[67]": "hohm_1993",
    "[68]": "victor_1969",
    "[69]": "leonard_1975",
    "[70]": "gerhart_1975",
    "[71]": "eggarter_1975",
    "[136]": "thakkar_1981",
    "[161]": "sauer_1990",
}
_CELL = re.compile(r"^(?P<mantissa>\d+\.\d+|\d+)(?:\((?P<power>\d+)\))?$")
_REFERENCE = re.compile(r"\s*\[[\d,]+\]$")


def _number(cell: str) -> str:
    """Return a printed cell's number, ``mantissa(power of ten)`` in e-notation."""
    match = _CELL.match(cell)
    if match is None:
        msg = f"not a number: {cell!r}"
        raise ValueError(msg)
    power = match["power"]
    return f"{match['mantissa']}e{power}" if power else match["mantissa"]


def _table_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id in TABLE_3:
        gas = ""
        for cells in json_table(DIRECTORY, table_id)[1:]:
            gas = cells[0] or gas  # the gas is printed on its block's first row only
            if gas not in GASES:
                continue
            material = GASES[gas]
            kind, reference = cells[1], cells[2]
            printed = " | ".join(cells)
            for order, cell in zip(ORDERS, cells[3:11], strict=True):
                if not cell:
                    continue
                if kind == "Expt":
                    # The refined spectrum's S(-2) is the polarizability it was
                    # renormalised to (Table 3 note a): its own file, reference
                    # stripped.
                    stem = "s-2" if order == "-2" else "moments"
                    value = _number(_REFERENCE.sub("", cell))
                else:
                    stem = SOURCES[reference]
                    value = _number(cell)
                point = Point(order, order, value, "", printed)
                files.setdefault(f"{material}/olney_1997__{stem}_{material}.csv", [])
                files[f"{material}/olney_1997__{stem}_{material}.csv"].append(point)
    return files


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for table_id in TABLE_1:
        for cells in json_table(DIRECTORY, table_id)[3:]:
            if cells[0] not in GASES:
                continue
            material = GASES[cells[0]]
            printed = " | ".join(cells[:2])
            if material == "CO":
                # "12.64 [8]13.18^{aa}": the published spectrum's value, then (note aa)
                # the value after the energy-scale correction of Section 3.1.
                published, corrected = cells[1].split("]")
                values = {
                    "table1_CO_published": _REFERENCE.sub("", published + "]"),
                    "table1_CO": corrected.removesuffix("^{aa}"),
                }
            else:
                values = {f"table1_{material}": _REFERENCE.sub("", cells[1])}
            for stem, value in values.items():
                point = Point("static", "inf", _number(value), "", printed)
                files[f"{material}/olney_1997__{stem}.csv"] = [point]
    return files


def extract() -> dict[str, list[Point]]:
    """Return every data file of the source, keyed by its ``<gas>/<file>`` name."""
    return _table_3() | _table_1()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hohm (1994), Mol. Phys. 81, 157: polarizabilities and B_R at 325-633 nm.

Table 1 (``tbl_43f3f21f6009``, no CSV rendering: read from its JSON cells, checked
against the image) gives the mean polarizability alpha(omega_k, T) in atomic units at
the five vacuum wavelengths of the table's note a) (632.99, 611.971, 594.096, 543.516
and 325.130 nm; the text, p. 159, prints 632.990). Each gas has two rows: this work
(first cell marked ``^b``; for CH4 and the hydrocarbons ``^h``, ref. [4], Hohm's 1993
paper) with the absolute uncertainty of the last digits in parentheses, and a row
calculated from the Cauchy moments of a DOSD (``^c`` [21] Kumar and Meath 1985, ``^d``
[22] Zeiss et al. 1977, ``^e`` [23] Thomas and Meath 1977, ``^j`` [28] Jhanwar and
Meath 1982; the others are not curated species). The decimal point is printed "·";
"—" is no value.

Table 2 (``tbl_83390c4d1bfd``, CSV) gives B_R(omega_k, T) in 1e-12 m^6 mol^-2
(= cm^6/mol^2) at omega_1, omega_3, omega_4 and omega_5. The omega_1 column is the
reference value B_R^(0) that the others are measured against (cell marker ``^a`` [8]
Burns et al. 1986, ``^b`` [12] Hohm's heuristic estimate); the omega_3-omega_5 columns
are this work. An uncertainty in parentheses without a decimal point is in the last
digits of the value ("1·57(58)" is 1.57 +- 0.58); one with a decimal point is absolute
("6·23(1·55)"); "(-)" is none. The gas cells carry footnote marks (``Kr^f``,
``Xe^{(1)}``, ``CO_{2}^{(2)}``: two independent measurements). The last column, B*,
is a fit constant and is not extracted. The library swaps the table's caption and its
footnotes (the footnotes appear as the caption): the image (p. 162) shows the caption
above the table and the footnotes a-f below it; nothing in the cells depends on it.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Hohm_1994_FrequencyDependenceSecond"
TABLE_1 = "tbl_43f3f21f6009"
TABLE_2 = "tbl_83390c4d1bfd"
WAVELENGTHS = ("632.990", "611.971", "594.096", "543.516", "325.130")  # nm, vacuum
TABLE_2_COLUMNS = {2: "632.990", 3: "594.096", 4: "543.516", 5: "325.130"}
GASES_1 = {"Ar": "Ar", "Kr": "Kr", "Xe": "Xe", "O_2": "O2", "N_2O": "N2O"}
GASES_1 |= {"CO_2": "CO2", "CH_4": "CH4"}
# The calculated row's marker names its DOSD (note b of Table 1).
DOSD = {"c": "kumar_1985", "d": "zeiss_1977", "e": "thomas_1977", "j": "jhanwar_1982"}
GASES_2 = {
    "Ar": "Ar",
    "Kr^f": "Kr",
    "Xe^{(1)}": "Xe",
    "Xe^{(2)}": "Xe",
    "O_2": "O2",
    "N_2O": "N2O",
    "CO_{2}^{(1)}": "CO2",
    "CO_{2}^{(2)}": "CO2",
    "CH_4": "CH4",
}
RUN = {"Xe^{(1)}": "_1", "Xe^{(2)}": "_2", "CO_{2}^{(1)}": "_1", "CO_{2}^{(2)}": "_2"}
_MARK = re.compile(r"^\^(?:\{(?P<b>[a-z])\}|(?P<p>[a-z]))")
_CELL = re.compile(r"^(?P<value>-?[\d.]+)(?:\((?P<sigma>[\d.]+|-)\))?$")
_NONE = {"—", "-", ""}


def _mark(cell: str) -> tuple[str, str]:
    """Split a leading footnote mark (``^b``, ``^{a}``) from a cell."""
    text = cell.strip()
    match = _MARK.match(text)
    if match is None:
        return "", text
    return match["b"] or match["p"], text[match.end() :]


def _value(cell: str) -> tuple[str, str]:
    """Return a cell's value and absolute sigma, as printed ("·" read as ".")."""
    text = cell.replace("·", ".").strip()
    match = _CELL.match(text)
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value = number(match["value"])
    digits = match["sigma"] or ""
    if digits in {"", "-"}:
        return value, ""
    if "." in digits:  # an absolute uncertainty, "6·23(1·55)"
        return value, digits
    decimals = len(value.partition(".")[2])  # last-digits uncertainty, "1·57(58)"
    return value, f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _alpha(cells: list[str], printed: str) -> list[Point]:
    """Return the five polarizability cells of a Table 1 row, skipping "—"."""
    points = []
    for wavelength, cell in zip(WAVELENGTHS, cells[2:7], strict=True):
        _, text = _mark(cell)
        if text in _NONE:
            continue
        value, sigma = _value(text)
        points.append(Point(wavelength, wavelength, value, sigma, printed))
    return points


def _table_1() -> dict[str, list[Point]]:
    """Return this work's and the DOSD rows of Table 1, per gas."""
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_1)[1:]
    for first, second in zip(rows[::2], rows[1::2], strict=True):
        if first[0] not in GASES_1:
            continue  # NH3 and the hydrocarbons other than CH4 are not curated
        gas = GASES_1[first[0]]
        mark, _ = _mark(first[2])
        if mark not in {"b", "h"}:
            msg = f"Table 1 {gas}: measured row marked {mark!r}"
            raise ValueError(msg)
        files[f"{gas}/hohm_1994b__alpha_{gas}.csv"] = _alpha(first, " | ".join(first))
        mark, _ = _mark(second[2])
        name = f"{gas}/hohm_1994b__{DOSD[mark]}_{gas}.csv"
        files[name] = _alpha(second, " | ".join(second))
    return files


def _table_2() -> dict[str, list[Point]]:
    """Return B_R of Table 2: the omega_1 reference and this work's omega_3-omega_5."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        if cells[0] not in GASES_2:
            continue  # NH3 and the hydrocarbons other than CH4 are not curated
        gas = GASES_2[cells[0]]
        run = RUN.get(cells[0], "")
        printed = " | ".join(cells)
        points = []
        for column, wavelength in TABLE_2_COLUMNS.items():
            mark, text = _mark(cells[column])
            if text in _NONE:
                continue
            value, sigma = _value(text)
            point = Point(wavelength, wavelength, value, sigma, printed)
            if column == 2:  # the reference, once per gas (Xe and CO2 repeat it)
                if mark not in {"a", "b"}:
                    msg = f"Table 2 {gas}: reference marked {mark!r}"
                    raise ValueError(msg)
                files.setdefault(f"{gas}/hohm_1994b__B_R_reference_{gas}.csv", [point])
            else:
                points.append(point)
        files[f"{gas}/hohm_1994b__B_R{run}_{gas}.csv"] = points
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table 1 polarizabilities and the Table 2 B_R, per gas."""
    return _table_1() | _table_2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

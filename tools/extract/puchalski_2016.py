# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Puchalski et al. (2016): helium polarizability, dispersion coefficients.

- Table II (``tbl_047bf9b3647c``): the nonrelativistic dispersion coefficients
  alpha_k^(0), k = 0, 2, 4, 6, of 4He in atomic units (a0^3 t0^k, not reduced units),
  "Present" with the last-digit uncertainty in parentheses, and Bhatia and Drachman's
  values converted by the authors (column "BD", footnote letters a and b stripped).
  alpha_k = 2 <z R^(k+1) z> = S(-k-2) (Eq. 10 and the sum over states), so each row is
  written as the Cauchy moment S(-k-2), the order in the wavelength column.
- Table VII (``tbl_f97b1bb8eafc``, cells file only): the static total and the total
  alpha(6329.908 A), the value of Ref. [13] (Piszczatowski et al. 2015) and the value
  the authors compute from Bhatia and Drachman's correctly converted coefficients. The
  component rows are not extracted (they are in the source record's formulas).
- Table IX (``tbl_a088e2718911``, cells file only): molar polarizabilities
  A_eps = 4 pi N_A alpha / 3 in cm^3/mol, static and at 6329.908 A, present and
  reproduced.

The table images agree with the cells: the JSON cells drop the thin spaces of
"0.517 254 19(10)" (Table IX, Ref. [2]) and write the footnote marks as ``^{a}``.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Puchalski_2016_TheoreticalDeterminationPolarizability"
TABLE_II = "tbl_047bf9b3647c"
TABLE_VII = "tbl_f97b1bb8eafc"
TABLE_IX = "tbl_a088e2718911"
PREFIX = "He/puchalski_2016__"
_VALUE = re.compile(r"^(-?[\d.]+)(?:\((\d+)\))?(?:\^\{?[a-f]\}?|[a-f])?$")
# Table II rows alpha_k^(0) -> Cauchy moment order -k-2.
ALPHA = "\N{GREEK SMALL LETTER ALPHA}"
ORDERS = {f"{ALPHA}{k}(0)": str(-k - 2) for k in (0, 2, 4, 6)}


def split(cell: str) -> tuple[str, str]:
    """Return the value and the last-digits uncertainty of ``1.234(5)^{a}``.

    The uncertainty is returned in the value's unit (``1.234(5)`` gives ``0.005``),
    or "" when none is printed; footnote marks are dropped.
    """
    text = cell.strip().replace(" ", "").replace("\N{THIN SPACE}", "")
    match = _VALUE.match(text)
    if match is None:
        msg = f"cannot read {cell!r}"
        raise ValueError(msg)
    value, digits = match.group(1), match.group(2)
    if digits is None:
        return number(value), ""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def _table_ii() -> dict[str, list[Point]]:
    present, bhatia = [], []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        order = ORDERS[cells[0].strip()]
        printed = " | ".join(cells)
        value, sigma = split(cells[2])
        present.append(Point(order, order, value, sigma, printed))
        value, _ = split(cells[1])
        bhatia.append(Point(order, order, value, "", printed))
    return {
        f"{PREFIX}cauchy_nonrelativistic.csv": present,
        f"{PREFIX}cauchy_bhatia_drachman.csv": bhatia,
    }


def _table_vii() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_VII)
    # Row 7 is the static "total"; rows 21-23 the totals at 6329.908 A.
    labels = [(r[0], r[1]) for r in rows]
    assert labels[7] == ("", "total"), labels[7]
    assert labels[21] == ("Total", "present"), labels[21]
    assert labels[22][1] == "Ref. [13]", labels[22]
    assert labels[23][1] == "BD^{f}", labels[23]

    def point(row: int, name: str, wavelength: str) -> Point:
        value, sigma = split(rows[row][2])
        return Point(name, wavelength, value, sigma, " | ".join(rows[row]))

    return {
        f"{PREFIX}alpha.csv": [
            point(7, "static", "inf"),
            point(21, "6329.908", "6329.908"),
        ],
        f"{PREFIX}alpha_piszczatowski_2015.csv": [point(22, "6329.908", "6329.908")],
        f"{PREFIX}alpha_bhatia_drachman.csv": [point(23, "6329.908", "6329.908")],
    }


def _table_ix() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_IX)
    names = {
        (0, "present"): ("molar_polarizability", "static", "inf"),
        (1, "Ref. [2]"): ("molar_lach_2004", "static", "inf"),
        (2, "expt. [8]"): ("molar_schmidt_2007", "static", "inf"),
        (3, "present"): ("molar_polarizability", "6329.908", "6329.908"),
        (4, "expt. [41,42]"): ("molar_achtermann", "6329.908", "6329.908"),
        (5, "expt. [43]"): ("molar_birch_1991", "6329.908", "6329.908"),
    }
    files: dict[str, list[Point]] = {}
    for (index, label), (name, row, wavelength) in names.items():
        cells = rows[index]
        assert cells[1] == label, cells
        value, sigma = split(cells[2])
        point = Point(row, wavelength, value, sigma, " | ".join(cells))
        files.setdefault(f"{PREFIX}{name}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables II, VII (totals) and IX."""
    return _table_ii() | _table_vii() | _table_ix()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ford & Browne (1973): the dynamic polarizability of H2 by sum over states.

All tables were checked against their images; Tables II and VII are JSON only.

- Table VII (``tbl_8218ee3b07f5``): alpha(omega) of H2, a0^3, at 24 vacuum wavelengths
  in angstroms (the static limit printed as an infinity sign, stored as ``inf``):
  "This work" (thermal average at 293 K, footnote a) and the "Semi-empirical" column
  (Victor & Dalgarno), one file each. The "Expt." column (Victor & Dalgarno's
  compilation, and 12.2 at 1215.7 A from Gill & Heddle re-reduced with Ford's
  depolarization ratio) and the gamma(omega) block are not extracted.
- Table VII continued (``tbl_f90a64648330``): the Cauchy coefficients alpha_k of
  alpha(omega) = sum_k alpha_k omega^(2k), k = 0-3, "This work" and "From experimental
  refractive index data". They are the dipole sums S(-2k - 2), stored as
  ``cauchy_moment`` with rows k' = -2, -4, -6, -8 (the moment order).
- Table II (``tbl_8a57694d1337``): the static alpha_00^JJ, v = 0, J = 0-4, "This work"
  (rows ``J0``..``J4``, all static), and its thermal average at 300 K (one file each).
  The J = 0 cells carry the unit "a_{0}^{3}", stripped. The KW column (Kolos &
  Wolniewicz 1969) and the anisotropies are not extracted.
- Table XI (``tbl_31dd70cf3437``): results at Lyman alpha (1215.7 A), thermal average
  over rotational and vibrational states at 300 K: alpha_00 and the Rayleigh
  cross sections Q_00^c (elastic) and Q_00^T (with the rotational Raman lines), one
  file each.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Ford_1973_DirectResolventOperator"
TABLE_II = "tbl_8a57694d1337"
TABLE_VII = "tbl_8218ee3b07f5"
TABLE_VII_CAUCHY = "tbl_f90a64648330"
TABLE_XI = "tbl_31dd70cf3437"
LYMAN_ALPHA = "1215.7"  # A, Table XI's caption


def _wavelength(cell: str) -> str:
    cell = cell.strip()
    return "inf" if cell == "\N{INFINITY}" else number(cell)


def _table_vii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_VII)
    # The alpha block runs from its header row to the gamma block's header row.
    alpha, gamma = "\N{GREEK SMALL LETTER ALPHA}(", "\N{GREEK SMALL LETTER GAMMA}("
    start = next(i for i, c in enumerate(grid) if c[0].startswith(alpha)) + 1
    end = next(i for i, c in enumerate(grid) if c[0].startswith(gamma))
    this_work, semiempirical = [], []
    for cells in grid[start:end]:
        wavelength = _wavelength(cells[0])
        printed = " | ".join(cells)
        this_work.append(Point(wavelength, wavelength, number(cells[1]), "", printed))
        if cells[2]:  # no semi-empirical value at 6943.0 and 4880.0 A
            semiempirical.append(
                Point(wavelength, wavelength, number(cells[2]), "", printed)
            )
    return {
        "H2/ford_1973__thermal_293K.csv": this_work,
        "H2/ford_1973__victor_dalgarno_semiempirical.csv": semiempirical,
    }


def _cauchy() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {
        "H2/ford_1973__cauchy_this_work.csv": [],
        "H2/ford_1973__cauchy_refractive.csv": [],
    }
    for cells in json_table(DIRECTORY, TABLE_VII_CAUCHY)[2:]:
        order = str(-2 * int(cells[0]) - 2)  # alpha_k = S(-2k - 2)
        printed = " | ".join(cells)
        for column, path in enumerate(files, start=1):
            files[path].append(Point(order, order, number(cells[column]), "", printed))
    return files


def _table_ii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_II)
    states = []
    average = []
    for cells in grid[3:9]:
        printed = " | ".join(cells)
        value = number(cells[1].removesuffix("a_{0}^{3}"))
        if cells[0].startswith("thermal average"):
            average.append(Point("300K", "inf", value, "", printed))
        else:
            j = cells[0].removeprefix("J = ")
            states.append(Point(f"J{int(j)}", "inf", value, "", printed))
    return {
        "H2/ford_1973__static_states.csv": states,
        "H2/ford_1973__static_300K.csv": average,
    }


_XI_ALPHA = re.compile(r"\\alpha_\{00\} = (?P<value>[\d.]+)a_\{0\}\^\{3\}")
_XI_Q = re.compile(
    r"Q_\{00\}\^\{(?P<kind>[cT])\} = (?P<mantissa>[\d.]+) \\times 10\^\{-24\}"
)


def _table_xi() -> dict[str, list[Point]]:
    files = {}
    for (cell,) in table(DIRECTORY, TABLE_XI):
        if alpha := _XI_ALPHA.search(cell):
            files["H2/ford_1973__lyman_alpha_300K.csv"] = [
                Point(LYMAN_ALPHA, LYMAN_ALPHA, alpha["value"], "", cell)
            ]
        elif q := _XI_Q.search(cell):
            stem = {"c": "rayleigh_elastic", "T": "rayleigh_total"}[q["kind"]]
            files[f"H2/ford_1973__{stem}.csv"] = [
                Point(LYMAN_ALPHA, LYMAN_ALPHA, q["mantissa"], "", cell)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables II, VII and XI's polarizabilities and cross sections."""
    return _table_vii() | _cauchy() | _table_ii() | _table_xi()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

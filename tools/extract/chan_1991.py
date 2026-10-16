# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper & Brion (1991): helium oscillator strengths by dipole (e,e).

- Table II (``tbl_bd7d4d0ec5c7``): the low-resolution (1 eV FWHM) absolute
  differential oscillator strengths df/dE (1e-2 /eV) at 24.6-180 eV, three column
  pairs read pair by pair (print order).
- Table III (``tbl_ae85cd966ffd``, cells file only): discrete oscillator strengths
  1 1S -> n 1P, n = 2-7, theory and experiment, one file per author; rows and the
  wavelength column are the principal quantum number n (the table gives no transition
  energies). Experimental cells print the uncertainty in parentheses, "0.280(0.007)".
  Korolyov and Odintsov have two n = 2 values on consecutive rows (the second row's
  label is empty, image checked): rows "2-1" and "2-2". The "Total to ionization
  threshold" column is not extracted (not a line). Printed misprints are corrected by
  errata in the source record, not here.

The images agree with the cells.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1991_AbsoluteOpticalOscillator"
TABLE_II = "tbl_bd7d4d0ec5c7"
TABLE_III = "tbl_ae85cd966ffd"
PREFIX = "He/chan_1991__"
LEVEL_CROSSING = "Hanle"  # codespell:ignore hanle
FRY_WILLIAMS = f"Fry and Williams (Ref. 12)(Lifetime: {LEVEL_CROSSING} effect)"
AUTHORS = {
    "Fernley, Taylor, and Seaton (Ref. 9)": "fernley_1987",
    "Schiff, Pekeris, and Accad (Ref. 8)": "schiff_1971",
    "Weiss (Ref. 7)": "weiss_1967",
    "Green, Johnson, and Kolchin (Ref. 5)": "green_1967",
    "Dalgarno and Parkinson (Ref. 6)": "dalgarno_parkinson_1967",
    "Wiese, Smith, and Glennon (Ref. 10)": "wiese_1966",
    "Schiff and Pekeris (Ref. 4)": "schiff_pekeris_1964",
    "Dalgarno and Stewart (Ref. 3)": "dalgarno_stewart_1960",
    "Present work [HR dipole (e,e)]": "present",
    "Tsurubuchi, Watanabe, and Arikawa (Ref. 20)(Self absorption)": "tsurubuchi_1989",
    "Westerveld and Van Eck (Ref. 19)(Self absorption)": "westerveld_1977",
    "Backx et al. (Ref. 47)(Electron impact)^a": "backx_1975",
    "Burger and Lurio (Ref. 13)(Lifetime: Level-crossing)": "burger_lurio_1971",
    "DeJongh and Van Eck (Ref. 18)(Self absorption)^b": "de_jongh_1971",
    "Lassetre, Skerbele, and Dillon (Ref. 30)(Electron impact)": "lassettre_1970",
    "Martinson and Bickel (Ref. 14)(Lifetime: Beam foil)": "martinson_1969",
    FRY_WILLIAMS: "fry_williams_1969",
    "Lincke and Griem (Ref. 15)(Plasmas emission profile)": "lincke_griem_1966",
    "Korolyov and Odintsov (Ref. 16)(Beam emission profile)": "korolyov_1965",
    "Kuhn and Vaughan (Ref. 17)(Resonance broadening emission profile)": "kuhn_1964",
    "Geiger (Ref. 63)(Electron impact)": "geiger_1963",
}
_CELL = re.compile(r"^([\d.]+)(?:\(([\d.]+)\))?$")


def _table_ii() -> list[Point]:
    rows = table(DIRECTORY, TABLE_II)[1:]
    return [
        Point(
            cells[2 * k].strip(),
            number(cells[2 * k]),
            number(cells[2 * k + 1]),
            "",
            f"{cells[2 * k].strip()} | {cells[2 * k + 1].strip()}",
        )
        for k in range(3)
        for cells in rows
        if cells[2 * k].strip()
    ]


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    name = ""
    for cells in json_table(DIRECTORY, TABLE_III)[3:]:
        if cells[0].startswith("("):
            continue  # the "(b) Experiment" heading
        second = not cells[0]  # Korolyov and Odintsov's second n = 2 row
        name = name if second else AUTHORS[cells[0]]
        printed = " | ".join(cells)
        for n, cell in zip(range(2, 8), cells[1:7], strict=True):
            if not cell:
                continue
            match = _CELL.match(cell)
            assert match is not None, cell
            value, sigma = match.group(1), match.group(2) or ""
            row = (
                f"{n}-2"
                if second
                else (f"{n}-1" if name == "korolyov_1965" else str(n))
            )
            point = Point(row, str(n), number(value), sigma, printed)
            files.setdefault(f"{PREFIX}f_{name}.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table II (df/dE) and Table III (discrete f)."""
    return {f"{PREFIX}df_dE.csv": _table_ii()} | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

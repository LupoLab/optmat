# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wahlstrand et al. (2018), supplement: the fits of the ionisation yield.

Table 1 of the Supplemental Material (``supplement_01``, ``tbl_668d3655bf2e``, its p.
11; the text calls it Table S1) lists, per gas (Ar, Kr, Xe, N2, O2), c_1 and m_fit of
the fit Y = c_1 I_0^m to the measured yield, Eq. (S1), the multiphoton order m_MPI, and
two columns of n_2. Per gas, one file of c_1 (``yield_fit_c1_<gas>``) and one of m
(``yield_fit_m_<gas>``). m_MPI is not a fitted value and the n_2 columns are not
ionisation data: they are kept in the ``printed`` cells only.

The table image was checked, every cell of the first four columns (2026-09-28). The
text layer garbles the header's symbols (the image has "c_1 [(cm^2/TW)^m]", "m_fit",
"m_MPI"; the library has the exponent of the unit as "-m" and both m headers alike),
so the columns are fixed by position and the header is checked as the library has it.
A value printed "2.05 x 10^-13" is "2.05\N{MULTIPLICATION SIGN}10-13" in the text layer
(the superscript is lost) and is written 2.05e-13.
"""

import re

from _library import Point, table, write

DIRECTORY = "Wahlstrand_2018_BoundElectronNonlinearity/supplement_01"
TABLE = "tbl_668d3655bf2e"
GASES = ("Ar", "Kr", "Xe", "N2", "O2")
COLUMNS = 6
_POWER = re.compile(
    r"^(?P<mantissa>\d+\.\d+)\N{MULTIPLICATION SIGN}10-(?P<exponent>\d+)$"
)
_PLAIN = re.compile(r"^\d+\.\d+$")


def extract() -> dict[str, list[Point]]:
    """Return c_1 and m of the yield fit per gas."""
    grid = table(DIRECTORY, TABLE)
    header = grid[0]
    if (
        len(header) != COLUMNS
        or header[0].strip() != "Gas"
        or "(cm2/TW)" not in header[1]
        or header[2] != header[3]
        or not header[4].startswith(header[5][:10])
        or tuple(row[0].strip() for row in grid[1:]) != GASES
    ):
        msg = f"unexpected Table 1 layout {grid!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        gas = cells[0].strip()
        coefficient = _POWER.match(cells[1].strip())
        exponent = cells[2].strip()
        if coefficient is None or _PLAIN.match(exponent) is None:
            msg = f"Table 1, {gas}: unreadable cells {cells[1]!r}, {cells[2]!r}"
            raise ValueError(msg)
        printed = " | ".join(cell.strip() for cell in cells)
        value = f"{coefficient['mantissa']}e-{coefficient['exponent']}"
        files[f"{gas}/wahlstrand_2018__yield_fit_c1_{gas}.csv"] = [
            Point("c1", "c1", value, "", printed)
        ]
        files[f"{gas}/wahlstrand_2018__yield_fit_m_{gas}.csv"] = [
            Point("m", "m", exponent, "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

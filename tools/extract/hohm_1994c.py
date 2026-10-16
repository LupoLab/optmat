# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hohm & Trümper (1994): the polarizability of xenon against temperature.

Both tables have CSV renderings; both images checked (no differences).

* Table 1 (``tbl_080bc44cea27``): alpha(omega_k, T) in atomic units at the four
  frequencies omega_1-omega_4 of the table's note ("omega_1 = 0.07198116 au, ...": the
  columns carry only their labels, so the frequencies are mapped here from the note, as
  printed) and the extrapolated alpha(omega -> 0), at nine temperatures. Each value is
  printed with its statistical rms error on the last digits ("27.808(17)"). One file per
  temperature for the four measured frequencies (``alpha_<T>_Xe``, rows and abscissa =
  omega in hartree, as printed in the note), and one file of the extrapolated static
  values (``alpha_static_T_Xe``, row = the temperature in K as printed without its own
  uncertainty, abscissa 0).
* Table 2 (``tbl_5107ed7c9cd1``): the fit alpha(omega_k, T) = b0(omega) + b1(omega) T^2
  per frequency. The b0 column (the fit's T = 0 intercept, with its uncertainty) is a
  dataset (``b0_Xe``); b1 and the theoretical b1 of Eq. (16) are in the record's
  formulas.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Hohm_1994_TemperatureDependenceDipole"
TABLE_1 = "tbl_080bc44cea27"
TABLE_2 = "tbl_5107ed7c9cd1"
# Table 1's note: "omega_1 = 0.07198116 au, omega_2 = 0.07669359 au, omega_3 =
# 0.08383075 au, and omega_4 = 0.1401389 au" (the columns alpha(omega_1) to
# alpha(omega_4)).
OMEGAS = ("0.07198116", "0.07669359", "0.08383075", "0.1401389")
_UNCERTAIN = re.compile(r"^([\d.]+)\(([\d.]+)\)$")


def _split(cell: str) -> tuple[str, str]:
    """Return "27.808(17)" as ("27.808", "0.017") and "3.91(1.2)" as ("3.91", "1.2")."""
    match = _UNCERTAIN.match(cell.strip())
    if match is None:
        msg = f"no uncertainty in {cell!r}"
        raise ValueError(msg)
    value, digits = match[1], match[2]
    if "." in digits:  # an uncertainty printed in the value's own unit
        return number(value), number(digits)
    decimals = len(value.partition(".")[2])
    return number(value), f"{int(digits) / 10**decimals:.{decimals}f}"


def _table1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_1)
    if grid[0][0] != "T(K)" or grid[0][5] != "\\alpha(\\omega\\to 0)":
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    static = []
    for cells in grid[1:]:
        temperature, _ = _split(cells[0])
        printed = " | ".join(cells)
        points = []
        for omega, cell in zip(OMEGAS, cells[1:5], strict=True):
            value, sigma = _split(cell)
            points.append(Point(omega, omega, value, sigma, printed))
        files[f"Xe/hohm_1994c__alpha_{temperature}_Xe.csv"] = points
        static.append(Point(temperature, "0", *_split(cells[5]), printed))
    files["Xe/hohm_1994c__alpha_static_T_Xe.csv"] = static
    return files


def _table2() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_2)
    if not grid[0][1].startswith("b_{0}"):
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[1:]:
        omega = number(cells[0])
        points.append(Point(omega, omega, *_split(cells[1]), " | ".join(cells)))
    return {"Xe/hohm_1994c__b0_Xe.csv": points}


def extract() -> dict[str, list[Point]]:
    """Return Tables 1 and 2."""
    return _table1() | _table2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

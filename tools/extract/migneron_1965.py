# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Migneron & Levinger (1965): helium n - 1 at 60-700 nm and sum rules S(p).

Table I (``tbl_a097df9cffbb``) is printed as two pairs of columns, wavelength (A) and
``n - 1`` with its factor ``x 10^-6``; both pairs are read. Table II
(``tbl_f60413766558``) gives the oscillator-strength sums S(p), p = -2..2, in units of
Ry^p, by energy region, their total, and "other calculation" values (footnote e:
Essen's and Johnston et al.'s measured polarizability for S(-2); f: Pekeris's
ground-state wave function for the others). Only the totals and the other values are
extracted, one file per p (the datasets convert Ry^p to atomic units with their scale,
2^-p); the footnote letters are stripped. The table images agree with the cells.
"""

from _library import Point, number, table, write

DIRECTORY = "Migneron_1965_IndexRefractionSum"
TABLE_I = "tbl_a097df9cffbb"
TABLE_II = "tbl_f60413766558"
PREFIX = "He/migneron_1965__"
# Table II column -> the order p and its file tag.
ORDERS = {2: ("-2", "m2"), 3: ("-1", "m1"), 4: ("0", "0"), 5: ("1", "1"), 6: ("2", "2")}


def _table_i() -> list[Point]:
    points = []
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        for pair in (0, 1):
            wavelength, value = cells[2 * pair].strip(), cells[2 * pair + 1].strip()
            if wavelength:
                points.append(
                    Point(wavelength, wavelength, number(value), "", " | ".join(cells))
                )
    return sorted(points, key=lambda p: float(p.x))


def _table_ii() -> dict[str, list[Point]]:
    rows = {cells[0].strip(): cells for cells in table(DIRECTORY, TABLE_II)[1:]}
    total, other = rows["Total"], rows["Other calculation"]
    files = {}
    for column, (order, tag) in ORDERS.items():
        files[f"{PREFIX}sum_rule_{tag}.csv"] = [
            Point(order, order, number(total[column]), "", " | ".join(total))
        ]
        value, _, mark = other[column].strip().partition("^")
        source = "essen" if mark == "e" else "pekeris"
        files[f"{PREFIX}{source}_sum_rule_{tag}.csv"] = [
            Point(order, order, number(value), "", " | ".join(other))
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table I and the totals and comparison values of Table II."""
    return {f"{PREFIX}theory.csv": _table_i()} | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

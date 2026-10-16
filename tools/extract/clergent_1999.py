# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Clergent, Durou & Laurens (1999): n - 1 of Ar, N2 and CO2 at 632.8 nm.

Table 5 (``tbl_4e3479144d9d``, CSV; image checked) compares (n - 1) x 10^6 at 101.325
kPa and 293.15 K: the row "Haufand Grigull" (Hauf and Grigull 1970, interpolated from
the International Critical Tables), the row "our results" (this work, from Table 4's mu
through Eq. 1) with its uncertainty, and the relative difference (not extracted).

Structural fix: the image prints "263.3 ± 0.8"; the CSV has "263.3 ( 0.8" (the "±"
extracted as "("). The cell is split at that character.

One file per gas and row (this work; Hauf and Grigull), each with one point at
632.8 nm.
"""

from _library import Point, number, table, write

DIRECTORY = "Clergent_1999_RefractiveIndexVariations"
TABLE_5 = "tbl_4e3479144d9d"
GASES = ("Ar", "N2", "CO2")


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``263.3 ( 0.8`` (printed 263.3 ± 0.8)."""
    value, sign, sigma = cell.partition("(")
    if not sign:
        msg = f"no uncertainty in {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return this work's and Hauf and Grigull's n - 1 per gas."""
    rows = table(DIRECTORY, TABLE_5)
    if rows[0] != ["", *GASES]:
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    hauf, ours = rows[1], rows[2]
    if hauf[0] != "Haufand Grigull" or ours[0] != "our results":
        msg = f"unexpected rows {hauf[0]!r}, {ours[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    printed_ours, printed_hauf = " | ".join(ours), " | ".join(hauf)
    for column, gas in enumerate(GASES, start=1):
        value, sigma = split(ours[column])
        files[f"{gas}/clergent_1999__refractivity_{gas}.csv"] = [
            Point("our results", "632.8", value, sigma, printed_ours)
        ]
        files[f"{gas}/clergent_1999__hauf_grigull_{gas}.csv"] = [
            Point("Hauf and Grigull", "632.8", number(hauf[column]), "", printed_hauf)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

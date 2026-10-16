# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Loriot et al. (2009): n2 of N2, O2 and Ar at 800 nm (Tables 1 and 2).

Table 2 (``tbl_732962a00e08``) prints the fitted coefficients of the intensity
expansion n_kerr = n2 I + n4 I^2 + ... + n10 I^5 per gas; only its n2 row is a
third-order value and is extracted (the higher orders are recorded in the source's
notes). Its n2 row repeats Table 1's "This work" row, with one more digit for Ar
(2.01 +- 0.19 against 2.0 +- 0.2); the paper says n2 was held at the Table 1 value
when fitting the higher orders (p. 5). The erratum (Opt. Express 18, 3011, 2010) divides
exactly this row by 2; the correction is applied as errata in the source record.

Table 1 (``tbl_90188e88dc6d``) also prints Nibbering et al.'s (1997) values in its
"n2 (Ref [2])" row: extracted as a reproduced dataset.

Both table images were checked: the CSV cells equal the printed ones. Air (the last
column) is not a curated species and is skipped.
"""

from _library import Point, number, table, write

DIRECTORY = "Loriot_2009_MeasurementHighOrder"
TABLE_1 = "tbl_90188e88dc6d"
TABLE_2 = "tbl_732962a00e08"
GASES = {"N2": 1, "O2": 2, "Ar": 3}  # column index in both tables; column 4 is air


def _plus_minus(cell: str) -> tuple[str, str]:
    """Split a printed ``a ± b`` cell into its two numbers, as printed."""
    value, sigma = cell.split("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return this work's n2 (Table 2) and Nibbering's (Table 1) per gas."""
    rows_2 = table(DIRECTORY, TABLE_2)
    (n2_row,) = [row for row in rows_2 if row[0].startswith("n2 ")]
    rows_1 = table(DIRECTORY, TABLE_1)
    (nibbering_row,) = [row for row in rows_1 if row[0].startswith("n2 (Ref")]
    files: dict[str, list[Point]] = {}
    for gas, column in GASES.items():
        value, sigma = _plus_minus(n2_row[column])
        files[f"{gas}/loriot_2009__n2_{gas}.csv"] = [
            Point("n2", "800", value, sigma, " | ".join(n2_row))
        ]
        value, sigma = _plus_minus(nibbering_row[column])
        files[f"{gas}/loriot_2009__nibbering_1997_{gas}.csv"] = [
            Point("n2_ref2", "800", value, sigma, " | ".join(nibbering_row))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

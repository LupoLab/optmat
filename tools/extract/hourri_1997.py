# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hourri, St-Arnaud & Bose (1997): the pure-gas rows of Table I.

Table I (``tbl_d52ac5f743cd``, CSV) lists A_eps (cm^3/mol), B_eps (cm^6/mol^2), C_eps
(cm^9/mol^3) and B_p of CO2-SF6 mixtures at 322.85 K by mole fraction of CO2. Only the
pure-gas rows are taken, both reprints from the literature: mole fraction 0 (SF6, A_eps
from Ref. 4, Hosticka and Bose 1974) and 1 (CO2, A_eps, B_eps and C_eps from Ref. 10,
Bose and Cole 1970). The mixtures are not curated species.

Structural fix: the CSV has lost the "±" and the minus signs. The image prints
"16.543±0.017^a", "7.350±0.006^c", "50.7±0.9^c" and "-2515±235^c"; the CSV has
"16.54360.017a", "7.35060.006c", "50.760.9c" and "225156235c" ("±" read as "6", the
minus as "2", the footnote letter run on). Each cell used is checked to be exactly the
garbled string and replaced by the image's reading, below.
"""

from _library import Point, table, write

DIRECTORY = "Hourri_1997_DielectricPressureVirial"
TABLE_I = "tbl_d52ac5f743cd"

# raw CSV cell -> (value, uncertainty, footnote) as the table image prints them.
IMAGE = {
    "16.54360.017a": ("16.543", "0.017", "a"),
    "7.35060.006c": ("7.350", "0.006", "c"),
    "50.760.9c": ("50.7", "0.9", "c"),
    "225156235c": ("-2515", "235", "c"),
}


def reading(cell: str) -> tuple[str, str]:
    """Return the image's value and uncertainty for a garbled cell.

    The reading only re-segments the cell: garbling it back ("-" to "2", "±" to "6",
    the footnote letter run on) must give the raw cell, digit for digit.
    """
    if cell not in IMAGE:
        msg = f"cell {cell!r} is not the one checked against the image"
        raise ValueError(msg)
    value, sigma, footnote = IMAGE[cell]
    if value.replace("-", "2") + "6" + sigma + footnote != cell:
        msg = f"the reading of {cell!r} changes its digits"
        raise ValueError(msg)
    return value, sigma


def extract() -> dict[str, list[Point]]:
    """Return SF6's A_eps and CO2's A_eps, B_eps and C_eps (reprints)."""
    rows = table(DIRECTORY, TABLE_I)
    sf6, co2 = rows[1], rows[7]
    if sf6[0] != "0" or co2[0] != "1":
        msg = f"unexpected pure rows {sf6[0]!r}, {co2[0]!r}"
        raise ValueError(msg)
    printed_sf6, printed_co2 = " | ".join(sf6), " | ".join(co2)
    return {
        "SF6/hourri_1997__A_eps_hosticka_SF6.csv": [
            Point("pure SF6", "inf", *reading(sf6[1]), printed_sf6)
        ],
        "CO2/hourri_1997__A_eps_bose_cole_CO2.csv": [
            Point("pure CO2", "inf", *reading(co2[1]), printed_co2)
        ],
        "CO2/hourri_1997__B_eps_bose_cole_CO2.csv": [
            Point("pure CO2", "inf", *reading(co2[2]), printed_co2)
        ],
        "CO2/hourri_1997__C_eps_bose_cole_CO2.csv": [
            Point("pure CO2", "inf", *reading(co2[3]), printed_co2)
        ],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bejot et al. (2010, hollow-core fibre): the n2 of Ar at 1.83 um (Table I).

Table I (``tbl_0df301bb6004``) prints, in its last row, n2 ... n10 and the fibre loss
alpha, each as ``<number>p`` with p = P / 1 bar; the rows above hold the powers of ten
and the units (the CSV keeps them as rows 1 and 2, as the image shows). Only n2 is a
third-order value. The trailing ``p`` is the pressure factor and is stripped; the value
is per bar.
"""

from _library import Point, number, table, write

DIRECTORY = "Bejot_2010_MechanismHollowCore"
TABLE = "tbl_0df301bb6004"


def extract() -> dict[str, list[Point]]:
    """Return the n2 column (per bar of Ar)."""
    header, powers, units, values = table(DIRECTORY, TABLE)
    assert header[0] == "n_{2}" and powers[0] == "10^{-24}", (header, powers)
    assert units[0].startswith(r"\text{m}^{2}"), units
    printed = values[0].strip()
    assert printed.endswith("p"), printed
    return {
        "Ar/bejot_2010_hollow__n2_Ar.csv": [
            Point("n2", "1.83", number(printed.removesuffix("p")), "", printed)
        ]
    }


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

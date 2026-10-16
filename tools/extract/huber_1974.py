# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Huber & Tondello (1974): helium n - 1 at 92-190 nm (Table I).

Table I gives, per wavelength, the experimental n - 1 at "NPT", a theoretical value
(Victor et al.) and their difference. Only the first row prints the factor
``\times 10^{-5}``; the rest print the mantissa. The theory column is extracted as a
separate, reproduced dataset.
"""

from _library import Point, number, table, write

DIRECTORY = "Huber_1974_RefractiveIndexHe"
TABLE = "tbl_6b16f76840f9"


def extract() -> dict[str, list[Point]]:
    """Return the experimental and the theoretical column as data files."""
    rows = table(DIRECTORY, TABLE)[1:]
    experiment = []
    theory = []
    for index, cells in enumerate(rows, start=1):
        printed = " | ".join(cells)
        experiment.append(Point(str(index), cells[0], number(cells[1]), "", printed))
        theory.append(Point(str(index), cells[0], number(cells[2]), "", printed))
    return {
        "He/huber_1974__experiment.csv": experiment,
        "He/huber_1974__victor_theory.csv": theory,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

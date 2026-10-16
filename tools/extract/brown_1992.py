# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Brown, Rich & Nibler (1992), J. Mol. Spectrosc. 151, 482: O2 S lines.

Table I (``tbl_ad5d63e4a91f``, CSV, p. 486; checked against its image): the Delta N =
Delta J = +2 (^SS) triplet lines S(J,N) of the pure rotational Raman spectrum of 16O2,
N = 1-19, as the peak maxima of the CARS spectrum ("Observed") and as the transition
wavenumbers from the spectral simulation ("Corrected"). One file per column; rows are
the printed labels S(J,N) (J and N of the lower level). The corrected - calculated
difference and the two intensity columns stay in ``printed``.
"""

from _library import Point, number, table, write

DIRECTORY = "Brown_1992_HighResolutionRotational"
TABLE_I = "tbl_ad5d63e4a91f"


def extract() -> dict[str, list[Point]]:
    """Return the observed and corrected line files of Table I."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0][:3] != ["S(J,N)", "Observedcm^{-1}", "Correctedcm^{-1}"]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        printed = " | ".join(cells)
        for column, name in ((1, "table1_observed"), (2, "table1_corrected")):
            files.setdefault(f"O2/brown_1992__{name}.csv", []).append(
                Point(cells[0], cells[0], number(cells[column]), "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

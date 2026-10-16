# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Mironov et al. (2025): argon's C_n*l in three approximations (Table III).

Table III (``tbl_bcc8c2bd802b``, CSV; image checked, p. 17) lists, per argon ion charge
0-17 and configuration, C_n*l "exact" (from the Hartree-Fock asymptotic tail), in the
Hartree approximation (Eq. (11)) and in the ADK approximation (l -> l* = n* - 1), with
the relative errors of the last two. Only charges 0 (Ar, 3p^6) and 1 (Ar+, 3p^5) are
curated materials; the coefficient is that of the 3p orbital, so the row is l = 1.
The neutral atom's Hartree value, "1.0*", 'is set manually' (p. 17) and is not a
calculated coefficient: it is not written. The relative errors are not written.
"""

from _library import Point, number, table, write

DIRECTORY = "Mironov_2025_StrongFieldIonization"
TABLE = "tbl_bcc8c2bd802b"
TARGETS = {"0": "Ar", "1": "Ar+"}
CONFIGS = {"0": "3p6", "1": "3p5"}
COLUMNS = {2: "exact", 3: "hartree", 5: "adk"}
WAVE = "1"  # the 3p orbital


def extract() -> dict[str, list[Point]]:
    """Return C_n*l of Ar and Ar+ in each approximation."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        charge = cells[0].strip()
        if charge not in TARGETS:
            continue
        if cells[1].strip() != CONFIGS[charge]:
            msg = f"unexpected configuration in {cells!r}"
            raise ValueError(msg)
        target = TARGETS[charge]
        printed = " | ".join(cells)
        for column, name in COLUMNS.items():
            if cells[column].strip().endswith("*"):  # set by hand, not computed
                continue
            files[f"{target}/mironov_2025__{name}_{target}.csv"] = [
                Point(WAVE, WAVE, number(cells[column]), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

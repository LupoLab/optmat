# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Awasthi et al. (2008): I_p and C_l of H2 at R = 1.4 a0 by method (Table I).

Table I (``tbl_4bbfbf50c79d``, CSV; image checked, p. 6): per electronic-structure
method (HF, HF [6], X-alpha [17], LB94, TSLDA, LDA, SAE-CI, CI) the vertical ionisation
potential (eV) and, for the independent-particle models, the coefficients C_0, C_2,
C_4 (m = 0). One file per method and quantity: ``ip_<method>`` and ``cl_<method>``;
blank cells (no C_l for LDA and the CI models, no C_4 for HF [6]) are not values.
"""

from _library import Point, number, table, write

DIRECTORY = "Awasthi_2008_SingleActiveElectron"
TABLE = "tbl_4bbfbf50c79d"
METHODS = {
    "HF": "hf",
    "HF [6]": "hf_kjeldsen",
    "X\\alpha[17]": "xalpha_tong",
    "LB94": "lb94",
    "TSLDA": "tslda",
    "LDA": "lda",
    "SAE-CI": "sae_ci",
    "CI": "ci",
}
WAVES = {2: "0", 3: "2", 4: "4"}


def extract() -> dict[str, list[Point]]:
    """Return I_p and C_l of each method."""
    files: dict[str, list[Point]] = {}
    grid = table(DIRECTORY, TABLE)
    if grid[0] != ["Method", "I_p(eV)", "C_0", "C_2", "C_4"]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[1:]:
        name = METHODS[cells[0].strip()]
        printed = " | ".join(cells)
        files[f"H2/awasthi_2008__ip_{name}.csv"] = [
            Point("I_p", "", number(cells[1]), "", printed)
        ]
        waves = [
            Point(wave, wave, number(cells[column]), "", printed)
            for column, wave in WAVES.items()
            if cells[column].strip()
        ]
        if waves:
            files[f"H2/awasthi_2008__cl_{name}.csv"] = waves
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

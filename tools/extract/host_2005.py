# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Host et al. (2005): coupled-cluster ESHG gamma of Ne, Ar and Kr.

Tables II (Ne, ``tbl_91f67ab85870``), VI (Ar, ``tbl_1441d7637357``) and XI (Kr,
``tbl_51b4b06ed8e7``), cells files only; images checked, every cell agrees. Columns:
gamma_par for CCS, CC2, CCSD, CC3 and "Expt.", then gamma_perp for CCS, CC2, CCSD, CC3;
rows omega = 0.0, 0.042 823, 0.088 559, 0.093 368 a.u. (the thin space in the printed
omega is dropped in the row id and wavelength). The experimental cells hold "value \pm
sigma" (or "value+-sigma", or a bare value) and are split into value and sigma.

The small-basis FCI benchmarks (Tables I, V, X), the basis-set studies (III, IV, VIII,
IX) and the relativistic corrections (VII, XII) are summarised in the record's notes,
not extracted.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Host_2005_FrequencyDependentHyperpolarizabilities"
TABLES = {"Ne": "tbl_91f67ab85870", "Ar": "tbl_1441d7637357", "Kr": "tbl_51b4b06ed8e7"}
HEADER = ["", "CCS", "CC2", "CCSD", "CC3", "Expt.^{a}", "CCS", "CC2", "CCSD", "CC3"]
COLUMNS = {
    1: "ccs_par",
    2: "cc2_par",
    3: "ccsd_par",
    4: "cc3_par",
    6: "ccs_perp",
    7: "cc2_perp",
    8: "ccsd_perp",
    9: "cc3_perp",
}
_EXPT = re.compile(r"^([\d.]+)\s*(?:(?:\\pm|±)\s*([\d.]+))?$")


def extract() -> dict[str, list[Point]]:
    """Return each method's gamma_par and gamma_perp, and the experimental column."""
    files: dict[str, list[Point]] = {}
    for gas, table_id in TABLES.items():
        grid = json_table(DIRECTORY, table_id)
        if grid[1] != HEADER:
            msg = f"unexpected header in {table_id}: {grid[1]!r}"
            raise ValueError(msg)
        prefix = f"{gas}/host_2005__"
        for cells in grid[2:]:
            omega = cells[0].replace(" ", "")
            float(omega)
            printed = " | ".join(cells)
            for column, name in COLUMNS.items():
                value = cells[column]
                float(value)
                files.setdefault(f"{prefix}{name}_{gas}.csv", []).append(
                    Point(omega, omega, value, "", printed)
                )
            expt = _EXPT.match(cells[5])
            if expt is None:
                msg = f"unparsed experimental cell {cells[5]!r}"
                raise ValueError(msg)
            files.setdefault(f"{prefix}expt_{gas}.csv", []).append(
                Point(omega, omega, expt.group(1), expt.group(2) or "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

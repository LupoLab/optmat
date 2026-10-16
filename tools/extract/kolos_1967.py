# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kolos & Wolniewicz (1967): static polarizabilities of H2 and D2 per state.

Both tables are JSON only (merged header cells) and were checked against their images.

- Table III (``tbl_008b76c50f4e``): the diagonal matrix elements alpha_{vJ,vJ} (a.u.)
  for v = 0, J = 0-31 and J = 0, v = 1-8, of H2, HD and D2 (with the anisotropies
  gamma). One file per gas, H2 and D2 (HD is not a curated species; gamma has no
  quantity); every value is static, so each row is a rovibrational state, named
  ``v<v>J<J>``, at wavelength ``inf``. D2's J = 14 prints "9.8344" (image checked): an
  erratum in the source record, not a structural fix.
- Table VI (``tbl_3187f42a4fee``): "This calculation" alpha_00 in cubic angstroms, two
  values in one merged cell, "0.8023" (the J = 0 state, Table III's 5.4139 a.u.) and
  "0.8045" with footnote d, "Rotational average for 293 K". Only the average is
  extracted (the J = 0 value repeats Table III). The Ishiguro et al. column and the
  experimental column are not extracted (source notes).
Tables I (convergence), II (alpha(R), clamped nuclei), IV (off-diagonal elements) and
V (isotope effects, relative differences) are not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Kolos_1967_PolarizabilityHydrogenMolecule"
TABLE_III = "tbl_008b76c50f4e"
TABLE_VI = "tbl_3187f42a4fee"
COLUMNS = {"H2": 2, "D2": 6}  # the alpha columns of Table III
# Table VI's merged cell: the J = 0 value then the 293 K average with footnote d.
_MERGED = re.compile(r"^(?P<j0>\d\.\d{4})(?P<average>\d\.\d{4})\^\{d\}$")


def extract() -> dict[str, list[Point]]:
    """Return the per-state static polarizabilities and the 293 K average."""
    rows = json_table(DIRECTORY, TABLE_III)[2:]
    files = {}
    for gas, column in COLUMNS.items():
        files[f"{gas}/kolos_1967__static_states.csv"] = [
            Point(f"v{c[0]}J{c[1]}", "inf", number(c[column]), "", " | ".join(c))
            for c in rows
        ]
    grid = json_table(DIRECTORY, TABLE_VI)
    (cells,) = [c for c in grid if c[0] == "\\alpha_{00}"]
    merged = _MERGED.match(cells[2])
    if merged is None:
        msg = f"unexpected Table VI cell {cells[2]!r}"
        raise ValueError(msg)
    files["H2/kolos_1967__normal_293K.csv"] = [
        Point("293K", "inf", merged["average"], "", " | ".join(cells))
    ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bishop & Pipin (1989): the dynamic second hyperpolarizability of helium.

Table II (``tbl_7d807148ca60``) and Table III (``tbl_ec80e5c576f2``), cells files only
(merged headers); both images checked, every cell agrees. Columns: gamma_zzzz for
"Kerr", "DFWM", "ESHG", "THG", then gamma_xxzz for "Kerr", "DFWM", "ESHG" (no THG
column: for THG gamma_xxzz = gamma_zzzz / 3). Bishop's "Kerr" is the dc Kerr effect
gamma(-w; w, 0, 0) and his "DFWM" is gamma(-w; -w, w, w), the intensity-dependent index
(p. 3549): they are written as the processes ``dc_kerr`` and ``kerr``.

Table II's first column is the wavelength in angstrom (used by the authors as a vacuum
wavelength, see the record); Table III's is omega in atomic units (hartree). The
static row (omega = 0.0) is kept in each dataset.
"""

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Bishop_1989_ImprovedDynamicHyperpolarizabilities"
TABLE_II = "tbl_7d807148ca60"
TABLE_III = "tbl_ec80e5c576f2"
PREFIX = "He/bishop_1989__"
# Column index -> dataset name (process as in the vocabulary, then the component).
COLUMNS = {
    1: "dc_kerr_zzzz",
    2: "kerr_zzzz",
    3: "eshg_zzzz",
    4: "thg_zzzz",
    5: "dc_kerr_xxzz",
    6: "kerr_xxzz",
    7: "eshg_xxzz",
}
HEADER = ["", "Kerr", "DFWM", "ESHG", "THG", "Kerr", "DFWM", "ESHG"]


def _table(table_id: str, tag: str, first: str) -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, table_id)
    if grid[1] != HEADER or not grid[0][0].startswith(first):
        msg = f"unexpected layout of {table_id}: {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"{PREFIX}{tag}_{name}.csv": [] for name in COLUMNS.values()
    }
    for cells in grid[2:]:
        abscissa = cells[0]
        float(abscissa)
        printed = " | ".join(cells)
        for column, name in COLUMNS.items():
            value = cells[column]
            float(value)
            files[f"{PREFIX}{tag}_{name}.csv"].append(
                Point(abscissa, abscissa, value, "", printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables II (by wavelength) and III (by omega), one file per column."""
    return _table(TABLE_III, "omega", "\\omega") | _table(
        TABLE_II, "lambda", "\\lambda"
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

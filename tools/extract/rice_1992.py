# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rice (1992): SCF and MP2 dynamic gamma of Ar, Kr and Ne.

Tables I (Ar, ``tbl_0346004ae107``), II (Kr, ``tbl_8eb3fba4bcde``) and III (Ne,
``tbl_0f061820910a``), cells files only; images checked, every value cell agrees. Each
table has rows alpha-bar, gamma^K(-w; w, 0, 0) and gamma_par(-2w; w, w, 0) for SCF and
for MP2, at omega = 0, 0.01, ..., 0.10 a.u. (the header also gives the wavelength in
angstrom). Only the gamma rows are taken (alpha-bar is linear data).

Structural fix (Table I): the library merged the two header rows (omega in a.u. and the
wavelength in angstrom) into one row of cells such as "0.0145 563" and "0∞". The
image shows the same header as Tables II and III (0, 0.01, ..., 0.10 over infinity,
45 563, ..., 4556): omega is split off the front of each merged cell and checked
against Table II's header.

The printed Kr SCF gamma^K values at 0.09 and 0.10 a.u. (2618, 2541) are swapped in the
table: corrected by errata in the record, not here.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Rice_1992_FrequencyDependentHyperpolarizabilities"
TABLES = {"Ar": "tbl_0346004ae107", "Kr": "tbl_8eb3fba4bcde", "Ne": "tbl_0f061820910a"}
ROWS = {
    "\\gamma^{K}(-\\omega;\\omega,0,0)": "dc_kerr",
    "\\gamma_{\\parallel}(-2\\omega;\\omega,\\omega,0)": "eshg",
}
_MERGED = re.compile(r"^(0(?:\.\d\d)?)(.+)$")


def _omegas(grid: list[list[str]], gas: str) -> list[str]:
    if gas != "Ar":
        return grid[1][1:]
    merged = [_MERGED.match(cell) for cell in grid[1][1:]]
    if any(m is None for m in merged):
        msg = f"unparsed merged header {grid[1]!r}"
        raise ValueError(msg)
    return [m.group(1) for m in merged if m is not None]


def extract() -> dict[str, list[Point]]:
    """Return gamma^K and gamma_par (ESHG), SCF and MP2, per gas."""
    reference = json_table(DIRECTORY, TABLES["Kr"])[1][1:]
    files: dict[str, list[Point]] = {}
    for gas, table_id in TABLES.items():
        grid = json_table(DIRECTORY, table_id)
        omegas = _omegas(grid, gas)
        if omegas != reference:
            msg = f"{gas}: frequencies {omegas!r} differ from Table II's"
            raise ValueError(msg)
        method = ""
        for cells in grid[2:]:
            if cells[0] in {"SCF", "MP2"}:
                method = cells[0].lower()
                continue
            process = ROWS.get(cells[0])
            if process is None:  # alpha-bar
                continue
            points = []
            for omega, value in zip(omegas, cells[1:], strict=True):
                float(value)
                points.append(Point(omega, omega, value, "", " | ".join(cells)))
            files[f"{gas}/rice_1992__{method}_{process}_{gas}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Beetar et al. (2020): third-order simulation inputs (Table S1).

The supplement's Table S1 is split by the library into two tables without CSV (merged
cells), read from their JSON cells (images checked):

- ``tbl_d1e82e48818a`` ("1. Density Matrix Model Simulations"): n2 of Ar, N2 and N2O
  printed as ``0.97 x Pressure`` (10^-19 cm^2/W per bar), the polarizability anisotropy
  of N2 and N2O printed as ``2 x 9.3`` (10^-25 cm^3; the footnote says Delta alpha was
  doubled to match experiment), rotational constants, dephasing times and spin weights.
- ``tbl_696db144fb76`` ("2. Capillary Propagation Simulations"): the "n2 ... for long
  pulse" of N2 and N2O (``3.2 x Pressure``), whose Ar entry repeats the first table's.

The printed multiplication sign is written x here. The ``x Pressure`` suffix and the
``2 x`` prefix are stripped; the values are per bar and the undoubled anisotropy. The
wavelength is the simulations' 1.03 um.
"""

from _json_table import json_table
from _library import Point, number, write

TIMES = "\N{MULTIPLICATION SIGN}"
DIRECTORY = "Beetar_2020_MultioctaveSupercontinuumGeneration/supplement_01"
DENSITY_MATRIX = "tbl_d1e82e48818a"
CAPILLARY = "tbl_696db144fb76"
GASES = {"Ar": "Ar", "N_2": "N2", "N_2O": "N2O", "N2": "N2", "N2O": "N2O"}


def _per_pressure(cell: str) -> str:
    """Return the number of a printed ``<v> x Pressure`` cell."""
    value, suffix = cell.split(TIMES)
    assert suffix.strip() == "Pressure", cell
    return number(value)


def _doubled(cell: str) -> str:
    """Return the anisotropy of a printed ``2 x <v>`` cell (the undoubled value)."""
    factor, value = cell.split(TIMES)
    assert factor.strip() == "2", cell
    return number(value)


def _rows(grid: list[list[str]]) -> dict[str, list[str]]:
    """Index a table's rows by their first cell."""
    return {row[0]: row for row in grid}


def extract() -> dict[str, list[Point]]:
    """Return electronic n2, Delta alpha and long-pulse n2 per gas."""
    first = json_table(DIRECTORY, DENSITY_MATRIX)
    rows = _rows(first)
    gases = [GASES[g] for g in rows[""][1:]]
    assert gases == ["Ar", "N2", "N2O"], rows[""]
    files: dict[str, list[Point]] = {}
    n2_row = rows["n_2(x 10-19cm2/W)"]
    for gas, cell in zip(gases, n2_row[1:], strict=True):
        files[f"{gas}/beetar_2020__n2_{gas}.csv"] = [
            Point("n2", "1.03", _per_pressure(cell), "", " | ".join(n2_row))
        ]
    alpha_row = rows["Δα (10-25 cm3)*"]
    for gas, cell in zip(gases[1:], alpha_row[2:], strict=True):
        files[f"{gas}/beetar_2020__delta_alpha_{gas}.csv"] = [
            Point("delta_alpha", "1.03", _doubled(cell), "", " | ".join(alpha_row))
        ]
    second = _rows(json_table(DIRECTORY, CAPILLARY))
    (long_key,) = [key for key in second if "for long pulse" in key]
    long_row = second[long_key]
    assert long_row[1] == n2_row[1], (long_row, n2_row)  # Ar: the same 0.97 x Pressure
    for gas, cell in zip(gases[1:], long_row[2:], strict=True):
        files[f"{gas}/beetar_2020__n2_long_{gas}.csv"] = [
            Point("n2_long", "1.03", _per_pressure(cell), "", " | ".join(long_row))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

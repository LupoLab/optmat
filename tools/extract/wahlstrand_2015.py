# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wahlstrand et al. (2015): Table I, H2 and D2 nonlinear coefficients (p. 5).

Columns: n2,elec^s (1e-20 cm2/W) this experiment and Ref. [35]; Delta alpha (1e-25 cm3)
this experiment and Ref. [30]; d alpha/dQ (1e-16 cm2, not extracted: no quantity); and
the long-pulse n2,rot, n2,vib, n2,tot at 800 nm (1e-20 cm2/W). Cells print
``value \pm error`` or a bare value. The library has no CSV for this table; its JSON
cells match the table image (checked). One single-row data file per column and gas.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Wahlstrand_2015_AbsoluteMeasurementUltrafast"
TABLE = "tbl_64c4c34f908d"
GASES = {"H_2": "H2", "D_2": "D2"}
COLUMNS = {
    1: "n2_elec",
    2: "n2_shelton_rice",
    3: "dalpha",
    4: "dalpha_bridge_buckingham",
    8: "n2_rot",
    9: "n2_vib",
    10: "n2_tot",
}


def _value_error(cell: str) -> tuple[str, str]:
    r"""Split ``value \pm error`` (or a bare value) into the printed numbers."""
    if r"\pm" in cell:
        value, error = cell.split(r"\pm")
        return number(value), number(error)
    return number(cell), ""


def extract() -> dict[str, list]:
    """Return one data file per extracted column and gas."""
    files: dict[str, list] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, name in COLUMNS.items():
            value, error = _value_error(cells[column])
            path = f"{gas}/wahlstrand_2015__{name}_{gas}.csv"
            # Bridge and Buckingham's anisotropy has no stated wavelength; the rest
            # refers to the 800-nm pump.
            wavelength = "not stated" if name == "dalpha_bridge_buckingham" else "800"
            files[path] = [Point(wavelength, wavelength, value, error, printed)]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

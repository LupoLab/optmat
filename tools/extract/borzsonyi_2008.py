# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Borzsonyi et al. (2008): specific GD, GDD and TOD at 800 nm (Table 1).

Table 1 (``tbl_15ed1378c327``, published by the library only as JSON cells because of
its two-row header) gives, per gas, the specific group delay (SGD, fs/m), group-delay
dispersion (SGDD, fs^2/m) and third-order dispersion (STOD, fs^3/m) at 800 nm "deducted
for 1000 mbar" (p. 5), each printed as ``value +- error``, the laboratory temperature
of the measurement, and Lorentz-Lorenz constants. Checked against the table image: the
cells are complete and the minus signs of the helium and neon STOD are hyphens in the
extraction and minus signs in print.

Air is not a curated species and is skipped. The Lorentz-Lorenz constants are not
extracted (see the source record). Each gas and quantity is its own dataset, because
the temperature (the dataset's state) differs from gas to gas.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Borzsonyi_2008_DispersionMeasurementInert"
TABLE = "tbl_15ed1378c327"

GASES = {
    "Nitrogen": "N2",
    "Helium": "He",
    "Neon": "Ne",
    "Argon": "Ar",
    "Krypton": "Kr",
    "Xenon": "Xe",
}
QUANTITIES = {1: "sgd", 2: "sgdd", 3: "stod"}


def _value_error(cell: str) -> tuple[str, str]:
    """Split ``"947.4 +- 1.0"`` (a printed plus-minus sign) into value and error."""
    value, error = cell.split("\N{PLUS-MINUS SIGN}")
    return number(value), number(error)


def extract() -> dict[str, list[Point]]:
    """Return one single-point file per gas and quantity."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        name = cells[0]
        if name not in GASES:  # Air
            continue
        gas = GASES[name]
        for column, quantity in QUANTITIES.items():
            value, error = _value_error(cells[column])
            point = Point(name, "800", value, error, " | ".join(cells))
            files[f"{gas}/borzsonyi_2008__{quantity}_{gas}.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

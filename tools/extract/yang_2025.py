# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Yang, Stone & Egan (2025): static molar polarizabilities of Ne, Ar, N2, He.

Table I (``tbl_ca4a321ccfce``, cells file only; image checked, cells agree) compares,
per gas, A_eps (cm^3/mol) and the first Cauchy coefficient A_2 (10^-4 cm^3/mol) of
Eq. (4), A_R(nu) = A_eps + A_2 nu_bar^2 + A_4 nu_bar^4 + A_6 nu_bar^6 + A_mu with
nu_bar = nu/(10^14 Hz), from theory, earlier experiment and this work. The gas is
named on the first row of its block; the reference ("Theoretical [20]", "Experimental
[17]", "This work") in the last column; "..." marks an empty cell. The parenthesised
last-digit standard uncertainty goes to the sigma column (none for N2's "Experimental
[22]" row, printed without one).

Only the A_eps column becomes data files (static molar polarizability, wavelength
"inf"); A_2 is not a quantity of the package and is recorded as formulas in the source
record.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Yang_2025_DemonstrationDispersionGas"
TABLE_I = "tbl_ca4a321ccfce"
GASES = {"Ne": "Ne", "Ar": "Ar", "N_{2}": "N2", "He": "He"}
# (gas, reference as printed) -> dataset name.
DATASETS = {
    ("Ne", "Theoretical [20]"): "lesiuk_2020",
    ("Ne", "Experimental [17]"): "gaiser_2018",
    ("Ne", "This work"): "A_eps",
    ("Ar", "Theoretical [21]"): "lesiuk_2023",
    ("Ar", "Experimental [17]"): "gaiser_2018",
    ("Ar", "This work"): "A_eps",
    ("N2", "Experimental [22]"): "kumar_meath_1992",
    ("N2", "Experimental [23]"): "schmidt_moldover_2003",
    ("N2", "This work"): "A_eps",
    ("He", "Theoretical [8,14,24]"): "puchalski",
    ("He", "Experimental [17]"): "gaiser_2018",
}
_VALUE = re.compile(r"^(?P<number>\d+\.\d+)(?:\((?P<digits>\d+)\))?$")


def _split(cell: str) -> tuple[str, str]:
    """Return a value and its last-digit uncertainty (or "") in the value's unit."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    value, digits = match["number"], match["digits"]
    if digits is None:
        return number(value), ""
    decimals = len(value.split(".")[1])
    return number(value), f"{int(digits) * 10.0**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and reference with a printed A_eps."""
    rows = json_table(DIRECTORY, TABLE_I)
    assert rows[0][1].startswith("A_{\\epsilon}"), rows[0]
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in rows[1:]:
        gas = GASES.get(cells[0].strip(), gas) if cells[0].strip() else gas
        if cells[1].strip() == "...":
            continue  # helium "This work": A_eps fixed to theory, only A_2 printed
        name = DATASETS[gas, cells[3].strip()]
        value, sigma = _split(cells[1])
        files[f"{gas}/yang_2025__{name}_{gas}.csv"] = [
            Point("static", "inf", value, sigma, " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kirouac & Bose (1973): dielectric virial coefficients of N2O.

Table I (``tbl_18fe8b762246``, CSV; image checked, cells agree) gives A_eps
(cm^3/mole), B_eps (cm^6/mole^2) and C_eps (cm^9/mole^3) at 6.5, 30.1 and 75.1 degrees
C, each as ``value +- uncertainty`` (the uncertainty printed without its leading zero,
``8.370 +- .006``; kept as printed). A_eps goes to one file with the temperature in K
(degrees C + 273.15) as the row id; each B_eps to its own file (one temperature per
dataset). C_eps has no quantity here and is left to the source notes.
"""

from _library import Point, number, table, write

DIRECTORY = "Kirouac_1973_DielectricPressureVirial"
TABLE_I = "tbl_18fe8b762246"
_PM = "\N{PLUS-MINUS SIGN}"


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``8.370 +- .006``, as printed."""
    value, sign, sigma = cell.partition(_PM)
    assert sign, cell
    return number(value), number(sigma)


def kelvin(celsius: str) -> str:
    """Return a printed Celsius temperature in kelvin, to the printed decimals."""
    return f"{float(number(celsius)) + 273.15:.2f}"


def extract() -> dict[str, list[Point]]:
    """Return A_eps (one file, three temperatures) and B_eps (one file each)."""
    rows = table(DIRECTORY, TABLE_I)
    assert rows[0][:3] == ["Temp (°C)", "Aε (cm3/mole)", "Bε (cm6/mole2)"], rows[0]
    files: dict[str, list[Point]] = {}
    a_eps = []
    for cells in rows[1:]:
        temperature = kelvin(cells[0])
        printed = " | ".join(cells)
        a_eps.append(Point(temperature, "inf", *split(cells[1]), printed))
        files[f"N2O/kirouac_1973__B_eps_{temperature}.csv"] = [
            Point(temperature, "inf", *split(cells[2]), printed)
        ]
    files["N2O/kirouac_1973__A_eps.csv"] = a_eps
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

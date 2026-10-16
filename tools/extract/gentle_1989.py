# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Gentle, Laver & Ritchie (1989): the dc Kerr effect of CO2 at 632.8 nm.

Both tables checked against their images (no differences in the cells).

* Table I (``tbl_18459ce921c2``, CSV): per temperature (9, 299.2-489.5 K) the number of
  pressures, the pressure range, the second density virial coefficient B used, A_K
  (10^27 m^5 V^-2 mol^-1) and B_K (10^32; the header prints m^8 V^-2 mol^-1, see the
  record), each "value +- error". Taken: A_K (``table1_A_K_CO2``; the printed errors
  are fit errors and stay in ``printed`` only, the dataset carries the stated 2 %
  accuracy) and B_K (``table1_B_K_CO2``, with its printed error as sigma). Temperature
  as printed is the row id and the abscissa.
* Table II (``tbl_1b1d9fc6dd0c``, CSV): gamma^K (``table2_gamma_K_CO2``), the static
  anisotropy Delta alpha^0 (``table2_anisotropy_static_CO2``) and the optical Delta
  alpha of ref 9 that the analysis uses (``table2_anisotropy_optical_CO2``). The
  intercept and slope are the record's formula; Delta alpha Delta alpha^0 is in the
  notes. Footnote markers "^a", "^b" are dropped.
"""

from _library import Point, number, table, write

DIRECTORY = "Gentle_1989_SecondHyperpolarizabilityStatic"
TABLE_1 = "tbl_18459ce921c2"
TABLE_2 = "tbl_1b1d9fc6dd0c"
WAVELENGTH = "632.8"  # nm, the He-Ne line (Fig. 1 caption, Table I caption)
PLUS_MINUS = "\N{PLUS-MINUS SIGN}"
TABLE_2_ROWS = {  # start of the property cell -> (dataset, wavelength)
    "10^{60}\\gamma^K": ("gamma_K_CO2", WAVELENGTH),
    "10^{40}\\Delta\\alpha/": ("anisotropy_optical_CO2", WAVELENGTH),
    "10^{40}\\Delta\\alpha^0": ("anisotropy_static_CO2", "inf"),
}


def _split(cell: str, separator: str) -> tuple[str, str]:
    """Split "1.508 +- 0.004" into its value and error, as printed."""
    value, error = cell.split(separator)
    return number(value), number(error.strip().removesuffix("^a").removesuffix("^b"))


def _table1() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_1)
    header = grid[0]
    a_k, b_k = header[5], header[6]
    if not (a_k.startswith("10^{27} A_K") and b_k.startswith("10^{32} B_K")):
        msg = f"unexpected Table I header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        temperature = number(cells[0])
        printed = " | ".join(cells)
        value, _ = _split(cells[5], PLUS_MINUS)  # the fit error stays in printed
        files.setdefault("CO2/gentle_1989__table1_A_K_CO2.csv", []).append(
            Point(temperature, temperature, value, "", printed)
        )
        value, error = _split(cells[6], PLUS_MINUS)
        files.setdefault("CO2/gentle_1989__table1_B_K_CO2.csv", []).append(
            Point(temperature, temperature, value, error, printed)
        )
    return files


def _table2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        for start, (name, wavelength) in TABLE_2_ROWS.items():
            if cells[0].startswith(start):
                value, error = _split(cells[1], "\\pm")
                files[f"CO2/gentle_1989__table2_{name}.csv"] = [
                    Point(wavelength, wavelength, value, error, " | ".join(cells))
                ]
    if len(files) != len(TABLE_2_ROWS):
        msg = "Table II: a property row is missing"
        raise ValueError(msg)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table I and Table II files."""
    return _table1() | _table2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

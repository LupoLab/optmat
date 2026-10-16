# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton & Buckingham (1982): ESHG hyperpolarizabilities, Tables II and III.

Table II (``tbl_09cd9a18fa27``, CSV; image checked) gives 10^63 gamma (C^4 m^4 J^-3) at
514.5 nm; Table III (``tbl_e6603d44b34e``, CSV; image checked) gives gamma_He
(calculated, the calibration) and gamma_CH4 at seven wavelengths, the 694.3 nm CH4 value
marked e) "From Ref. 4" (Miller & Ward 1977). Files: ``table2_<gas>`` for Ar, H2, N2,
CO2, SF6 (the CH4 row of Table II repeats Table III's 514.5 nm row and is taken from
Table III only; He's 2.950 is the calibrant, in ``table3_He_calc``; CH3F, CHF3, CF4 are
not curated), ``table3_He_calc``, ``table3_CH4`` (the six argon-ion lines) and
``table3_CH4_miller_ward`` (694.3 nm). Row ids are the printed air wavelengths. The
phase-matching densities (linear optics) are not extracted (see the record).
"""

from _library import Point, number, table, write

DIRECTORY = "Shelton_1982_OpticalSecondHarmonic"
TABLE2 = "tbl_09cd9a18fa27"
TABLE3 = "tbl_e6603d44b34e"
TABLE2_GASES = {"Ar": "Ar", "H_2": "H2", "N_2": "N2", "CO_2": "CO2", "SF_6": "SF6"}


def _split(cell: str) -> tuple[str, str]:
    cleaned = cell.replace("\\pm", "\N{PLUS-MINUS SIGN}").removesuffix("^e")
    value, error = cleaned.split("\N{PLUS-MINUS SIGN}")
    return number(value), number(error)


def extract() -> dict[str, list[Point]]:
    """Return the Table II gammas and the Table III dispersion of He and CH4."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE2)[1:]:
        if cells[0] not in TABLE2_GASES:
            continue
        gas = TABLE2_GASES[cells[0]]
        files[f"{gas}/shelton_1982__table2_{gas}.csv"] = [
            Point("514.5", "514.5", *_split(cells[1]), " | ".join(cells))
        ]
    he, ch4, miller_ward = [], [], []
    for cells in table(DIRECTORY, TABLE3)[1:]:
        printed = " | ".join(cells)
        he.append(Point(cells[0], cells[0], number(cells[1]), "", printed))
        point = Point(cells[0], cells[0], *_split(cells[3]), printed)
        (miller_ward if cells[3].endswith("^e") else ch4).append(point)
    files["He/shelton_1982__table3_He_calc.csv"] = he
    files["CH4/shelton_1982__table3_CH4.csv"] = ch4
    files["CH4/shelton_1982__table3_CH4_miller_ward.csv"] = miller_ward
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Carusotto, Perrone & Polacco (1992): dc Kerr constants of N2, H2, D2.

Table I (``tbl_bd417c3b63d6``, CSV; image checked) gives Kerr constants K (Delta n = K
rho E^2, rho in amagat) of N2, H2 and D2, in the column headed "10^17 cm^2/statV^2
amagat" (the text writes x 10^-17; the values are 10^-17 esu per amagat): the present
work at 295 K (value +- statistical +- systematic) and earlier values. The present rows
go to ``kerr_<gas>`` (sigma: the first, statistical, uncertainty as printed in the
table; the systematic one is in the record); the rows of Ref. 3 (Carusotto et al.
1985, at 514.5 nm) to ``kerr_<gas>_ref3``. The rows of Refs. 5 and 2 (high pressure,
temperature ranges) are recorded in the notes only. Table II (``tbl_6fc9aecd71e0``,
CSV; image checked): the H2 and D2 values of gamma (esu x 10^-37), Delta alpha and
Delta alpha_0 (cm^3 x 10^-25) used in Eq. (12), one file per gas and quantity.
"""

from _library import Point, number, table, write

DIRECTORY = "Carusotto_1992_MeasurementKerrConstant"
TABLE1 = "tbl_bd417c3b63d6"
TABLE2 = "tbl_6fc9aecd71e0"
GASES = {"N_2": "N2", "H_2": "H2", "D_2": "D2"}
TABLE2_ROWS = {
    "\\gamma \\times 10^{-37} (esu)": "gamma",
    "\\Delta\\alpha \\times 10^{-25} (cm^{3})": "anisotropy_optical",
    "\\Delta\\alpha_{0} \\times 10^{-25} (cm^{3})": "anisotropy_static",
}


def extract() -> dict[str, list[Point]]:
    """Return the present and Ref. 3 Kerr constants and the Table II values."""
    rows = table(DIRECTORY, TABLE1)[1:]
    files: dict[str, list[Point]] = {}
    for index, cells in enumerate(rows):
        if index >= 6:  # rows of Refs. 5 and 2: in the record's notes
            break
        gas = GASES[cells[0]]
        suffix = "" if index < 3 else "_ref3"
        parts = cells[3].split("\N{PLUS-MINUS SIGN}")
        wavelength = "632.8" if index < 3 else "514.5"
        files[f"{gas}/carusotto_1992__kerr_{gas}{suffix}.csv"] = [
            Point(
                f"{cells[1]} atm, {cells[2]} K",
                wavelength,
                number(parts[0]),
                number(parts[1]),
                " | ".join(cells),
            )
        ]
    for cells in table(DIRECTORY, TABLE2)[1:]:
        name = TABLE2_ROWS.get(cells[0])
        if name is None:
            continue  # B_0 and eta: in the record's notes
        for column, gas in ((1, "H2"), (2, "D2")):
            wavelength = "inf" if name == "anisotropy_static" else "632.8"
            files[f"{gas}/carusotto_1992__{name}_{gas}.csv"] = [
                Point(name, wavelength, number(cells[column]), "", " | ".join(cells))
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

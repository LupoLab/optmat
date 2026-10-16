# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1986): the CH4/N2 polarizability-dispersion ratio by ESHG.

Table I for CH4 (``tbl_9090f1c5cfc2``, image checked) gives the fundamental's air
wavelength (nm), its vacuum wavenumber (cm^-1), the ratio
``Delta alpha_CH4 / Delta alpha_N2`` (Delta alpha = alpha(2 omega) - alpha(omega),
equal to the phase-match density ratio rho_N2 / rho_CH4) with its +-, and
hyperpolarizability values (not extracted), 476.5-671.7 nm. The row id is the printed
air wavelength; the abscissa is the printed vacuum wavenumber. The 694.3 and 457.9 nm
rows (gamma values of other authors only) have no ratio and are skipped.

Third-order data: Table I's gamma_CH4/gamma_N2 (``gamma_CH4_over_N2``) and
gamma_CH4 (``gamma_CH4``, column d, 1e-63 C^4 m^4 J^-3); the last column's values of
other authors, marked e) (Miller & Ward, 694.3 nm: ``gamma_CH4_miller_ward``) and f)
(Shelton & Buckingham 1982: ``gamma_CH4_shelton_1982``), the marks stripped. Table II
(``tbl_0ba11cfe8ccf``, CSV; image checked): gamma_CH4 by other processes, one file per
process (abscissa: the printed wavelength, for FWM the pump of "694,975" and "532,683").
"""

from _library import Point, number, table, write

DIRECTORY = "Shelton_1986_HyperpolarizabilityDispersionMeasured_37a53a"
TABLE = "tbl_9090f1c5cfc2"
TABLE2 = "tbl_0ba11cfe8ccf"
TABLE2_ROWS = {
    "dc Kerr^{c}": ("dc_kerr", "ref 21 Buckingham and Orr"),
    "FWM^{d}": ("cars", "ref 22 Hauchecorne et al."),
    "FWM^{e}": ("cars", "ref 23 Rado as corrected in ref 24, vibration restored"),
    "FWM^{f}": ("cars", "ref 24 Lundeen et al., vibration restored"),
    "THG^{g}": ("thg", "ref 25 Ward and Elliott"),
}


def extract() -> dict[str, list[Point]]:
    """Return the dispersion ratios and the gamma data as data files."""
    points = []
    for cells in table(DIRECTORY, TABLE)[1:]:
        if not cells[2].strip():
            continue
        value, error = cells[2].split("\N{PLUS-MINUS SIGN}")
        points.append(
            Point(
                number(cells[0]),
                number(cells[1]),
                number(value),
                number(error),
                " | ".join(cells),
            )
        )
    files = {"N2/shelton_1986_ch4__density_N2_over_CH4.csv": points}
    stem = "CH4/shelton_1986_ch4__"
    for cells in table(DIRECTORY, TABLE)[1:]:
        printed = " | ".join(cells)
        label, wavenumber = number(cells[0]), number(cells[1])
        for column, name in ((3, "gamma_CH4_over_N2"), (4, "gamma_CH4")):
            if cells[column]:
                value, error = cells[column].split("\N{PLUS-MINUS SIGN}")
                files.setdefault(f"{stem}{name}.csv", []).append(
                    Point(label, wavenumber, number(value), number(error), printed)
                )
        if cells[5]:
            other = cells[5]
            name = "miller_ward" if other.endswith("^e") else "shelton_1982"
            value, error = other[:-2].split("\N{PLUS-MINUS SIGN}")
            files.setdefault(f"{stem}gamma_CH4_{name}.csv", []).append(
                Point(label, wavenumber, number(value), number(error), printed)
            )
    for cells in table(DIRECTORY, TABLE2)[1:]:
        process, row = TABLE2_ROWS[cells[0]]
        value, error = cells[3].split("\N{PLUS-MINUS SIGN}")
        wavelength = number(cells[1].split(",")[0])
        files.setdefault(f"{stem}table2_{process}.csv", []).append(
            Point(row, wavelength, number(value), number(error), " | ".join(cells))
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

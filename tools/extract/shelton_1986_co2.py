# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1986): the CO2/N2 polarizability-dispersion ratio by ESHG.

Table I for CO2 (``tbl_ba9d1caa914e``, image checked) gives the fundamental's air
wavelength (nm), its vacuum wavenumber (cm^-1), the ratio
``Delta alpha_CO2 / Delta alpha_N2`` (Delta alpha = alpha(2 omega) - alpha(omega),
equal to the phase-match density ratio rho_N2 / rho_CO2) with its +-, and
hyperpolarizability values (not extracted), 488.0-692.9 nm. The row id is the printed
air wavelength; the abscissa is the printed vacuum wavenumber. Every row has a
ratio.

Third-order data: Table I's gamma_CO2/gamma_N2 (``gamma_CO2_over_N2``) and
gamma_CO2 in 1e-63 C^4 m^4 J^-3 (``gamma_CO2``), same rows and abscissa. Table II
(``tbl_c8ef15b5d66e``, CSV; image checked) collects gamma_CO2 from other processes and
authors, one file per process (row id: the footnote letter's source; abscissa: the
printed wavelength in nm, for FWM the first, pump, wavelength of "694, 975" and "532,
683").
"""

from _library import Point, number, table, write

DIRECTORY = "Shelton_1986_HyperpolarizabilityDispersionMeasured"
TABLE = "tbl_ba9d1caa914e"
TABLE2 = "tbl_c8ef15b5d66e"
# Table II rows: (process label as printed, file, row id)
TABLE2_ROWS = {
    "dc Kerr^{\\mathrm{c}}": ("dc_kerr", "ref 21 Buckingham et al. 1970"),
    "FWM^{\\mathrm{d}}": ("cars", "ref 23 Rado as corrected in ref 24"),
    "FWM^{\\mathrm{e}}": ("cars", "ref 24 Lundeen et al."),
    "ESHG^{\\mathrm{f}}": ("eshg", "ref 22 Ward and Miller"),
    "ESHG^{\\mathrm{g}}": ("eshg", "ref 13 Shelton and Buckingham"),
    "THG^{\\mathrm{h}}": ("thg", "ref 25 Ward and New"),
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
    files = {"N2/shelton_1986_co2__density_N2_over_CO2.csv": points}
    ratios, gammas = [], []
    for cells in table(DIRECTORY, TABLE)[1:]:
        printed = " | ".join(cells)
        for column, target in ((3, ratios), (4, gammas)):
            value, error = cells[column].split("\N{PLUS-MINUS SIGN}")
            target.append(
                Point(
                    number(cells[0]),
                    number(cells[1]),
                    number(value),
                    number(error),
                    printed,
                )
            )
    files["CO2/shelton_1986_co2__gamma_CO2_over_N2.csv"] = ratios
    files["CO2/shelton_1986_co2__gamma_CO2.csv"] = gammas
    for cells in table(DIRECTORY, TABLE2)[1:]:
        process, row = TABLE2_ROWS[cells[0]]
        value, error = cells[3].split("\N{PLUS-MINUS SIGN}")
        wavelength = number(cells[1].split(",")[0])
        path = f"CO2/shelton_1986_co2__table2_{process}.csv"
        files.setdefault(path, []).append(
            Point(row, wavelength, number(value), number(error), " | ".join(cells))
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

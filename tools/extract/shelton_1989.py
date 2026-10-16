# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Shelton (1989): the Ne/He linear-dispersion ratios of Table I.

Table I (``tbl_96462ca3eb80``, cells file only; image checked) lists nonlinear-optical
measurements for Ne; the ESHG rows also give Delta alpha_Ne / Delta alpha_He with
Delta alpha = alpha(2 nu) - alpha(nu), "determined from phase-match density
measurements" (the Ref. 18 row "calculated from refractive index data"). To first
order in density that is rho_He / rho_Ne at equal phase mismatch, the package's
``phase_match_density_ratio`` of He with reference gas Ne, so the printed number is
kept as is (as in ``shelton_1990``). One file per source of the rows: the present
work ("a"), Shelton and Lu (Ref. 15), Finn and Ward (Ref. 18). Rows and abscissa are
the vacuum wavenumbers nu (cm^-1): the lambda column of the last Ref. 15 row prints
488.8 nm where nu = 20487 cm^-1 (and nu_L^2 = 6 nu^2) is the 488.0 nm argon line; the
record notes it. Rows whose ratio is "..." (Kerr, THG) are skipped; the gamma ratios
are not refractive data.

Third-order data: the gamma_Ne/gamma_He column, one file per source
(``gamma_Ne_over_He_<source>``: the dc Kerr rows of Ref. 17, whose two numbers divide
the same gamma_Ne by the measured and by the ab initio gamma_He, rows "15798 He
measured" and "15798 He ab initio"; the present ESHG rows; Refs. 18 and 15; the THG rows
of Refs. 19 and 20). Row ids are the vacuum wavenumbers, as above.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Shelton_1989_AnomalousHyperpolarizabilityDispersion"
TABLE = "tbl_96462ca3eb80"
SOURCES = {"a": "present", "15": "shelton_lu_1988", "18": "finn_ward_1971"}
GAMMA_SOURCES = {
    ("Kerr", "17"): "buckingham_dunmur_kerr",
    ("ESHG", "a"): "present",
    ("ESHG", "18"): "finn_ward_1971",
    ("ESHG", "15"): "shelton_lu_1988",
    ("THG", "19"): "lehmeier_thg",
    ("THG", "20"): "ward_new_thg",
}


def extract() -> dict[str, list[Point]]:
    """Return the Delta alpha and gamma ratios, one file per source."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        ratio = cells[5].strip()
        if ratio == "...":
            continue
        value, error = ratio.split("\N{PLUS-MINUS SIGN}")
        wavenumber = cells[3].strip()
        name = SOURCES[cells[1].strip()]
        path = f"He/shelton_1989__{name}_He_over_Ne.csv"
        files.setdefault(path, []).append(
            Point(
                wavenumber,
                number(wavenumber),
                number(value),
                number(error),
                " | ".join(cells),
            )
        )
    process = ""
    kerr_rows = iter(["He measured", "He ab initio"])
    for cells in json_table(DIRECTORY, TABLE)[1:]:
        process = cells[0].strip() or process
        wavenumber = cells[3].strip()
        name = GAMMA_SOURCES[process, cells[1].strip()]
        row = wavenumber if process != "Kerr" else f"{wavenumber} {next(kerr_rows)}"
        value, error = cells[6].split("\N{PLUS-MINUS SIGN}")
        path = f"Ne/shelton_1989__gamma_Ne_over_He_{name}.csv"
        files.setdefault(path, []).append(
            Point(
                row, number(wavenumber), number(value), number(error), " | ".join(cells)
            )
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

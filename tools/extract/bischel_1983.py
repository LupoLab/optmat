# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bischel and Black (1983): the single-resonance fits of their Table 1.

Table 1 (the library's unlabelled ``tbl_377daa8e19f8``, p. 7; its caption, 'Single
resonance fit to Raman scattering cross sections', and the fitted form are printed above
it; the CSV cells match the image) gives per gas the Raman shift (cm^-1), A (10^-28 cm^2
sr^-1), nu_i (10^4 cm^-1) and lambda_i (nm). The shifts are written as data (their CH4
misprint is an erratum in the record). A, nu_i and lambda_i are the coefficients of the
record's formulas; ``extract`` checks that the recorded strings are the printed cells,
so the formulas cannot drift from the table.
"""

import tomllib

from _library import DATA, Point, number, table, write

DIRECTORY = "Bischel_1983_WavelengthDependenceRaman"
TABLE_1 = "tbl_377daa8e19f8"
HEADER = (
    "Gas",
    "Raman Shift(cm^{-1})",
    "A(x10^{-28})cm^{2} sr^{-1}",
    "\\nu_{1}(x10^{4} cm^{-1})",
    "\\lambda_{1}(nm)",
)
GASES = {
    "H_{2}": "H2",
    "D_{2}": "D2",
    "N_{2}": "N2",
    "O_{2}": "O2",
    "CH_{4} (v_{1})": "CH4",
}
TRANSITIONS = {
    "H2": "v=0-1 Q branch",
    "D2": "v=0-1 Q branch",
    "N2": "v=0-1 Q branch",
    "O2": "v=0-1 Q branch",
    "CH4": "nu1 Q branch",
}


def _formulas() -> dict[str, dict[str, str]]:
    """Return the coefficients of the record's formulas, by material."""
    with (DATA / "sources" / "bischel_1983.toml").open("rb") as handle:
        record = tomllib.load(handle)
    return {entry["material"]: entry["coefficients"] for entry in record["formulas"]}


def extract() -> dict[str, list[Point]]:
    """Return the Raman shifts, after checking the formulas against the table."""
    rows = table(DIRECTORY, TABLE_1)
    if tuple(rows[0]) != HEADER:
        msg = f"Table 1 header {rows[0]!r}"
        raise ValueError(msg)
    formulas = _formulas()
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        # The formula's coefficients are the printed cells times the column factors.
        expected = {
            "A": f"{number(cells[2])}e-28",
            "nu_i": f"{number(cells[3])}e4",
            "lambda_i": number(cells[4]),
        }
        recorded = {key: formulas[gas][key] for key in expected}
        if recorded != expected:
            msg = f"{gas}: the record's coefficients {recorded} are not {expected}"
            raise ValueError(msg)
        label = TRANSITIONS[gas]
        files[f"{gas}/bischel_1983__shift_{gas}.csv"] = [
            Point(label, label, number(cells[1]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

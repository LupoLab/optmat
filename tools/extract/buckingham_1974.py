# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Buckingham & Graham (1974): the first and second refractivity virials.

Two tables (cells files only; images checked, cells agree):

- ``tbl_0cfd7fe493cd`` holds Table 4 (rows 0-9, representative single runs, not
  extracted) and Table 5 (from row 10): per gas, 10^6 A_R (m^3/mol = cm^3/mol), the
  Landolt-Bornstein A_R (not extracted: no wavelength), the mean temperature, 10^12 B_R
  (m^6/mol^2 = cm^6/mol^2; the header's "m^6 mol^-3" is a misprint), its standard
  deviation, the drift term B_Delta, the "estimated uncertainty" (their sum with the
  distortion term, taken as sigma), the calculated classical B_F and the literature
  B_R. CHF3 is not curated. Wavelength: the He-Ne line, lambda_0 = 633.0 nm in vacuo
  (p. 280).
- Table 1 (``tbl_8d9922d080a7``), the literature compilation: only Bose and Cole's
  (1970) dielectric A_eps and B_eps of argon (row "Bose & Cole 1970; Sutter 1972", in
  parentheses, infinite wavelength) are extracted; the other curated rows at a single
  wavelength or at infinity are orcutt_1967's values, identical and already curated,
  and the rows "447.1 to 667.8" nm have no single wavelength.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Buckingham_1974_DensityDependenceRefractivity"
TABLES_4_5 = "tbl_0cfd7fe493cd"
TABLE_1 = "tbl_8d9922d080a7"
WAVELENGTH = "633.0"  # nm, "lambda_0 = 633.0 nm" in vacuo (p. 280)
GASES = {
    "Ne": "Ne",
    "Ar": "Ar",
    "N_{2}": "N2",
    "CO_{2}": "CO2",
    "CH_{4}": "CH4",
    "SF_{6}": "SF6",
}
_PARENTHESISED = re.compile(r"^\((?P<value>-?[\d.]+) ± (?P<sigma>[\d.]+)\)$")


def _path(gas: str, name: str) -> str:
    return f"{gas}/buckingham_1974__{name}_{gas}.csv"


def _table_5() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLES_4_5)
    start = next(i for i, r in enumerate(rows) if r[0].startswith("TABLE 5"))
    assert rows[start + 1][0] == "gas", rows[start + 1]
    files: dict[str, list[Point]] = {}
    for cells in rows[start + 2 :]:
        gas = GASES.get(cells[0].strip())
        if gas is None:
            continue  # CHF3
        printed = " | ".join(cells)
        files[_path(gas, "A_R")] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[1]), "", printed)
        ]
        files[_path(gas, "B_R")] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[4]), number(cells[7]), printed)
        ]
        files[_path(gas, "B_F_theory")] = [
            Point(WAVELENGTH, WAVELENGTH, number(cells[8]), "", printed)
        ]
    return files


def _bose_cole() -> dict[str, list[Point]]:
    for cells in json_table(DIRECTORY, TABLE_1):
        if cells[5].strip().startswith("Bose & Cole 1970"):
            printed = " | ".join(cells)
            files = {}
            names = ((3, "A_eps_bose_cole_1970"), (4, "B_eps_bose_cole_1970"))
            for column, name in names:
                match = _PARENTHESISED.match(cells[column].strip())
                assert match is not None, cells[column]
                files[_path("Ar", name)] = [
                    Point("inf", "inf", match["value"], match["sigma"], printed)
                ]
            return files
    msg = "Bose & Cole row not found"
    raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return Table 5's A_R, B_R and B_F, and Bose and Cole's argon values."""
    return _table_5() | _bose_cole()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

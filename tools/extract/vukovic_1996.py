# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Vukovic, Woolsey and Scelsi (1996): SF6 refractivity at 632.99 and 1300 nm.

Tables 1 (``tbl_dc762adbe049``, 632.99 nm) and 2 (``tbl_786fca828745``, 1300 nm) give,
for SF6 and SOF2 (not curated), the refractivity x 10^6 at 101 325 Pa and the
measurement temperature, the same converted to 25 C, the first refractivity virial
coefficient A_R (10^-6 m^3/mol, i.e. cm^3/mol) and the coefficients A and B of
n - 1 = A (p/RT) + B (p/RT)^2. The cells run the value, its +- and the temperature
together (``708.05 ± 0.5at 22.9 °C``). Images checked.

Table 3 (``tbl_d5adcde5ec76``, image checked) compares SF6 at 632.99 nm with earlier
work: Buckingham and Graham's and St-Arnaud and Bose's A_R and B_R are extracted as
reproduced values. Its A_R header prints ``10^-5 m^3 mol^-1`` but the values are in
10^-6 m^3/mol, as in Table 1 (``11.336`` in both): the unit is the dataset's. The
refractivity and A 'calculated' by the authors from others' coefficients, and Watson
and Ramaswamy's value (curated as ``watson_1936``) are not extracted.
"""

import re

from _library import Point, table, write

DIRECTORY = "Vukovic_1996_RefractivitiesWavelengths632"
TABLES = {"632.99": "tbl_dc762adbe049", "1300": "tbl_786fca828745"}
COMPARISON = "tbl_d5adcde5ec76"
_VALUE = re.compile(r"^\s*([\d.]+)\s*(?:±\s*([\d.]+))?")
PREFIX = "SF6/vukovic_1996__"


def _value(cell: str) -> tuple[str, str]:
    """Return the printed value and its +- (or "") at the start of a cell."""
    match = _VALUE.match(cell)
    if match is None:
        msg = f"no value in {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2) or ""


def _this_work() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for wavelength, table_id in TABLES.items():
        (cells,) = [r for r in table(DIRECTORY, table_id) if r[0] == "SF_{6}"]
        printed = " | ".join(cells)
        measured, measured_sigma = _value(cells[1])
        files[f"{PREFIX}measured_{wavelength.split('.')[0]}.csv"] = [
            Point(wavelength, wavelength, measured, measured_sigma, printed)
        ]
        at_25, sigma_25 = _value(cells[2])
        files.setdefault(f"{PREFIX}refractivity_25C.csv", []).append(
            Point(wavelength, wavelength, at_25, sigma_25, printed)
        )
        a_r, a_r_sigma = _value(cells[3])
        files.setdefault(f"{PREFIX}A_R.csv", []).append(
            Point(wavelength, wavelength, a_r, a_r_sigma, printed)
        )
    return files


def _comparison() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = {r[0]: r for r in table(DIRECTORY, COMPARISON)[1:]}
    for author, dataset in (
        ("Buckingham and Graham (1974)", "buckingham_graham"),
        ("St-Arnaud and Bose (1979)", "st_arnaud_bose"),
    ):
        cells = rows[author]
        printed = " | ".join(cells)
        a_r, a_r_sigma = _value(cells[3])
        files[f"{PREFIX}{dataset}_A_R.csv"] = [
            Point("632.99", "632.99", a_r, a_r_sigma, printed)
        ]
        b_r, b_r_sigma = _value(cells[4])
        files[f"{PREFIX}{dataset}_B_R.csv"] = [
            Point("632.99", "632.99", b_r, b_r_sigma, printed)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return SF6's values of Tables 1, 2 and 3 as data files."""
    return _this_work() | _comparison()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

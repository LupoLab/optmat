# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hohm & Kerl (1990), Part II: alpha(omega, T) of Ar, H2 and O2.

All five tables have CSV renderings; each was checked against its image (no
differences; decimal points printed as raised dots, turned into points here).

* Tables 1-3 (Ar ``tbl_8e74d3a3b9bc``, H2 ``tbl_5e447ea27a86``, O2
  ``tbl_c07dae7af5ae``): per temperature the single-term fit parameters of Eq. (18)
  (recorded as formulas in the source record) and the static polarizability
  alpha(0, T) = C0/C1 (10^30 alpha/m^3, angstrom^3). The six rows with a temperature
  uncertainty in parentheses ("295.60 (22)") are this work (``alpha0_T_<gas>``, row =
  the temperature in K as printed, without its uncertainty); the rows marked [8], [16]
  and [17] are the authors' re-evaluation of other workers' dispersion formulas
  (one file each, row = its temperature); the rows with a blank temperature carry
  only a static polarizability from another source, marked "[n]" (one file each, row
  "ref"; the parenthesised digits of "1.64110 (24) [19]" are its uncertainty).
* Tables 4 and 5 (H2 ``tbl_54ded087f133``, O2 ``tbl_a1ae6c41a18b``): the dynamic
  polarizabilities extrapolated to T = 0 K (this work, with the parenthesised
  uncertainty on the last digits), the ab initio values they are compared with, and
  the derivative alpha'_e = r_e (d alpha/d r)_e of this work and of the ab initio
  work. Rows are the wavelengths as printed ("∞" -> "inf"); the ratio column a1/a2 is
  not extracted (the record's notes).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Hohm_1990_InterferometricMeasurementsDipole_d234f3"
DOT = "\N{MIDDLE DOT}"
INFINITY = "\N{INFINITY}"
FIT_TABLES = {
    "Ar": "tbl_8e74d3a3b9bc",
    "H2": "tbl_5e447ea27a86",
    "O2": "tbl_c07dae7af5ae",
}
ZERO_K_TABLES = {"H2": "tbl_54ded087f133", "O2": "tbl_a1ae6c41a18b"}
# The file suffix of each literature row: its reference number -> the cited work.
REFERENCES = {
    "8": "kerl_1982",
    "16": "cuthbertson_1910",
    "17": "cuthbertson_1913",
    "18": "kumar_1985",
    "19": "bose_1970",
    "20": "orcutt_1967",
    "21": "zeiss_1977",
    "22": "bishop_1980",
}
ZERO_K_COLUMNS = {  # column -> (file stem, quantity kind)
    "H2": {1: "alpha_0K", 2: "bishop_1980_0K", 3: "dalpha", 4: "adamov_1965_dalpha"},
    "O2": {1: "alpha_0K", 2: "rijks_1989_0K", 3: "dalpha", 4: "inoue_1985_dalpha"},
}
_THIS_WORK = re.compile(r"^([\d.]+) \((\d+)\)$")
_CITED_T = re.compile(r"^([\d.]+) \[(\d+)\]$")
_CITED_VALUE = re.compile(r"^([\d.]+)(?: \((\d+)\))? \[(\d+)\]$")
_UNCERTAIN = re.compile(r"^([\d.]+) \((\d+)\)$")


def _clean(cell: str) -> str:
    return cell.replace(DOT, ".").strip()


def _last_digits(value: str, digits: str) -> str:
    """Return the uncertainty "(57)" on the last digits of "0.80320" as "0.00057"."""
    decimals = len(value.partition(".")[2])
    return f"{int(digits) / 10**decimals:.{decimals}f}"


def _fits(gas: str) -> dict[str, list[Point]]:
    grid = table(DIRECTORY, FIT_TABLES[gas])
    if not (grid[0][0] == "T/K" and "alpha(0" in grid[0][5]):
        msg = f"unexpected header in {FIT_TABLES[gas]}: {grid[0]!r}"
        raise ValueError(msg)
    prefix = f"{gas}/hohm_1990b__"
    files: dict[str, list[Point]] = {}
    for raw in grid[1:]:
        cells = [_clean(cell) for cell in raw]
        printed = " | ".join(raw)
        own = _THIS_WORK.match(cells[0])
        cited = _CITED_T.match(cells[0])
        if own:
            temperature = number(own[1])
            files.setdefault(f"{prefix}alpha0_T_{gas}.csv", []).append(
                Point(temperature, "inf", number(cells[5]), "", printed)
            )
        elif cited:
            temperature, reference = number(cited[1]), cited[2]
            name = f"{prefix}{REFERENCES[reference]}_{gas}.csv"
            files[name] = [Point(temperature, "inf", number(cells[5]), "", printed)]
        elif cells[0] == "" and cells[5]:
            match = _CITED_VALUE.match(cells[5])
            if match is None:
                msg = f"unexpected literature cell {cells[5]!r}"
                raise ValueError(msg)
            value, digits, reference = match[1], match[2], match[3]
            sigma = _last_digits(value, digits) if digits else ""
            name = f"{prefix}{REFERENCES[reference]}_{gas}.csv"
            files[name] = [Point("ref", "inf", number(value), sigma, printed)]
        elif cells[0] == "" and cells[6]:
            continue  # a C6 value only (no quantity; the record's notes)
        else:
            msg = f"unexpected row {raw!r}"
            raise ValueError(msg)
    return files


def _zero_kelvin(gas: str) -> dict[str, list[Point]]:
    grid = table(DIRECTORY, ZERO_K_TABLES[gas])
    if grid[0][0] != "λ/nm" or "a_{1}/a_{2}" not in grid[0][5]:
        msg = f"unexpected header in {ZERO_K_TABLES[gas]}: {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for raw in grid[1:]:
        cells = [_clean(cell) for cell in raw]
        wavelength = "inf" if cells[0] == INFINITY else number(cells[0])
        for column, stem in ZERO_K_COLUMNS[gas].items():
            cell = cells[column]
            match = _UNCERTAIN.match(cell)
            value, sigma = (
                (match[1], _last_digits(match[1], match[2])) if match else (cell, "")
            )
            files.setdefault(f"{gas}/hohm_1990b__{stem}_{gas}.csv", []).append(
                Point(wavelength, wavelength, number(value), sigma, " | ".join(raw))
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return Tables 1-5."""
    files: dict[str, list[Point]] = {}
    for gas in FIT_TABLES:
        files |= _fits(gas)
    for gas in ZERO_K_TABLES:
        files |= _zero_kelvin(gas)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

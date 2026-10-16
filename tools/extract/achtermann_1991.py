# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Achtermann, Magnus & Bose (1991): refractivity virial coefficients.

- Table II (``tbl_dffe72791eb0`` and its continuation ``tbl_5293ebc0b026``, cells files
  only): per gas and reference, the temperature, A_R (cm^3/mol), B_R (cm^6/mol^2), C_R
  and D_R, with "+-" uncertainties. A_R goes to ``A_R_<reference>_<gas>`` (rows are
  the temperatures), B_R to ``B_R_<reference>_<T>_<gas>`` (one file per temperature,
  as the dataset carries the temperature). C_R and D_R go to ``C_R_...`` and
  ``D_R_...`` the same way (the quantity's order 3 and 4; OEP 4). The ethylene and
  ethane rows are skipped (not species). The continuation table's first row (B_R =
  0.036, "Theory (Ref. 36)", no gas or temperature) belongs to H2 (Orcutt and Cole's
  calculated B_eps, as the caption says), checked on the images. Kirouac and Bose's
  helium B_R carries the footnote mark "a" (their misprinted uncertainty, corrected
  in print here), stripped. Temperatures printed "..." (not stated) become the row id
  "unknown".
- Table III (``tbl_b69cbca038de``, cells file only): the static dielectric virial
  coefficients B_eps of other authors (with their reference numbers) and the authors'
  extrapolations (mark "a"), one file per gas, reference and temperature. The B_R
  column (this work, repeated, and estimates) and the derived B_OR, B_FIR, B_IR
  columns are not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Achtermann_1991_RefractivityVirialCoefficients"
TABLE_II = ("tbl_dffe72791eb0", "tbl_5293ebc0b026")
TABLE_III = "tbl_b69cbca038de"
GASES = {
    "CH_4": "CH4",
    "CO_2": "CO2",
    "SF_6": "SF6",
    "H_2": "H2",
    "N_2": "N2",
    "He": "He",
    "Ar": "Ar",
    "C_2H_4": None,
    "C_2H_6": None,
}
# Table II reference label -> dataset tag.
REFERENCES = {
    "This work": "this_work",
    "Achtermann et al. (Ref. 25)": "achtermann_1986",
    "Montixi and Coulon (Ref. 26)": "montixi_coulon_1983",
    "Olson (Ref. 9)": "olson_1975",
    "Buckingham and Graham (Ref. 12)": "buckingham_graham_1974",
    "St-Arnaud and Bose (Ref. 13)": "st_arnaud_bose_1976",
    "St-Arnaud and Bose (Ref. 14)": "st_arnaud_bose_1979",
    "Burns et al. (Ref. 6)": "burns_1986",
    "Jaeschke (Ref. 28)": "jaeschke_1986",
    "Diller (Ref. 8)": "diller_1968",
    "Montixi et al. (Ref. 7)": "montixi_1983",
    "Orcutt and Cole (Ref. 36)": "orcutt_1967",
    "Kirouac and Bose (Ref. 55)": "kirouac_bose_1976",
    "Coulon et al. (Ref. 29)": "coulon_1981",
    "Theory (Ref. 6)": "theory_burns_1986",
    "Theory (Ref. 7)": "theory_montixi_1983",
    "Theory (Ref. 12)": "theory_buckingham_graham_1974",
    "Theory (Ref. 13)": "theory_st_arnaud_bose_1976",
    "Theory (Ref. 14)": "theory_st_arnaud_bose_1979",
    "Theory (Ref. 15)": "theory_bose_1988",
    "Theory (Ref. 36)": "theory_orcutt_1967",
    "Theory (Ref. 59)": "theory_dacre_1982",
}
# The caption's exceptions to the vacuum wavelength 632.99 nm.
WAVELENGTHS = {
    "olson_1975": "546.2",
    "diller_1968": "546.2",
    "orcutt_1967": "inf",
    "theory_orcutt_1967": "inf",
    "kirouac_bose_1976": "inf",
    "theory_dacre_1982": "inf",
}
TABLE_III_REFERENCES = {
    "34": "bose_cole_1970",
    "45": "bose_1972",
    "46": "bose_cole_1971",
    "48": "huot_bose",
    "49": "hosticka_bose_1974",
}
_REF = re.compile(r"^(.*?)\s*\(Ref\. (\d+)\)$")


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``6.08 +- 0.10`` (marks dropped)."""
    value, _, sigma = cell.partition("\N{PLUS-MINUS SIGN}")
    return number(value.strip().removesuffix("a")), sigma.strip().removesuffix("a")


def _table_ii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas: str | None = None
    for table_id in TABLE_II:
        for cells in json_table(DIRECTORY, table_id)[1:]:
            if cells[0]:
                gas = GASES[cells[0]]
            if cells[6] == "Theory (Ref. 36)" and not cells[1]:
                gas = "H2"  # the continuation row, see the module docstring
            if gas is None:
                continue
            tag = REFERENCES[cells[6]]
            temperature = cells[1] if cells[1] not in {"...", ""} else "unknown"
            wavelength = WAVELENGTHS.get(tag, "632.99")
            printed = " | ".join(cells)
            if cells[2]:
                value, sigma = split(cells[2])
                files.setdefault(
                    f"{gas}/achtermann_1991__A_R_{tag}_{gas}.csv", []
                ).append(Point(temperature, wavelength, value, sigma, printed))
            for column, coefficient in ((3, "B_R"), (4, "C_R"), (5, "D_R")):
                if cells[column]:
                    value, sigma = split(cells[column])
                    name = f"{coefficient}_{tag}_{temperature}_{gas}"
                    path = f"{gas}/achtermann_1991__{name}.csv"
                    assert path not in files, path
                    files[path] = [
                        Point(temperature, wavelength, value, sigma, printed)
                    ]
    return files


def _table_iii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    gas: str | None = None
    for cells in json_table(DIRECTORY, TABLE_III)[1:]:
        if cells[0]:
            gas = GASES[cells[0]]
        if gas is None or not cells[2]:
            continue
        match = _REF.match(cells[2])
        if match:
            text, tag = match.group(1), TABLE_III_REFERENCES[match.group(2)]
        else:
            text, tag = cells[2], "extrapolated"
        value, sigma = split(text)
        path = f"{gas}/achtermann_1991__B_eps_{tag}_{cells[1]}_{gas}.csv"
        assert path not in files, path
        files[path] = [Point(cells[1], "inf", value, sigma, " | ".join(cells))]
    return files


def extract() -> dict[str, list[Point]]:
    """Return A_R, B_R, C_R and D_R of Table II and B_eps of Table III."""
    return _table_ii() | _table_iii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

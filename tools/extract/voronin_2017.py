# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Check Voronin & Zheltikov (2017): the recorded GSE terms against Table 1.

Table 1 (``tbl_0899b125f94c``, JSON only; checked against the image) lists the 15
bands of the generalized Sellmeier equation: A_1r, A_2r, lambda_1r, lambda_2r and the
density N_r (10^15 cm^-3, printed once per species in a merged cell). The library's
cells lose the superscripts of ``4.051 x 10^-6`` (the text reads ``4.051∙10-6``);
the script reads the exponent back and writes the number in e-notation with the
digits unchanged. The record keeps the CO2 rows (r = 1-4) and the N2, O2 and Ar
electronic rows (r = 12-14) as formulas without a form; water is not a curated
species. This script writes no data file: ``extract()`` asserts that every recorded
coefficient equals the table's cell, so the library test re-checks the formulas.
"""

import re
import tomllib

from _json_table import json_table
from _library import DATA, Point

DIRECTORY = "Voronin_2017_GeneralizedSellmeierEquation"
TABLE_1 = "tbl_0899b125f94c"
_POWER = re.compile(r"^(?P<mantissa>\d+\.\d+)∙\s*10(?P<exponent>-\d+)$")
# The rows each formula takes, and the suffix of its coefficient names.
FORMULAS = {
    "gse_CO2_infrared": {"1": "_r1", "2": "_r2", "3": "_r3", "4": "_r4"},
    "gse_N2": {"12": ""},
    "gse_O2": {"13": ""},
    "gse_Ar": {"14": ""},
}


def _number(cell: str) -> str:
    """Return a cell's number in e-notation, the printed digits unchanged."""
    text = cell.strip()
    match = _POWER.match(text)
    if match:
        return f"{match['mantissa']}e{match['exponent']}"
    float(text)  # a plain number, or an error
    return text


def rows() -> dict[str, list[str]]:
    """Return Table 1's rows by band number r."""
    table = json_table(DIRECTORY, TABLE_1)
    header = table[0]
    assert header[0] == "r" and header[3].startswith("A1r"), header
    assert header[5].startswith("λ1r") and header[7].startswith("Nr"), header
    return {cells[0]: cells for cells in table[1:]}


def extract() -> dict[str, list[Point]]:
    """Assert that the record's formula coefficients equal Table 1; return no files."""
    record = tomllib.loads(
        (DATA / "sources" / "voronin_2017.toml").read_text(encoding="utf-8")
    )
    formulas = {f["id"]: f for f in record["formulas"]}
    table = rows()
    assert sorted(table, key=int) == [str(r) for r in range(1, 16)]
    assert set(formulas) == set(FORMULAS)
    for formula_id, bands in FORMULAS.items():
        coefficients = formulas[formula_id]["coefficients"]
        expected = {}
        for band, suffix in bands.items():
            cells = table[band]
            names = ("A1", "A2", "lambda1", "lambda2")
            for name, cell in zip(names, cells[3:7], strict=True):
                expected[f"{name}{suffix}"] = _number(cell)
        # N_r is printed once per species, in the species' first row (10^15 cm^-3).
        first = table[next(iter(bands))]
        expected["N"] = f"{_number(first[7])}e15"
        assert coefficients == expected, (formula_id, coefficients, expected)
    return {}


def main() -> None:
    """Run the check."""
    extract()
    print("voronin_2017: the recorded GSE terms equal Table 1")


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Check Margoliash & Meath (1978): the recorded pseudospectra against Table I.

Table I (``tbl_0daae09a7a5d``, JSON only, checked against the image) lists ten
pseudostates (epsilon_i, f_i) per species in blocks of four species, each number
printed as mantissa(power of ten). The source record keeps the pseudospectra of H2,
N2, O2, N2O and CH4 as formulas (``oscillator`` form, terms in e-notation with the
printed digits). This script writes no data file: ``extract()`` asserts that every
recorded term equals the library's cell, so the library test re-checks the formulas.
"""

import re
import tomllib

from _json_table import json_table
from _library import DATA, Point

DIRECTORY = "Margoliash_1978_PseudospectralDipoleOscillator"
TABLE_I = "tbl_0daae09a7a5d"
SPECIES = {"H_2": "H2", "N_2": "N2", "O_2": "O2", "N_2O": "N2O", "CH_4": "CH4"}
_CELL = re.compile(r"^(?P<mantissa>\d+\.\d+)\((?P<exponent>-?\d+)\)$")


def _e(cell: str) -> str:
    """Return a printed ``mantissa(exponent)`` cell in e-notation."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    return f"{match['mantissa']}e{match['exponent']}"


def pseudospectra() -> dict[str, tuple[list[str], list[str]]]:
    """Return Table I's (energies, strengths) per curated species, in e-notation."""
    blocks: dict[str, tuple[list[str], list[str]]] = {}
    names: list[str] = []
    for cells in json_table(DIRECTORY, TABLE_I):
        if any(_CELL.match(c.strip()) for c in cells if c):
            for i, name in enumerate(names):
                if name in SPECIES and cells[2 * i]:
                    energies, strengths = blocks.setdefault(SPECIES[name], ([], []))
                    energies.append(_e(cells[2 * i]))
                    strengths.append(_e(cells[2 * i + 1]))
        elif cells[0] and "epsilon" not in cells[0]:
            names = cells[0::2]  # a row of species names heads each block
    return blocks


def extract() -> dict[str, list[Point]]:
    """Assert that the record's formula terms equal Table I; return no files."""
    record = tomllib.loads(
        (DATA / "sources" / "margoliash_1978.toml").read_text(encoding="utf-8")
    )
    formulas = {f["material"]: f for f in record["formulas"]}
    table = pseudospectra()
    assert set(table) == set(SPECIES.values()) == set(formulas)
    for gas, (energies, strengths) in table.items():
        assert len(energies) == len(strengths) == 10, gas
        # Every formula's printed field lists the cells verbatim; the terms, where a
        # form is given, are the same numbers in e-notation.
        printed = formulas[gas]["printed"]
        for energy, strength in zip(energies, strengths, strict=True):
            pair = f"({energy.replace('e', '(')}), {strength.replace('e', '(')})"
            assert pair in printed, (gas, pair)
        if "terms" in formulas[gas]:
            assert formulas[gas]["terms"]["E"] == energies, gas
            assert formulas[gas]["terms"]["f"] == strengths, gas
    return {}


def main() -> None:
    """Run the check."""
    extract()
    print("margoliash_1978: the recorded pseudospectra equal Table I")


if __name__ == "__main__":
    main()

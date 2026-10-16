# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Buckingham & Dunmur (1968): dc Kerr second hyperpolarizabilities, Table 1.

Table 1 (``tbl_23a492df71b3``, CSV; image checked) gives, per gas at 23 C, the second
and third virial coefficients used for the densities, 10^36 gamma (esu) and 10^14 A_K
(esu, the molar Kerr constant at low density). The paper prints decimal points as raised
dots ("0·027"), kept in the ``printed`` column and written as "." in the numbers. One
file per gas for gamma (``gamma_<gas>``, row "632.8": the He-Ne line, the only
wavelength); A_K, a restatement of gamma by Eq. (6), and the virial coefficients stay in
the record's notes. Table 2 (second Kerr virial coefficients, polarizabilities) is also
recorded in the notes only.
"""

from _library import Point, number, table, write

DIRECTORY = "Buckingham_1968_KerrEffectInert"
TABLE = "tbl_23a492df71b3"
GASES = {
    "helium": "He",
    "neon": "Ne",
    "argon": "Ar",
    "krypton": "Kr",
    "xenon": "Xe",
    "sulphur hexafluoride": "SF6",
}
DOT = "\N{MIDDLE DOT}"


def extract() -> dict[str, list[Point]]:
    """Return one gamma file per gas."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = GASES[cells[0]]
        value, error = cells[3].replace(DOT, ".").split("\N{PLUS-MINUS SIGN}")
        files[f"{gas}/buckingham_1968__gamma_{gas}.csv"] = [
            Point("632.8", "632.8", number(value), number(error), " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

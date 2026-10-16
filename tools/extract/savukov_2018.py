# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Savukov (2018): CI+MBPT Ar I oscillator strengths from the ground state.

J. Phys. B 51, 065006 (2018), 7 pp. = PDF pages 2-8 (page 1 is the publisher's cover);
a vector PDF with a text layer.

- Table 5 (p. 4, ``tbl_09c26df39d6e``): the CI-MBPT absorption oscillator strengths
  ``f_CI-MBPT`` of the 12 transitions from the ground state 3p6 1S0 to the J = 1 odd
  levels 93750-121161 cm-1. Each row is the printed "NIST energy" (cm-1, unique); ``x``
  is the upper level's Racah label nl[K]J (a prime for the 2P1/2 core), that of the
  NIST level at the printed energy.

**Checks** (:func:`checks`; the script fails if one is missed): every printed NIST
energy is NIST's level energy (``external/nist_asd/Ar_I_levels.tsv``) of a J = 1 odd
level to the printed digits; the printed configuration (core and nl) is that level's
but for rows 9 and 10 (119847.77 and 119760.1725 cm-1), whose configurations are
swapped against NIST and Table 1; the column ``f_expt [30]`` is chan_1992's "present"
Ar values (TC-1012's neighbour record, ``Ar/chan_1992__present_Ar.csv``), row by row in
order, value and uncertainty; Tables 2-4's experimental energies of the 4p and 5p
levels are NIST's to the printed digits, and Tables 6-8's energies each identify one
NIST level.
"""

import csv
import re

from _library import DATA, Point, number, table, write

DIRECTORY = "Savukov_2018_CIMBPTCalculations"
TABLE_5 = "tbl_09c26df39d6e"
TABLES_2_4 = ("tbl_d4192b4591f8", "tbl_dde997c106e1", "tbl_5c91053e779e")
TABLES_6_8 = ("tbl_c887cd18ba52", "tbl_e482284f1b9b", "tbl_dd0aee94106d")
FILE_F = "Ar/savukov_2018__f_ground.csv"
SWAPPED = {"119847.77", "119760.1725"}  # Table 5 rows 9 and 10
_CONFIGURATION = re.compile(r"\((\d)/2\)\s*(\d+[spdf])$")


def _nist_levels() -> list[tuple[str, str, str, str]]:
    """NIST's Ar I levels as (core numerator, nl, Racah label, printed energy)."""
    path = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
    levels = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in list(csv.reader(handle, delimiter="\t"))[1:]:
            conf, term, j, _, energy = (cell.strip('"') for cell in row[:5])
            core = re.search(r"\(2P\*<(\d)/2>\)\.(\d+[spdf])$", conf)
            k = re.match(r"2\[(\d/2)\]", term)
            if core and k and energy:
                prime = "'" if core.group(1) == "1" else ""
                label = f"{core.group(2)}{prime}[{k.group(1)}]{j}"
                levels.append((core.group(1), core.group(2), label, energy))
    return levels


def _level(energy: str, odd: bool | None = None) -> tuple[str, str, str, str]:
    """Return the one NIST level at a printed energy (to its printed digits)."""
    printed = float(energy)
    digits = len([*energy.split("."), ""][1])
    tolerance = 0.5 * 10.0**-digits + 1e-9
    found = [
        level
        for level in _nist_levels()
        if abs(float(level[3]) - printed) <= tolerance
        and (odd is None or (level[1][-1] in "sd") == odd)
    ]
    if len(found) != 1:
        msg = f"{energy} cm-1: {found}"
        raise ValueError(msg)
    return found[0]


def _configuration(cell: str) -> tuple[str, str]:
    """Return a printed configuration's core numerator and nl (``3/2``, ``4s``)."""
    text = cell.replace("◦", "").replace(" ", "")
    match = re.search(r"P(\d)2\)(\d+[spdf])$", text) or _CONFIGURATION.search(text)
    if not match:
        msg = f"not a configuration: {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def extract() -> dict[str, list[Point]]:
    """Return Table 5's CI-MBPT oscillator strengths."""
    points = [
        Point(
            cells[1],
            _level(cells[1], odd=True)[2],
            number(cells[2]),
            "",
            " | ".join(cells),
        )
        for cells in table(DIRECTORY, TABLE_5)[1:]
    ]
    return {FILE_F: points}


def checks() -> None:
    """Check Table 5's energies, configurations and f_expt, and Tables 2-4, 6-8."""
    rows = table(DIRECTORY, TABLE_5)[1:]
    swapped = set()
    for cells in rows:
        core, shell, _, _ = _level(cells[1], odd=True)
        if _configuration(cells[0]) != (core, shell):
            swapped.add(cells[1])
    if swapped != SWAPPED:
        msg = f"Table 5 configurations off NIST: {sorted(swapped)}"
        raise ValueError(msg)
    chan = DATA / "values" / "absorption" / "Ar" / "chan_1992__present_Ar.csv"
    with chan.open(newline="", encoding="utf-8") as handle:
        present = list(csv.DictReader(handle))
    for cells, row in zip(rows, present, strict=False):
        expected = f"{row['value']}({row['sigma']})"
        if cells[3].replace(" ", "") != expected:
            msg = (
                f"Table 5 {cells[1]}: f_expt {cells[3]!r} against chan_1992 {expected}"
            )
            raise ValueError(msg)
    for table_id in TABLES_2_4:
        for cells in table(DIRECTORY, table_id)[1:]:
            if cells[1] != "0":
                core, shell, _, _ = _level(cells[1], odd=False)
                if _configuration(cells[0]) != (core, shell):
                    msg = f"{table_id} {cells[1]}: {cells[0]!r} is not NIST's level"
                    raise ValueError(msg)
    for table_id in TABLES_6_8:
        for cells in table(DIRECTORY, table_id)[1:]:
            _level(cells[0], odd=False)
            _level(cells[1], odd=True)


def main() -> None:
    """Check the tables and write the data file."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

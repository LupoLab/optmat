# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Werbowy (2026): static polarisabilities of the Ar 4s and 4p levels.

J. Quant. Spectrosc. Radiat. Transfer 362, 110041 (2026), 6 pages. The dc Stark shifts
of 21 Ar I 4s-4p lines in fields up to 500 kV/cm give the differences Delta alpha_0 of
the scalar polarisabilities of the two levels (Table 1, ``tbl_38bbff7ecd10``, p. 4);
the scalar polarisability alpha_0 of each of the 14 levels is then fitted so that its
differences reproduce the measured ones (Table 3, ``tbl_da793b53e5e2``, p. 5, column
'Fitted to the experiment', with its uncertainty in parentheses). Units as printed,
kHz/(kV/cm)^2, alpha/h with the level shift -(1/2) alpha F^2 (Eq. 1).

One ``level_polarizability`` file per level, static (``x`` = inf), the row the level's
energy as printed (digit-group spaces removed):

* ``alpha0_<level>``: the fitted scalar alpha_0 and its uncertainty, 14 levels;
* ``alpha2_<level>``: the tensor bound, 'The experimental values alpha_2 were found to
  be less than 10 kHz/(kV/cm)2' (Table 3 note a), for the 11 levels with J >= 1 (alpha_2
  vanishes for J = 0): the value 0 with 10 as a bound (no line splitting was seen).

The calculated alpha_0 and alpha_2 columns of Table 3 (NIST, CI-MBPT and measured line
strengths) are not extracted. Both tables have merged header cells and no CSV; they
are read from the JSON cells (``_json_table``).

:func:`checks` tests the paper against itself and NIST: every 'Fitted to the
experiment' Delta alpha_0 of Table 1 (30 lines) is the difference of the fitted alpha_0
of its upper and lower level, exactly; the abstract's ranges of alpha_0 (56-96 for the
4s levels, 95-327 for the 4p levels) are the table's; every level energy, which the
paper takes from NIST (Table 1 note a), is NIST's (``data/external/nist_asd/
Ar_I_levels.tsv``) to the printed three decimals.

Every cell taken was checked against the table images (2026-10-09).
"""

import csv
import re
from typing import Final

from _json_table import json_table
from _library import DATA, Point, number, write

DIRECTORY: Final = "Werbowy_2026_StaticStarkEffect"
LEVELS_TABLE: Final = "tbl_da793b53e5e2"
"""Table 3: alpha_0 and alpha_2 of the 4s and 4p levels."""
LINES_TABLE: Final = "tbl_38bbff7ecd10"
"""Table 1: Delta alpha_0 of the lines, calculated, fitted and measured."""
TENSOR_BOUND: Final = "10"
"""Table 3 note a: the experimental alpha_2 'less than 10 kHz/(kV/cm)2'."""
TENSOR_NOTE: Final = (
    "Table 3 note a: The experimental values alpha_2 were found to be less than 10"
    " kHz/(kV/cm)2"
)
ABSTRACT_RANGES: Final = {"4s": (56, 96), "4p": (95, 327)}
"""'in a range of 56 up to 96 ... and 95 up to 327 kHz/(kV/cm)2' (abstract)."""
_NUMBER: Final = re.compile(r"^\s*(-?\d+)\s*\((\d+)\)\s*$")


def _energy(cell: str) -> str:
    """Return a printed level energy without its digit-group spaces."""
    return number(re.sub(r"\s", "", cell))


def _levels() -> list[tuple[str, str, list[str]]]:
    """Return Table 3's level rows: the energy, the configuration (4s or 4p), cells."""
    rows = []
    configuration = ""
    for cells in json_table(DIRECTORY, LEVELS_TABLE)[3:]:
        if cells[0].startswith("3s^23p^5"):
            configuration = cells[0].removeprefix("3s^23p^5")
            continue
        rows.append((_energy(cells[0]), configuration, cells))
    if len(rows) != 14:
        msg = f"Table 3: {len(rows)} level rows, expected 14"
        raise ValueError(msg)
    return rows


def _fitted(cell: str) -> tuple[str, str]:
    """Return a 'value (uncertainty)' cell's two numbers as printed."""
    match = _NUMBER.match(cell)
    if match is None:
        msg = f"Table 3: cannot read {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def extract() -> dict[str, list[Point]]:
    """Return one file per level for alpha_0, and per J >= 1 level for alpha_2."""
    files = {}
    for energy, _, cells in _levels():
        name = energy.partition(".")[0]
        value, sigma = _fitted(cells[9])
        printed = " | ".join(cells)
        files[f"Ar/werbowy_2026__alpha0_{name}.csv"] = [
            Point(energy, "inf", value, sigma, printed)
        ]
        if cells[2] != "0":
            files[f"Ar/werbowy_2026__alpha2_{name}.csv"] = [
                Point(energy, "inf", "0", TENSOR_BOUND, f"{printed} | {TENSOR_NOTE}")
            ]
    return files


def _nist() -> list[float]:
    """Return NIST's Ar I level energies (cm^-1)."""
    path = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    energies = []
    for row in rows:
        try:
            energies.append(float(row["Level (cm-1)"]))
        except ValueError:
            continue
    return energies


def checks() -> None:
    """Check the fitted differences, the abstract's ranges and NIST; raise if not."""
    levels = _levels()
    fitted = {energy: int(_fitted(cells[9])[0]) for energy, _, cells in levels}
    lines = 0
    for cells in json_table(DIRECTORY, LINES_TABLE)[3:]:
        upper, lower = _energy(cells[1]), _energy(cells[2])
        difference = int(number(cells[7]))
        if difference != fitted[upper] - fitted[lower]:
            msg = f"Table 1 {cells[0]}: fitted {difference}, Table 3 differs"
            raise ValueError(msg)
        lines += 1
    if lines != 30:
        msg = f"Table 1: {lines} lines, expected 30"
        raise ValueError(msg)
    for configuration, (low, high) in ABSTRACT_RANGES.items():
        values = [fitted[e] for e, c, _ in levels if c == configuration]
        if (min(values), max(values)) != (low, high):
            msg = f"{configuration}: alpha_0 {min(values)}-{max(values)}"
            raise ValueError(msg)
    nist = _nist()
    for energy in fitted:
        if min(abs(float(energy) - level) for level in nist) > 0.0005:
            msg = f"level {energy} cm^-1 is not NIST's"
            raise ValueError(msg)


def main() -> None:
    """Check the paper, then write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

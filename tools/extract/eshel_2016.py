# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Eshel, Rice & Perram (2016): Ar s-p line broadening and shift in Ar-He.

Table 1 (p. 47, library table ``tbl_45fc5023b655``), 'Pressure broadening rate, 2 gamma
(MHz/Torr)', and Table 2 (p. 49, ``tbl_6d2606d4248c``), 'Pressure shifting rate, beta
(MHz/Torr)', 'Negative rates are red-shifts', print this work's rates for four Ar-He
mixtures (2.5, 5, 10 and 25 % Ar; 440 +- 20 K; 10-160 Torr total pressure) and the three
lines 801.699 nm (s5 -> p8), 800.836 nm (s4 -> p6) and 795.036 nm (s3 -> p4), as
'value +- 95 % confidence interval' ('-' for the 800.836 nm line at 25 % Ar). One file
per mixture and effect, the ``pressure_line_coefficient``; rows are the lines. The rows
quoted from the literature (refs. 21-24, 27) are not taken. The cells were compared with
the table images (2026-10-12): no misread cell.

**Checks** (:func:`checks`): the 10 % Ar rows must be the abstract's and the captions of
Figs. 6 and 7 (p. 40, p. 46): broadening 14.18 +- 0.65, 17.85 +- 0.78, 16.59 +- 1.22 and
shift 1.81 +- 0.30, 0.72 +- 0.32, 2.94 +- 0.48 MHz/Torr; every mixture row must state
440 +- 20 K and 10-160 Torr. Reported: each line's vacuum wavelength against NIST's (the
record notes).
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Eshel_2016_PressureBroadeningShift"
TABLES = {"broadening": "tbl_45fc5023b655", "shift": "tbl_6d2606d4248c"}
HEADINGS = {
    "broadening": r"Pressure broadening rate, 2\gamma (MHz/Torr)",
    "shift": r"Pressure shifting rate, \beta (MHz/Torr)",
}
LINES = {
    "801.699 nm(s_5 \\to p_8)": ("801.699", "801.699 nm (s5 -> p8)"),
    "800.836 nm(s_4 \\to p_6)": ("800.836", "800.836 nm (s4 -> p6)"),
    "795.036 nm(s_3 \\to p_4)": ("795.036", "795.036 nm (s3 -> p4)"),
}
"""Table column headings: the record's row id and label."""
MIXTURES = {"2.5% Ar": "mix2p5", "5% Ar": "mix5", "10% Ar": "mix10", "25% Ar": "mix25"}
"""This work's rows and the dataset id suffixes."""
_CELL = re.compile(r"^(-?\d+\.\d+) ± (\d+\.\d+)$")
TEXT_10PC = {
    "broadening": {
        "801.699": ("14.18", "0.65"),
        "800.836": ("17.85", "0.78"),
        "795.036": ("16.59", "1.22"),
    },
    "shift": {
        "801.699": ("1.81", "0.30"),
        "800.836": ("0.72", "0.32"),
        "795.036": ("2.94", "0.48"),
    },
}
"""The abstract (p. 40) and the captions of Figs. 6 and 7 (p. 46), 10 % Ar."""
NIST_VACUUM_NM = {"801.699": 801.6990, "800.836": 800.8359, "795.036": 795.0362}
"""NIST observed vacuum wavelengths (data/external/nist_asd/Ar_I_lines_nir.tsv)."""


def file(effect: str, mixture: str) -> str:
    """Return the data file of one effect and one mixture."""
    return f"Ar/eshel_2016__{effect}_{MIXTURES[mixture]}.csv"


def _grid(effect: str) -> list[list[str]]:
    """Return a table's rows, the heading, columns and this work's rows checked."""
    grid = json_table(DIRECTORY, TABLES[effect])
    if (
        grid[1][0] != HEADINGS[effect]
        or grid[0][4:7] != list(LINES)
        or [row[0] for row in grid[2:6]] != list(MIXTURES)
        or any(row[7] != "This work" for row in grid[2:6])
    ):
        msg = f"Table {effect} has changed: {grid[:6]}"
        raise ValueError(msg)
    return grid


def extract() -> dict[str, list[Point]]:
    """Return this work's rates, per effect and mixture, rows by line."""
    files: dict[str, list[Point]] = {}
    for effect in TABLES:
        grid = _grid(effect)
        for row in grid[2:6]:
            points = []
            for column, (line, label) in zip(grid[0][4:7], LINES.values(), strict=True):
                cell = row[grid[0].index(column)]
                if cell == "-":
                    continue
                match = _CELL.match(cell)
                if match is None:
                    msg = f"{effect} {row[0]}: unexpected cell {cell!r}"
                    raise ValueError(msg)
                points.append(
                    Point(
                        line,
                        label,
                        match.group(1),
                        match.group(2),
                        " | ".join([row[0], row[1], row[2], column, cell]),
                    )
                )
            files[file(effect, row[0])] = points
    return files


def checks() -> dict[str, str]:
    """Raise unless the tables agree with the abstract and captions; report more."""
    files = extract()
    failures = []
    for effect, expected in TEXT_10PC.items():
        got = {p.row: (p.value, p.sigma) for p in files[file(effect, "10% Ar")]}
        if got != expected:
            failures.append(f"{effect} 10 % Ar: {got} against the text's {expected}")
        for row in _grid(effect)[2:6]:
            if row[1] != "440 ± 20" or row[2] != "10-160":
                failures.append(f"{effect} {row[0]}: conditions {row[1:3]}")
    counts = {name: len(points) for name, points in files.items()}
    if sum(counts.values()) != 22:
        failures.append(f"expected 22 values, got {counts}")
    if failures:
        raise ValueError("; ".join(failures))
    return {
        line: f"printed {line} nm, NIST vacuum {nist:.4f} nm"
        for line, nist in NIST_VACUUM_NM.items()
    }


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

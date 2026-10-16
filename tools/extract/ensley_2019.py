# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ensley & Bambha (2019): n2 in the mid-wave infrared, and literature values.

The three tables are published only as JSON; every cell was checked against the table
images.

- Table 1 (``tbl_0d115aaf67aa``): semiconductors, band gap (eV, "taken from Ref.
  [19]") and n2 (1e-14 cm2/W) at 2 and 3.9 um. ZnSe at 3.9 um is "a previously
  published value from Ref. [24]" (caption) and is split into its own, reproduced file.
- Table 2 (``tbl_67677bf4d75b``, unlabelled in the library): fluoride crystals and
  glasses. STRUCTURAL FIX: the extracted header reads "h_2(10^{-6} cm^{-3}/W)", a
  mis-reading of the printed "n2 (10^-16 cm^2/W)" (image, and paper.md's display line
  above the table); the unit is set in the record, not read from the header. The
  columns are the Eq. (1) (Boling-Glass-Owyoung) prediction, then 532 nm, 800 nm,
  1.06 um and 2 um. A cell holds several values, each marked by a superscript letter
  naming its reference (the caption's list, a to m); an unmarked value is this work
  (800 nm and 2 um). A second mark "+", "#" or "*" means a measurement at 525 nm,
  810 nm or 1.03 um (caption), which then replaces the column's wavelength. The
  1.06-um column is written as 1060 nm.
- Table 3 (``tbl_a359645a6ee3``): chalcogenide glasses, n2 (1e-14 cm2/W) at 2 and
  3.9 um.

The pulse duration of each of the paper's own measurements, from the text (p. 37941),
is the per-row column ``pulse_fs``: 90 fs at 2 um, 240 fs at 3.9 um, and ~100 fs (the
Ti:sapphire amplifier's pulse) at 800 nm.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Ensley_2019_UltrafastNonlinearRefraction"
TABLE_1 = "tbl_0d115aaf67aa"
TABLE_2 = "tbl_67677bf4d75b"
TABLE_3 = "tbl_a359645a6ee3"

PULSE_FS = {"800": "100", "2000": "90", "3900": "240"}

SEMICONDUCTORS = {
    "ZnS": "ZnS",
    "ZnSe": "ZnSe",
    "GaAs": "GaAs",
    "InP": "InP",
    "GaSb": "GaSb",
    "Si": "Si",
    "Ge": "Ge",
}
FLUORIDES_AND_GLASSES = {
    r"\text{BaF}_{2}": "BaF2",
    r"\text{CaF}_{2}": "CaF2",  # codespell:ignore
    r"\text{MgF}_{2}": "MgF2",
    "LiF": "LiF",
    "Fused Silica": "SiO2",
    "BK7": "BK7",
}
CHALCOGENIDES = {
    "AMTIR-1 (Ge_{33}As_{12}Se_{55})": "AMTIR-1",
    "AMTIR-2 (AsSe)": "AMTIR-2",
    "AMTIR-4 (AsSe)": "AMTIR-4",
    "AMTIR-5 (AsSe)": "AMTIR-5",
    "AMTIR-6 (As_{2}S_{3})": "AMTIR-6",
}

# Table 2's columns after the Eq. (1) one, as wavelengths in nm.
TABLE_2_COLUMNS = {2: "532", 3: "800", 4: "1060", 5: "2000"}
# The caption's second marks: measured at 525 nm, 810 nm and 1.03 um.
MARKED_WAVELENGTH = {"+": "525", "#": "810", "*": "1030"}

_ENTRY = re.compile(r"^([-\d.]+)(?:\^\{\\text\{([a-m])(?:,([+#*]))?\}\})?$")


def _entries(cell: str) -> list[tuple[str, str, str]]:
    """Split a Table 2 cell into (value, reference letter or "", second mark or "")."""
    entries = []
    for part in cell.split(", "):
        match = _ENTRY.match(part.strip())
        if match is None:
            msg = f"unexpected Table 2 entry {part!r}"
            raise ValueError(msg)
        value, letter, mark = match.groups()
        entries.append((number(value), letter or "", mark or ""))
    return entries


def _own(wavelength: str, value: str, printed: str) -> Point:
    """Return one of the paper's own measurements, with its pulse duration."""
    extra = {"pulse_fs": PULSE_FS[wavelength]}
    return Point(wavelength, wavelength, number(value), "", printed, extra)


def _table_1() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_1)[2:]:
        material = SEMICONDUCTORS[cells[0]]
        printed = " | ".join(cells)
        # "1.12 (indirect)": the value; the label stays in printed.
        gap = number(cells[1].removesuffix(" (indirect)"))
        files[f"{material}/ensley_2019__eg_{material}.csv"] = [
            Point("1", "band gap", gap, "", printed)
        ]
        own = [_own("2000", cells[2], printed)]
        if material == "ZnSe":
            # Caption: ZnSe at 3.9 um is Ref. [24]'s published value.
            files["ZnSe/ensley_2019__n2_ZnSe_werner_2019.csv"] = [
                Point("3900", "3900", number(cells[3]), "", printed)
            ]
        else:
            own.append(_own("3900", cells[3], printed))
        files[f"{material}/ensley_2019__n2_{material}.csv"] = own
    return files


def _table_2() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_2)[2:]:
        material = FLUORIDES_AND_GLASSES[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/ensley_2019__eq1_{material}.csv"] = [
            Point("1", "not stated", number(cells[1]), "", printed)
        ]
        own = []
        literature: dict[str, list[Point]] = {}
        for column, wavelength in TABLE_2_COLUMNS.items():
            for value, letter, mark in _entries(cells[column]):
                if not letter:
                    own.append(_own(wavelength, value, printed))
                    continue
                x = MARKED_WAVELENGTH.get(mark, wavelength)
                points = literature.setdefault(letter, [])
                row = x if all(p.row != x for p in points) else f"{x}-2"
                points.append(Point(row, x, value, "", printed))
        files[f"{material}/ensley_2019__n2_{material}.csv"] = own
        for letter, points in literature.items():
            files[f"{material}/ensley_2019__lit_{material}_{letter}.csv"] = points
    return files


def _table_3() -> dict[str, list[Point]]:
    files = {}
    for cells in json_table(DIRECTORY, TABLE_3)[2:]:
        material = CHALCOGENIDES[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/ensley_2019__n2_{material}.csv"] = [
            _own("2000", cells[1], printed),
            _own("3900", cells[2], printed),
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the data files of Tables 1-3."""
    return {**_table_1(), **_table_2(), **_table_3()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

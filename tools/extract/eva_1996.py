# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Eva & Mann (1996): calorimetric absorption of UV-window materials (Table 1).

Table 1 is a wide table: material, thickness dz, E_gap, then three blocks, one per
wavelength, each opened by a label column ("351 nm:", "248 nm:", "193 nm:") that is
empty in every data row. The columns are, by index: 0 material, 1 dz [mm], 2 E_gap
[eV], 3 (label 351 nm), 4 alpha_3.5eV [1e-3/cm], 5 (label 248 nm), 6 alpha_5eV
[1e-3/cm], 7 beta_5eV [1e-10 cm/W], 8 beta_Lit [1e-10 cm/W], 9 (label 193 nm),
10 alpha_6.4eV [1e-3/cm], 11 beta_6.4eV [1e-10 cm/W].

Checked against the table image (tbl_48a4477a3382.jpg), row by row: every value sits
in the column the image puts it in, the label columns are empty, "-" marks a value
not measured and an empty cell nothing printed (CaF2 5 mm at 193 nm is the last row,
with nothing at 248 nm). Structural handling:

- The row id of a value is its table row (1-9) and wavelength, "6/248".
- "\approx" cells are approximate values (E_gap of crystalline quartz, alpha of
  Suprasil 2 at 193 nm, beta "\approx 0" of CaF2, LiF and MgF2 at 248 nm): the number
  is kept, the raw cell stays in ``printed``, and the record says what it is.
- The E_gap of crystalline quartz carries its citation, "[20]" (Weinberg et al.);
  the others come from Ref. [2] (Taylor et al.), as the caption says.
- The beta_Lit cell of Suprasil 2, ".45 (2.1 [19])", holds two values: Taylor's .45
  (the caption's default source) and 2.1 from Ref. [19] (Brimacombe et al.). They are
  split into their two reproduced datasets.
- "\leq" cells are upper bounds (Taylor's, for LiF and MgF2); kept and excluded.
"""

import re

from _library import Point, table, write

DIRECTORY = "Eva_1996_CalorimetricMeasurementTwo"
TABLE = "tbl_48a4477a3382"

# (material id, dataset suffix) of each table row, 1-9 in the printed order.
ROWS = {
    1: ("SiO2", "SiO2_suprasil1"),  # "Fused silica S1": Suprasil 1 (text, p. 2)
    2: ("CaF2", "CaF2"),  # 2 mm
    3: ("CaF2", "CaF2"),  # 5 mm
    4: ("LiF", "LiF"),
    5: ("MgF2", "MgF2"),
    6: ("SiO2", "SiO2_suprasil2"),  # "Fused silica S2": Suprasil 2 (text, p. 5)
    7: ("SiO2_quartz", "SiO2_quartz"),  # "Cryst. quartz", z-cut (p. 6)
    8: ("BaF2", "BaF2"),
    9: ("CaF2", "CaF2"),  # 5 mm, measured at 193 nm only
}
LABELS = {
    1: "Fused silica S1",
    2: "CaF_2",
    3: "CaF_2",
    4: "LiF",
    5: "MgF_2",
    6: "Fused silica S2",
    7: "Cryst. quartz",
    8: "BaF_2",
    9: "CaF_2",
}
ALPHA_COLUMNS = ((4, "351"), (6, "248"), (10, "193"))
BETA_COLUMNS = ((7, "248"), (11, "193"))
LITERATURE = 8

_NUMBER = re.compile(r"^(?:\\approx|\\leq)?\s*(?P<number>\d*\.?\d+)(?:\s*\[\d+\])?$")


def number(cell: str) -> str:
    r"""Return the number of a cell, dropping ``\approx``, ``\leq`` and a citation."""
    match = _NUMBER.match(cell.strip())
    if match is None:
        msg = f"not a number cell: {cell!r}"
        raise ValueError(msg)
    return match["number"]


def _given(cell: str) -> bool:
    """Return whether a cell prints a value ("-" is 'not measured', "" nothing)."""
    return cell.strip() not in {"", "-"}


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, the alphas, the betas and the literature betas."""
    rows = table(DIRECTORY, TABLE)
    header, body = rows[0], rows[1:]
    if [header[i] for i in (3, 5, 9)] != ["351 nm:", "248 nm:", "193 nm:"]:
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}

    def add(name: str, point: Point) -> None:
        files.setdefault(name, []).append(point)

    for index, cells in enumerate(body, start=1):
        if cells[0] != LABELS[index] or any(cells[i] for i in (3, 5, 9)):
            msg = f"row {index} is not where the image puts it: {cells!r}"
            raise ValueError(msg)
        material, suffix = ROWS[index]
        printed = " | ".join(cells)
        gap = "band_gap_SiO2" if material == "SiO2" else f"band_gap_{material}"
        add(
            f"{material}/eva_1996__{gap}.csv",
            Point(str(index), "", number(cells[2]), "", printed),
        )
        for column, wavelength in ALPHA_COLUMNS:
            if _given(cells[column]):
                add(
                    f"{material}/eva_1996__alpha_{suffix}.csv",
                    Point(
                        f"{index}/{wavelength}",
                        wavelength,
                        number(cells[column]),
                        "",
                        printed,
                    ),
                )
        for column, wavelength in BETA_COLUMNS:
            if _given(cells[column]):
                add(
                    f"{material}/eva_1996__beta_{suffix}.csv",
                    Point(
                        f"{index}/{wavelength}",
                        wavelength,
                        number(cells[column]),
                        "",
                        printed,
                    ),
                )
        literature = cells[LITERATURE].strip()
        if not literature:
            continue
        # ".45 (2.1 [19])": Taylor's value (the caption's default, Ref. [2]) and, in
        # parentheses, Brimacombe et al.'s (Ref. [19]); see the module docstring.
        split = re.fullmatch(
            r"(?P<taylor>\.\d+) \((?P<other>[\d.]+) \[19\]\)", literature
        )
        taylor = split["taylor"] if split else number(literature)
        add(
            f"{material}/eva_1996__taylor_248_{material}.csv",
            Point(f"{index}/248", "248", taylor, "", printed),
        )
        if split:
            add(
                f"{material}/eva_1996__brimacombe_248_{material}.csv",
                Point(f"{index}/248", "248", split["other"], "", printed),
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

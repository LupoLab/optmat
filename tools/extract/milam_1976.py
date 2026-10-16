# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Milam & Weber (1976): n2 and gamma of seven optical materials (Table I).

Table I (``tbl_9e8c390f5b4c``, p. 2499; JSON only: group headings span the row) gives,
per material, the type, n_D, the Abbe number, n2 in 10^-13 esu and gamma in 10^-16
cm^2/W, each with its uncertainty ("1.24 ± 0.12"), measured by time-resolved
interferometry with linearly polarized 125-ps, 1064-nm pulses (the caption). Group
headings ("Optical glasses", "Laser glass", ...) are rows of their own and are skipped.
KH2PO4 occupies two rows (merged material and type cells): n_D and the Abbe number of
the e ray on the first, of the o ray on the second, and one n2 and gamma for the
crystal on the first. The o-ray index is printed "1.509(0)" (the letter o), which is
kept in ``printed`` only.

The table image (``tables/tbl_9e8c390f5b4c.jpg``) matches every extracted cell.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Milam_1976_MeasurementNonlinearRefractive"
TABLE = "tbl_9e8c390f5b4c"
PUMP_NM = "1064"

# Table I "Type" cell -> registered material id.
MATERIALS = {
    "BK 7": "BK7",
    "FC 5": "FK-5",
    "4000": "SiO2",
    "FK 51": "FK-51",
    "ED-2": "ED-2",
    "FR-5": "FR-5",
    "Tetragonal": "KDP",
}
GROUPS = (
    "Optical glasses",
    "Laser glass",
    "Faraday rotator glass",
    "Nonlinear crystal",
)

_CELL = re.compile(r"^([\d.]+)\s*±\s*([\d.]+)$")


def _value(cell: str) -> tuple[str, str]:
    """Return the value and uncertainty of a cell ``v ± u``, as printed."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), number(match.group(2))


def extract() -> dict[str, list[Point]]:
    """Return gamma (cm^2/W) and n2 (esu) of each material."""
    rows = json_table(DIRECTORY, TABLE)
    header = rows[0]
    if header[4] != "n_2(10-13 esu)" or header[5] != "\\gamma(10-16 cm2/W)":
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    seen = []
    for index, cells in enumerate(rows[1:], start=1):
        if cells[0] in GROUPS or not cells[4]:
            continue
        material = MATERIALS[cells[1]]
        seen.append(material)
        printed = " | ".join(cells)
        if material == "KDP":
            # The merged KDP cell: the o-ray n_D and Abbe number are on the next row.
            printed += " || " + " | ".join(rows[index + 1])
        n2, n2_sigma = _value(cells[4])
        gamma, gamma_sigma = _value(cells[5])
        files[f"{material}/milam_1976__gamma_{material}.csv"] = [
            Point("1", PUMP_NM, gamma, gamma_sigma, printed)
        ]
        files[f"{material}/milam_1976__n2_{material}.csv"] = [
            Point("1", PUMP_NM, n2, n2_sigma, printed)
        ]
    if seen != list(MATERIALS.values()):
        msg = f"unexpected materials {seen!r}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

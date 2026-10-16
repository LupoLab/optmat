# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Adair, Chase & Payne (1987): n2 of 27 glasses by three-wave mixing.

Table 2 (``tbl_cd009634988d``, p. 879) gives, per glass, the n2 measured here by nearly
degenerate three-wave mixing at 1.064/1.071 um ("n2 meas. a", footnote a "This work"),
the literature value where there is one ("n2 meas.", footnote c: Ref. 13, Weber, Milam &
Smith 1978; footnote d: Ref. 7, Moran, She & Carman 1975) and the value calculated with
the Boling-Glass-Owyoung formula ("n2 calc. b", footnote b: Ref. 31, Stokowski, Saroyan
& Weber 1981), all in 10^-13 esu (the caption). The "Type" column names the glass
family on one row of each group only (the image shows the family label centred on its
group); it is not a value and is not used. A literature cell carries its footnote
letter as ``^c`` or ``^d``, which the script strips and uses to name the original.

The table image (``tables/tbl_cd009634988d.jpg``) matches every extracted cell.
Table 2 prints "E-133" where Table 1 (compositions) prints "K-133" for the same place
in the fluorophosphate group; both are the registered E-133 (the record says why).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Adair_1987_NonlinearRefractiveIndex"
TABLE = "tbl_cd009634988d"
PUMP_UM = "1.064"

# Table 2 sample name -> registered material id.
MATERIALS = {
    "E-115": "E-115",
    "E-131": "E-131",
    "E-132": "E-132",
    "E-133": "E-133",
    "K-1172": "K-1172",
    "LHG-5": "LHG-5",
    "LHG-6": "LHG-6",
    "Q-88": "Q-88",
    "P-108": "P-108",
    "5037": "Hoya-5037",
    "5038": "Hoya-5038",
    "FR-5": "FR-5",
    "ED-2": "ED-2",
    "Fused silica": "SiO2",
    "C2828": "C-2828",
    "ED-8": "ED-8",
    "S-7606": "S-7606",
    "LG-650": "LG-650",
    "ED-3": "ED-3",
    "Q-246": "Q-246",
    "BK-7": "BK7",
    "L-109": "L-109",
    "L-65": "L-65",
    "9028": "fluorozirconate-9028",
    "K-1261": "K-1261",
    "3151": "Kigre-3151",
    "SF-6": "SF-6",
}

# Footnote letter of a literature cell -> the dataset id prefix.
LITERATURE = {"c": "weber_1978", "d": "moran_1975"}

_LITERATURE_CELL = re.compile(r"^([\d.]+)\^([cd])$")


def extract() -> dict[str, list[Point]]:
    """Return the measured, literature and calculated n2 of each glass."""
    rows = table(DIRECTORY, TABLE)
    header, body = rows[0], rows[1:]
    if header[2] != "n_{2}^{meas.a}" or header[4] != "n_{2}^{calc.b}":
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    if [cells[0] for cells in body] != list(MATERIALS):
        msg = "unexpected rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        material = MATERIALS[cells[0]]
        printed = " | ".join(cells)
        files[f"{material}/adair_1987__n2_{material}.csv"] = [
            Point("1", PUMP_UM, number(cells[2]), "", printed)
        ]
        files[f"{material}/adair_1987__calc_{material}.csv"] = [
            Point("1", "not stated", number(cells[4]), "", printed)
        ]
        if cells[3]:
            match = _LITERATURE_CELL.match(cells[3])
            if match is None:
                msg = f"unexpected literature cell {cells[3]!r}"
                raise ValueError(msg)
            name = LITERATURE[match.group(2)]
            files[f"{material}/adair_1987__{name}_{material}.csv"] = [
                Point("1", PUMP_UM, number(match.group(1)), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

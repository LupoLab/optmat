# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kabacinski et al. (2019): n2 at 1030 nm and the literature values (Table 1).

Table 1 gives, per material (and polarization), this work's n2 at 1030 nm and, in one
cell, one or more literature values, each "value [± sigma] at <wavelength> nm [ref]".
The table image (``tables/tbl_46ff24325236.jpg``) prints each literature value on its
own line; the extracted cell joins them with spaces, so the script splits the cell with
a pattern and checks that nothing is left over. "7.39, 8.00 at 1064 nm [36]" is two
values at one wavelength. The two values "in Yb:YVO4" (Ref. [37]) are for a doped
crystal and are not extracted (source notes).

The unit is in the table footnote, garbled in the extraction ("10^-16 cm2 . W", with a
Unicode minus); the image reads "All values are in units of 10^-16 cm^2/W".
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Kabacinski_2019_NonlinearRefractiveIndex"
TABLE = "tbl_46ff24325236"

# The first cell as printed -> (material id, file-name suffix).
ROWS = {
    "Fused silica": ("SiO2", "SiO2"),
    "CaF2": ("CaF2", "CaF2"),
    "YAG": ("YAG", "YAG"),
    "BiBO E||X": ("BiBO", "BiBO_x"),
    "Calcite no": ("CaCO3", "CaCO3_o"),
    "Calcite ne": ("CaCO3", "CaCO3_e"),
    "YVO4 no": ("YVO4", "YVO4_o"),
    "YVO4 ne": ("YVO4", "YVO4_e"),
}

# Reference number -> dataset-id prefix: [17] Adair et al. (1989), [21] Flom et al.
# (2015), [33] Milam (1998), [34] Milam, Weber & Glass (1978 as printed), [35] Bliss,
# Speck & Simmons (1974), [36] Miller et al. (2008), [37] Selivanov et al. (2006).
REFERENCES = {
    "17": "adair_1989",
    "21": "flom_2015",
    "33": "milam_1998",
    "34": "milam_1978",
    "35": "bliss_1974",
    "36": "miller_2008",
}
DOPED = "37"

_ITEM = re.compile(
    r"(?P<values>[\d.]+(?:,\s*[\d.]+)*)(?:\s*±\s*(?P<sigma>[\d.]+))?"
    r"\s+at\s+(?P<nm>\d+)\s*nm\s*(?:in Yb:YVO4\s*)?\[(?P<ref>\d+)\]\s*"
)

# This work's wavelength (the caption: "Measured in This Work at 1030 nm").
THIS_WORK_NM = "1030"


def extract() -> dict[str, list[Point]]:
    """Return this work's value and each literature reference's values, per row."""
    files: dict[str, list[Point]] = {}
    for index, cells in enumerate(table(DIRECTORY, TABLE)[1:], start=1):
        material, suffix = ROWS[cells[0].strip()]
        printed = " | ".join(cells)
        value, sigma = cells[1].split("±")
        files[f"{material}/kabacinski_2019__this_work_{suffix}.csv"] = [
            Point(str(index), THIS_WORK_NM, number(value), number(sigma), printed)
        ]
        literature = cells[2].strip()
        items = list(_ITEM.finditer(literature))
        if "".join(item.group(0) for item in items).strip() != literature:
            msg = f"row {index}: unparsed literature cell {literature!r}"
            raise ValueError(msg)
        for item in items:
            if item.group("ref") == DOPED:
                continue
            name = f"{REFERENCES[item.group('ref')]}_{suffix}"
            points = files.setdefault(f"{material}/kabacinski_2019__{name}.csv", [])
            for k, text in enumerate(item.group("values").split(","), start=1):
                row = f"{index}" if "," not in item.group("values") else f"{index}.{k}"
                sigma = number(item.group("sigma")) if item.group("sigma") else ""
                points.append(
                    Point(row, item.group("nm"), number(text), sigma, printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

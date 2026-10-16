# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Zhang, Lu & Wang (2008): N2, O2, Ar, CO2 n - 1 at 800, 632.99 nm (Table 2).

Table 2 (``tbl_8d9995b38d24``, merged gas cells, so read from the JSON cells) gives per
gas the paper's value at 800 nm with its standard error and its extrapolation to 632.99
nm (no uncertainty, footnote b), then other authors' values marked ``[Ref.]``, all at
20 degrees C and 101 325 Pa. Air is skipped (not a species of the package). Each
reference becomes one dataset per gas, with rows ``800`` and ``632.99``.
"""

import json
import re

from _library import LIBRARY, Point, number, write

DIRECTORY = "Zhang_2008_PrecisionRefractiveIndex"
TABLE_2 = "tbl_8d9995b38d24"
GASES = {"CO2": "CO2", "N2": "N2", "O2": "O2", "Ar": "Ar"}
REFERENCES = {
    "25": "old_1971",
    "27": "bideau_mehu_1973",
    "14": "birch_1991",
    "28": "simmons_1978",
    "19": "peck_khanna_1966",
    "21": "griesmann_1999",
    "24": "peck_fisher_1964",
}
WAVELENGTHS = ("800", "632.99")
_CELL = re.compile(r"^([\d.]+)(?:\s*±\s*([\d.]+))?(?:\s*\[(\d+)\])?$")


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the paper's values and each reproduced reference, per gas."""
    files: dict[str, list[Point]] = {}
    gas = ""
    for cells in grid(TABLE_2)[1:]:
        gas = cells[0] or gas
        if gas not in GASES:
            continue  # air
        printed = " | ".join(cells)
        for wavelength, cell in zip(WAVELENGTHS, cells[1:], strict=True):
            if not cell:
                continue
            match = _CELL.match(cell)
            if match is None:
                msg = f"cannot read cell {cell!r}"
                raise ValueError(msg)
            value, sigma, ref = match.groups()
            if ref is None:
                name = "measured" if wavelength == "800" else "extrapolated"
            else:
                name = REFERENCES[ref]
            files.setdefault(f"{gas}/zhang_2008__{name}_{gas}.csv", []).append(
                Point(wavelength, wavelength, number(value), sigma or "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

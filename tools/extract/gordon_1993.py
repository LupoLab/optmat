# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Gordon et al. (1993): Ar and H2 n - 1 at 1107-1150 A (Tables I, II).

Both tables have merged cells and are read from the JSON cells. Table I: the library
joins the Ar row's experiment cell of 1109.45 A with the H2 row's "see Table II" (the
image shows ``8.2+-0.5`` for Ar and "see Table II" in the H2 row below, which the
library dropped); the wavelength and transition cells span the Ar and H2 rows. Table II:
the image prints the footnote mark ``a`` ("not observed") in the rows of R(4) and R(2),
which have no value; the library attaches it to the values above them (R(5), R(3)),
where it is stripped.
"""

import json

from _library import LIBRARY, Point, number, write

DIRECTORY = "Gordon_1993_UseCoherentPhase"
TABLE_I = "tbl_d009cc88bfad"
TABLE_II = "tbl_c084e004c082"


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the sigma of a ``7.5±0.4`` cell."""
    value, _, sigma = cell.partition("±")
    return number(value), sigma


def extract() -> dict[str, list[Point]]:
    """Return Table I's experiment and theory per gas, and Table II's H2 values."""
    files: dict[str, list[Point]] = {}
    wavelength = ""
    for cells in grid(TABLE_I)[1:]:
        wavelength = cells[0] or wavelength
        gas = cells[2]
        printed = " | ".join(cells)
        value, sigma = _split(cells[4].removesuffix("see Table II"))
        files.setdefault(f"{gas}/gordon_1993__experiment_{gas}.csv", []).append(
            Point(wavelength, wavelength, value, sigma, printed)
        )
        files.setdefault(f"{gas}/gordon_1993__theory_{gas}.csv", []).append(
            Point(wavelength, wavelength, number(cells[5]), "", printed)
        )
    anomalous = []
    for cells in grid(TABLE_II)[1:]:
        if not cells[3]:
            continue  # R(4), R(2): "not observed" (footnote a)
        value, sigma = _split(cells[3])
        anomalous.append(
            Point(cells[0], cells[0], value, sigma.rstrip("a"), " | ".join(cells))
        )
    files["H2/gordon_1993__anomalous_H2.csv"] = anomalous
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

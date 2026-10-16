# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Tang & Poon (1974): continued-factorization bounds, Ne-Xe (Tables I-VI, IX).

Table I: Cauchy moments S(-2)..S(-8) per gas (copied from Langhoff and Karplus, except
the second xenon row, Tang and Poon's recommended set, labelled ``(2)^c`` in the library
because the gas cell is merged). Tables II-V: per gas, the columns of bounds and means
(n - 1 x 10^4, NTP), cells like ``0.7649±0.0002`` split into value and spread. Table VI
(merged header; JSON cells): xenon from the recommended moments beside Kingston's,
Starkschall and Gordon's, Abjean et al.'s and Chashchina and Shreider's values; empty
cells are no value, and Chashchina and Shreider's 1500 A value is not extracted: their
xenon measurements stop at 1550 A (Leonard 1974, p. 4; Bideau-Mehu et al. 1981, p. 395),
so it is an extrapolation of their formula. Table IX: the 18 Cauchy moments S(-2n) of
argon from Chashchina et al.'s Sellmeier formula, printed as ``1.1044 \times 10``
(stored as ``1.1044e1``).
"""

import json
import re

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Tang_1974_ContinuedFactorizationMethod"
TABLE_I = "tbl_d9080d4b7bd8"
TABLES_II_V = {
    "Ne": ("tbl_59039d75fc92", ("upper41", "pade", "mean", "saxon_1973")),
    "Ar": ("tbl_6787ce6f0228", ("upper41", "pade", "mean", "kingston_1964")),
    "Kr": ("tbl_fa3fc33a90fd", ("upper41", "pade", "mean", "kingston_1964")),
    "Xe": ("tbl_cc7175832502", ("upper41", "pade", "lower42", "mean", "kingston_1964")),
}
TABLE_VI = "tbl_f64c046e444d"
TABLE_VI_COLUMNS = (
    "mean_recommended",
    "kingston_1964_table6",
    "starkschall_gordon_1971",
    "abjean_1971",
    "chashchina_shreider_1969",
)
TABLE_IX = "tbl_e64e93bfd6da"
_TIMES_TEN = re.compile(r"^([\d.]+) \\times 10(?:\^\{(\d+)\})?$")


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the spread of a ``9.093± 0.026`` cell."""
    value, _, spread = cell.partition("±")
    return number(value), spread.strip()


def _times_ten(cell: str) -> str:
    r"""Return ``1.1044 \times 10`` as ``1.1044e1`` and ``... 10^{2}`` as ``...e2``."""
    match = _TIMES_TEN.match(cell)
    if match is None:
        msg = f"not a mantissa x 10^k cell: {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2) or '1'}"


def extract() -> dict[str, list[Point]]:
    """Return the moments, the bound columns and the reproduced columns."""
    files: dict[str, list[Point]] = {}

    def add(gas: str, name: str, point: Point) -> None:
        files.setdefault(f"{gas}/tang_1974__{name}.csv", []).append(point)

    for cells in table(DIRECTORY, TABLE_I)[1:]:
        # "Xe (1)" is Langhoff and Karplus's set, "(2)^c" the recommended one.
        name = "cauchy_recommended_Xe" if cells[0].startswith("(2)") else None
        gas = "Xe" if name else cells[0].split()[0]
        for order, cell in zip(("-2", "-4", "-6", "-8"), cells[1:5], strict=True):
            point = Point(order, order, number(cell), "", " | ".join(cells))
            add(gas, name or f"cauchy_{gas}", point)

    for gas, (table_id, columns) in TABLES_II_V.items():
        for cells in table(DIRECTORY, table_id)[1:]:
            for name, cell in zip(columns, cells[1:], strict=True):
                value, spread = _split(cell)
                point = Point(cells[0], cells[0], value, spread, " | ".join(cells))
                add(gas, f"{name}_{gas}", point)

    for cells in grid(TABLE_VI)[2:]:
        for name, cell in zip(TABLE_VI_COLUMNS, cells[1:], strict=True):
            if not cell or (name == "chashchina_shreider_1969" and cells[0] == "1500"):
                continue
            value, spread = _split(cell)
            add(
                "Xe",
                f"{name}_Xe",
                Point(cells[0], cells[0], value, spread, " | ".join(cells)),
            )

    for cells in table(DIRECTORY, TABLE_IX)[1:]:
        order = f"-{2 * int(cells[0])}"
        add(
            "Ar",
            "cauchy_chashchina_Ar",
            Point(order, order, _times_ten(cells[1]), "", " | ".join(cells)),
        )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

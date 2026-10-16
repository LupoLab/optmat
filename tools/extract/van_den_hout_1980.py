# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract van den Hout et al. (1980), Physica 104A, 509: H2, D2 S-line widths, shifts.

Three tables exist only as JSON cells; each was checked against its image. Each is in
sections headed by a row naming the gas (the header row's other cells empty); a data
row's first cell is the temperature in K (``293^b``: the room-temperature values taken
from Keijser et al. 1974, note b). "-" or an empty cell is a line not given. HD's
sections are not taken (the record's notes).

* Table I (``tbl_51fa3e7377dd``, p. 526): HWHM broadening coefficients
  Delta nu_1/2 / rho, MHz m3/mol, "v ± s", of normal H2 and normal D2.
* Table II (``tbl_b01b59957931``, p. 528): the same, extrapolated to pure oH2, pH2, pD2
  and oD2 (the section header sits in the S0(2) column).
* Table V (``tbl_cd5b5a338e81``, p. 533): line-shift coefficients Delta nu_S / rho, MHz
  m3/mol, without uncertainties except the 293 K row (Cooper et al., note b). Values
  printed in parentheses are written without them; the record excludes those rows.

One file per table, gas and temperature; rows are the lines, "S0(J)".
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "VanDenHout_1980_BroadeningShiftRatational"
TABLE_I = "tbl_51fa3e7377dd"
TABLE_II = "tbl_b01b59957931"
TABLE_V = "tbl_cd5b5a338e81"
# Section headers as the cells print them, and the (material, file tag) they head.
SECTIONS = {
    "nH_2 (25.1\\% pH_2)": ("H2", "H2"),
    "nD_2 (66.7\\% oD_2)": ("D2", "D2"),
    "nH_2": ("H2", "H2"),
    "nD_2": ("D2", "D2"),
    "oH_2": ("H2", "oH2"),
    "pH_2": ("H2", "pH2"),
    "pD_2": ("D2", "pD2"),
    "oD_2": ("D2", "oD2"),
    "HD": ("HD", "HD"),
}
ABSENT = {"", "-"}


def _temperature(cell: str) -> str:
    """Return the temperature of a row, without the note mark (``293^b``)."""
    return number(cell.removesuffix("^b"))


def _tag(temperature: str) -> str:
    return temperature.replace(".", "p") + "K"


def _split(cell: str) -> tuple[str, str]:
    """Return (value, sigma) of a "v ± s", "v" or "(v)" cell."""
    value, _, sigma = cell.replace(" ", "").partition("±")
    value = value.removeprefix("(").removesuffix(")")
    return number(value), number(sigma) if sigma else ""


def _sections(
    table_id: str, lines: list[str], first: int
) -> dict[tuple[str, str, str], list[Point]]:
    """Return {(material, tag, temperature): points} for the H2 and D2 sections."""
    grid = json_table(DIRECTORY, table_id)
    found: dict[tuple[str, str, str], list[Point]] = {}
    section = None
    for cells in grid[first:]:
        header = [cell for cell in cells if cell]
        if len(header) == 1 and header[0] in SECTIONS:
            section = SECTIONS[header[0]]
            continue
        if section is None:
            msg = f"{table_id}: a data row before any section: {cells!r}"
            raise ValueError(msg)
        material, tag = section
        if material == "HD":
            continue  # not a material: the record's notes
        temperature = _temperature(cells[0])
        printed = " | ".join(cells)
        for line, cell in zip(lines, cells[1:], strict=False):
            if cell in ABSENT:
                continue
            value, sigma = _split(cell)
            found.setdefault((material, tag, temperature), []).append(
                Point(line, line, value, sigma, printed)
            )
    return found


def _files(
    stem: str, found: dict[tuple[str, str, str], list[Point]]
) -> dict[str, list[Point]]:
    return {
        f"{material}/van_den_hout_1980__{stem}_{tag}_{_tag(t)}.csv": points
        for (material, tag, t), points in found.items()
    }


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Tables I, II and V."""
    lines = [f"S0({j})" for j in range(6)]
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0] != ["T(K)", *(f"S_0({j})" for j in range(6))]:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    grid = json_table(DIRECTORY, TABLE_II)
    if "MHz m^3/mol" not in grid[0][1] or grid[1][1:] != lines[:5]:
        msg = f"unexpected Table II header {grid[:2]!r}"
        raise ValueError(msg)
    grid = json_table(DIRECTORY, TABLE_V)
    if grid[0] != ["T(K)", *(f"S_0({j})" for j in range(5))]:
        msg = f"unexpected Table V header {grid[0]!r}"
        raise ValueError(msg)
    return (
        _files("table1", _sections(TABLE_I, lines, 1))
        | _files("table2", _sections(TABLE_II, lines, 2))
        | _files("table5", _sections(TABLE_V, lines, 1))
    )


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

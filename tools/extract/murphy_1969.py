# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Murphy, Holzer and Bernstein (1969): gas-phase Raman intensities (Table I).

Table I is printed over five pages; the library splits it into ``tbl_ca70dea8e569``
(p. 212, the first part), ``tbl_30a71dee5856`` (p. 213), ``tbl_e29dfa66b07c`` (p. 214),
``tbl_d8cdb7b2c9ab`` (p. 215) and ``tbl_92207413a7ff`` (p. 216), in that order (the
page order of the PDF; each continuation repeats the header). Every curated molecule
(H2, D2, O2, N2, CO, CO2, CH4, SF6) is in the first part; the continuations hold only
halomethanes, hydrocarbons and other molecules outside the scope, and are not read.

STRUCTURAL FIX. The library has no CSV for the first part, and its cells merge the
sub-lines of each printed row without separators (``0.0700.073`` for H2's two
depolarization ratios, ``2.082.062.152.10`` for methane's four derivatives), attach the
reference numbers of all sub-lines to one cell, and lose some sub-lines altogether: the
third and fourth lines of CO2's I_j/I_1 and rho_nj columns (``[0.18Q]``,
``[0.019Q]``; ``0.054Q``, ``0.15Q``) and methane's two scattering activities of the
2914 cm-1 band (``207``, ref. 8; ``200``, ref. 19). The in-scope lines are therefore
transcribed below from the table image (``tbl_ca70dea8e569.jpg``, read at twice its
size), each with the reference printed on its line, and the script checks that the
transcribed sub-lines of every cell, concatenated, are exactly the library cell (with
the column factors removed); the lost lines are listed with ``cell = None``.

What the image shows and the cells lose or misread:
- the column factors are printed once, in the H2 row: ``[147]·N·10^-32`` (the
  scattering activity, cm^4 g^-1), ``[1.73]·N^(1/2)·10^-16`` (dalpha/dQ, cm^2 g^-1/2;
  the cells read the exponent as ``\frac{1}{4}``, the image and the dimensions give
  1/2) and ``1.23·10^-16`` (dalpha/dr, cm^2). They apply to the whole column. N is
  Avogadro's number, so the activity is in angstrom^4/amu and dalpha/dQ in
  angstrom^2/amu^(1/2) exactly (the source record shows this against the paper's own
  cross section);
- brackets mark 'calculated from original literature data' (by these authors) and are
  kept in ``printed``; a ``Q`` after a value means the Q branch only; the reference
  number carries a superscript ``n`` ('original data obtained using unpolarized
  incident irradiation');
- the CO2 activities are braced to the single dalpha/dQ and dalpha/dr, derived from
  their sum;
- methane's four derivative lines of the 2914 cm-1 band belong to refs. 17, 18, 8 and
  19 in that order (one per line), the activities 207 and 200 to the lines of refs. 8
  and 19.

Not transcribed (in the record's notes): the depolarization ratios printed ``6/7``
(the theoretical value of a depolarized band, not a number measured), methane's
``[<0.006]`` (an upper limit) and the ``Γ`` entries (a pointer to the references).
"""

from dataclasses import dataclass

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Murphy_1969_GasPhaseRaman"
TABLE_I = "tbl_ca70dea8e569"
EXCITATION = "4358"  # angstrom (the mercury line; p. 217), the abscissa of every file
# The column factors printed once, in the H2 row (see the module docstring).
FACTORS = (r"\cdot N\cdot 10^{-32}", r"\cdot N^{\frac{1}{4}}\cdot 10^{-16}")
FACTORS += (r"\cdot 10^{-16}",)
# Datasets whose abscissa is a label (the row); the others are at 4358 A.
LABELLED = ("shift_", "activity_")
# The library row of each curated molecule's first line, and its printed name.
GAS_ROWS = {1: "H_{2}", 4: "D_{2}", 5: "O_{2}", 6: "N_{2}", 7: "CO", 8: "CO_{2}"}
GAS_ROWS |= {15: "CH_{4}", 49: "SF_{6}"}


@dataclass(frozen=True)
class Line:
    """One printed sub-line of Table I, as the image shows it."""

    dataset: str
    row: str
    printed: str
    ref: str
    cell: tuple[int, int] | None
    """The (row, column) of the library cell that holds it; None where it is lost."""


def _lines() -> list[Line]:  # noqa: PLR0915 - one statement per transcribed line
    """Return the transcription of the in-scope lines of Table I (image checked)."""
    table: list[Line] = []

    def add(
        dataset: str, row: str, printed: str, ref: str, cell: tuple[int, int] | None
    ) -> None:
        table.append(Line(dataset, row, printed, ref, cell))

    # H2 (library row 1).
    add("shift_H2", "v=0-1", "4156", "1, 10", (1, 1))
    add("ratio_H2", "band (ref 1)", "[0.40]", "1", (1, 2))
    add("rho_H2", "band (ref 1)", "0.070", "1", (1, 3))
    add("rho_H2", "band (ref 10)", "0.073", "10", (1, 3))
    add("activity_H2", "band (ref 1)", "[147]", "1", (1, 4))
    add("dadQ_H2", "band (ref 1)", "[1.73]", "1", (1, 5))
    add("dadr_H2", "band (ref 1)", "1.23", "1", (1, 6))
    add("dadr_H2", "band (ref 10)", "1.3", "10", (1, 6))
    # D2 (row 4): the activity and derivatives are on the second line (ref. 10).
    add("shift_D2", "v=0-1", "2986", "1, 10", (4, 1))
    add("ratio_D2", "Q branch (ref 1)", "[0.32Q]", "1", (4, 2))
    add("ratio_D2", "band (ref 10)", "[0.35]", "10", (4, 2))
    add("rho_D2", "band (ref 1)", "0.07", "1", (4, 3))
    add("rho_D2", "band (ref 10)", "0.073", "10", (4, 3))
    add("activity_D2", "band (ref 10)", "[73]", "10", (4, 4))
    add("dadQ_D2", "band (ref 10)", "[1.22]", "10", (4, 5))
    add("dadr_D2", "band (ref 10)", "1.22", "10", (4, 6))
    # O2, N2, CO (rows 5, 6, 7): first line the Q branch, second the whole band.
    for gas, row, shift, q_ratio, ratio, q_rho, rho, activity, dq, dr in (
        (
            "O2",
            5,
            "1555",
            "[0.17Q]",
            "[0.24]",
            "0.115Q",
            "0.33",
            "[19.7]",
            "[0.52]",
            "[1.46]",
        ),
        (
            "N2",
            6,
            "2331",
            "[0.14Q]",
            "[0.17]",
            "0.055Q",
            "0.18",
            "[24.5]",
            "[0.66]",
            "[1.75]",
        ),
        (
            "CO",
            7,
            "2143",
            "[0.13Q]",
            "[0.165]",
            "0.104Q",
            "",
            "[23.1]",
            "[0.58]",
            "[1.50]",
        ),
    ):
        add(f"shift_{gas}", "v=0-1", shift, "1", (row, 1))
        add(f"ratio_{gas}", "Q branch (ref 1)", q_ratio, "1", (row, 2))
        add(f"ratio_{gas}", "band (ref 1)", ratio, "1", (row, 2))
        add(f"rho_{gas}", "Q branch (ref 1)", q_rho, "1", (row, 3))
        if rho:
            add(f"rho_{gas}", "band (ref 1)", rho, "1", (row, 3))
        add(f"activity_{gas}", "band (ref 1)", activity, "1", (row, 4))
        add(f"dadQ_{gas}", "band (ref 1)", dq, "1", (row, 5))
        add(f"dadr_{gas}", "band (ref 1)", dr, "1", (row, 6))
    # CO2 (rows 8-11): four Q branches; the braced derivatives come from their sum.
    for index, (shift, ratio, rho, activity) in enumerate(
        (
            ("1265", "[0.012Q]", "0.25Q", "[1.45]"),
            ("1285", "[0.125Q]", "0.084Q", "[10.5]"),
            ("1388", "[0.18Q]", "0.054Q", "[15.5]"),
            ("1409", "[0.019Q]", "0.15Q", "[2.15]"),
        )
    ):
        lost = index >= 2  # the cells keep only the first two lines of these columns
        add("shift_CO2", shift, shift, "1", (8 + index, 1))
        add(
            "ratio_CO2",
            f"{shift} Q branch (ref 1)",
            ratio,
            "1",
            None if lost else (8, 2),
        )
        add("rho_CO2", f"{shift} Q branch (ref 1)", rho, "1", None if lost else (8, 3))
        add("activity_CO2", shift, activity, "1", (8, 4))
    add("dadQ_CO2", "sum of the four bands (ref 1)", "[0.69]", "1", (8, 5))
    add("dadr_CO2", "sum of the four bands (ref 1)", "[1.95]", "1", (8, 6))
    # CH4 (rows 15-21).
    add("shift_CH4", "1306", "1306", "16", (15, 1))
    add("shift_CH4", "1535 (nu2)", "1535", "16, 18", (16, 1))
    add("shift_CH4", "2914 (nu1, A1)", "2914", "17, 18, 8, 19", (17, 1))
    add("shift_CH4", "3017", "3017", "16, 18", (21, 1))
    add("ratio_CH4_1535", "band (ref 16)", "[0.115]", "16", (16, 2))
    add("ratio_CH4_1535", "Q branch (ref 18)", "[0.031Q]", "18", (16, 2))
    add("activity_CH4_1535", "band (ref 16)", "[9.1]", "16", (16, 4))
    add("rho_CH4_nu1", "band (ref 17)", "0", "17", (17, 3))
    for ref, dq, dr in (("17", "2.08", "1.04"), ("18", "2.06", "1.03"),
                        ("8", "2.15", "1.08"), ("19", "2.10", "1.08")):  # fmt: skip
        add("dadQ_CH4_nu1", f"band (ref {ref})", dq, ref, (17, 5))
        add("dadr_CH4_nu1", f"band (ref {ref})", dr, ref, (17, 6))
    add("activity_CH4_nu1", "band (ref 8)", "207", "8", None)
    add("activity_CH4_nu1", "band (ref 19)", "200", "19", None)
    add("ratio_CH4_3017", "band (ref 16)", "[0.68]", "16", (21, 2))
    add("ratio_CH4_3017", "Q branch (ref 18)", "[0.26Q]", "18", (21, 2))
    add("activity_CH4_3017", "band (ref 16)", "[140]", "16", (21, 4))
    # SF6 (row 49).
    add("shift_SF6", "775", "775", "20", (49, 1))
    add("ratio_SF6", "band (ref 20)", "0.65", "20", (49, 2))
    add("rho_SF6", "band (ref 20)", "0", "20", (49, 3))
    add("activity_SF6", "band (ref 20)", "[22.7]", "20", (49, 4))
    add("dadQ_SF6", "band (ref 20)", "[0.71]", "20", (49, 5))
    add("dadr_SF6", "band (ref 20)", "1.23", "20", (49, 6))
    return table


def _bare(cell: str) -> str:
    """Return a library cell without the column factors printed in the H2 row."""
    for factor in FACTORS:
        cell = cell.replace(factor, "")
    return cell.strip()


def _check(lines: list[Line], grid: list[list[str]]) -> None:
    """Check that each library cell is the concatenation of its transcribed lines."""
    cells: dict[tuple[int, int], str] = {}
    for line in lines:
        if line.cell is not None:
            cells[line.cell] = cells.get(line.cell, "") + line.printed
    for (row, column), joined in cells.items():
        if _bare(grid[row][column]) != joined:
            msg = f"cell {row, column} is {grid[row][column]!r}, transcribed {joined!r}"
            raise ValueError(msg)
    for row, gas in GAS_ROWS.items():
        if grid[row][0] != gas:
            msg = f"row {row} is {grid[row][0]!r}, not {gas!r}"
            raise ValueError(msg)


def _value(printed: str) -> str:
    """Return the number of a printed entry, without its brackets and Q mark."""
    return number(printed.strip("[]").removesuffix("Q"))


def extract() -> dict[str, list[Point]]:
    """Return the in-scope lines of Table I as data files."""
    lines = _lines()
    _check(lines, json_table(DIRECTORY, TABLE_I))
    files: dict[str, list[Point]] = {}
    for line in lines:
        gas = line.dataset.split("_")[1]
        x = line.row if line.dataset.startswith(LABELLED) else EXCITATION
        point = Point(
            line.row, x, _value(line.printed), "", f"{line.printed} | ref. {line.ref}"
        )
        files.setdefault(f"{gas}/murphy_1969__{line.dataset}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Millot et al. (1992): O2 Q-branch widths and shifts, self and by N2.

J. Chem. Phys. 96, 961. All tables checked against their images. Widths are HWHM and
shifts per atm, in 10^-3 cm-1 atm-1 (the mean shifts per amagat).

* Table II (``tbl_23bcbc0cdebf``, CSV): O2 self-broadening against N at 295-1350 K,
  with Ouazzany et al.'s CARS values (note a), the pure anisotropic widths (b), Berard
  et al.'s pure rotational widths (c) and the first and second hot bands at 1350 K (d,
  e). The note marks are extracted as ``295\mathrm {K}^{\mathrm {a}}``; the 1350 K
  (fundamental) column holds only N = 39-55 and the 1350 K^d column N = 1-37, as the
  image shows.
* Table III (``tbl_7b0c5f5de7ea``, CSV): O2 broadened by N2.
* Table IV (``tbl_84280b026de1``, cells only): line shifts, O2-O2 and O2-N2, and their
  mean values per amagat (last row). Its caption is extracted as a table note and the
  running head ("Millot et al.: Collisional effects in Raman scattering") as the
  caption; the image and the text (p. 7) identify it as Table IV.
* Tables VI and VII (``tbl_6eda03db14cb``, ``tbl_6dcaeaf58dd1``, CSV): the
  semiclassical widths (VI: Q(N) and, at 295 K, S(N) of pure O2; VII: O2-N2).
* Table VIII (``tbl_9dc671631d22``, cells only): the semiclassical shifts, (A) pure O2
  and (B) O2-N2, the temperature cell spanning both rows.

Cells are "v ± s" (or "v \pm s"); empty cells are no value.
"""

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Millot_1992_CollisionalEffectsStimulated"
TABLE_II = "tbl_23bcbc0cdebf"
TABLE_III = "tbl_7b0c5f5de7ea"
TABLE_IV = "tbl_84280b026de1"
TABLE_VI = "tbl_6eda03db14cb"
TABLE_VII = "tbl_6dcaeaf58dd1"
TABLE_VIII = "tbl_9dc671631d22"


def _split(cell: str) -> tuple[str, str]:
    r"""Return (value, sigma) of a "v ± s", "v \pm s" or "v" cell."""
    value, _, sigma = cell.replace("\\pm", "±").replace(" ", "").partition("±")
    return number(value), number(sigma) if sigma else ""


def _columns(
    grid: list[list[str]], first: int, names: dict[str, str]
) -> dict[str, list[Point]]:
    """Return one file per named header column of N rows (row label in column 0)."""
    header = grid[0]
    if set(names) - set(header):
        msg = f"missing columns {set(names) - set(header)!r} in {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[first:]:
        row = f"N={number(cells[0])}"
        for head, name in names.items():
            cell = cells[header.index(head)]
            if not cell:
                continue
            value, sigma = _split(cell)
            files.setdefault(f"O2/millot_1992__{name}.csv", []).append(
                Point(row, row, value, sigma, " | ".join(cells))
            )
    return files


def _table_iv() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_IV)
    temperatures = ["295", "295", "446", "446", "990", "1317"]
    perturbers = ["O2", "N2", "O2", "N2", "O2", "O2"]
    heads = [cell.replace("_{2}", "2") for cell in grid[1][1:]]
    if [h.split("-")[1] for h in heads] != perturbers or grid[0][1::2][:3] != [
        "295 K",
        "446 K",
        "990 K",
    ]:
        msg = f"unexpected Table IV header {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        mean = cells[0].startswith("Mean value")
        for column, (temperature, perturber) in enumerate(
            zip(temperatures, perturbers, strict=True), start=1
        ):
            if not cells[column]:
                continue
            value, sigma = _split(cells[column])
            self_ = "self" if perturber == "O2" else "N2"
            if mean:
                name = f"table4_mean_{self_}"
                point = Point(
                    f"{temperature}K", temperature, value, sigma, " | ".join(cells)
                )
            else:
                name = f"table4_{self_}_{temperature}K"
                row = f"N={number(cells[0])}"
                point = Point(row, row, value, sigma, " | ".join(cells))
            files.setdefault(f"O2/millot_1992__{name}.csv", []).append(point)
    return files


def _table_viii() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_VIII)
    lines = [cell.replace(" ", "") for cell in grid[0][2:]]
    if lines != ["N=1", "N=5", "N=9", "N=13", "N=17", "N=21"]:
        msg = f"unexpected Table VIII header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    temperature = ""
    for cells in grid[1:]:
        temperature = cells[0].removesuffix(" K") or temperature
        system = {"A": "self", "B": "N2"}[cells[1]]
        files[f"O2/millot_1992__table8_{system}_{temperature}K.csv"] = [
            Point(line, line, number(cell), "", " | ".join(cells))
            for line, cell in zip(lines, cells[2:], strict=True)
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Table II, III, IV, VI, VII and VIII files."""
    files = _columns(
        table(DIRECTORY, TABLE_II),
        1,
        {
            "295 K": "table2_295K",
            "295\\mathrm {K}^{\\mathrm {a}}": "table2_ouazzany_295K",
            "295\\mathrm {K}^{\\mathrm {b}}": "table2_anisotropic_295K",
            "293\\mathrm {K}^{\\mathrm {c}}": "table2_berard_293K",
            "446 K": "table2_446K",
            "769 K": "table2_769K",
            "990 K": "table2_990K",
            "1317 K": "table2_1317K",
            "1350 K": "table2_1350K",
            "1350\\mathrm {K}^{\\mathrm {d}}": "table2_hot1_1350K",
            "1350\\mathrm {K}^{\\mathrm {e}}": "table2_hot2_1350K",
        },
    )
    files |= _columns(
        table(DIRECTORY, TABLE_III),
        1,
        {
            "295 K": "table3_N2_295K",
            "446 K": "table3_N2_446K",
            "990 K": "table3_N2_990K",
        },
    )
    files |= _table_iv()
    files |= _columns(
        table(DIRECTORY, TABLE_VI),
        1,
        {
            "295\\mathrm {K}^{\\mathrm {a}}": "table6_Q_295K",
            "295\\mathrm {K}^{\\mathrm {b}}": "table6_S_295K",
            "446 K": "table6_Q_446K",
            "769 K": "table6_Q_769K",
            "990 K": "table6_Q_990K",
            "1317 K": "table6_Q_1317K",
            "1350 K": "table6_Q_1350K",
        },
    )
    files |= _columns(
        table(DIRECTORY, TABLE_VII),
        1,
        {
            "295 K": "table7_N2_295K",
            "446 K": "table7_N2_446K",
            "990 K": "table7_N2_990K",
            "1315 K": "table7_N2_1315K",
            "2000 K": "table7_N2_2000K",
        },
    )
    return files | _table_viii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

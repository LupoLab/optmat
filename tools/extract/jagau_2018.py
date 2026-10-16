# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Jagau (2018), supplement: all-electron static ionisation rates of molecules.

The supplementary material (``supplement_01``) prints, against the field F (au), the
Stark shift Delta E and the rate Gamma (au) of the ground states of H2, CO, CH4 and O2
at the HF and correlated (CCSD; full CI for H2) levels, from complex-scaled Gaussian
basis functions. Taken, the rates only (the Stark shifts stay in the ``printed``
cells where they share a row):

* Table S3 (``tbl_f34ab4e8be6f``), H2 at R(HH) = 0.74 A, parallel and perpendicular,
  the full CI and HF columns of the standard basis aug-cc-pVQZ+3s3p3d3f (Table S1);
  the two other bases are the paper's convergence study;
* Table S6 (``tbl_aef780ef04ac``), CO at 13 angles phi (0 degrees: the field from C to
  O, the reverse of majety_2015's), at F = 0.06 and 0.09;
* Table S8 (``tbl_b453fabca2b8``), CH4 in three orientations A, B, C;
* Table S9 (``tbl_13bba2b740f3``), O2 parallel (0 degrees) and perpendicular (90);
* Table S10 (``tbl_9352914b0c27``), O2 at the angles between, at F = 0.06 and 0.10 (its
  0 and 90 degree rows are Table S9's and are not taken twice).

One file per molecule, level and angle or orientation. A cell printed as a dash (no
resonance at that level) is skipped. The tables have merged header cells, so they are
read from the ``.json`` cells; the script checks each header and each row label. Every
rate cell was checked against the table images (2026-09-29): they agree.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Jagau_2018_CoupledClusterTreatment/supplement_01"
DASH = "\N{EM DASH}"
DEGREE = "\N{DEGREE SIGN}"
PARALLEL, PERPENDICULAR = "(||)", "(\N{UP TACK})"


def _rate(cell: str) -> str | None:
    """Return a rate cell's number as printed, or None for a dash."""
    return None if cell == DASH else number(cell)


def _h2() -> dict[str, list[Point]]:
    """Return Table S3's full CI and HF rates of H2, parallel and perpendicular."""
    grid = json_table(DIRECTORY, "tbl_f34ab4e8be6f")
    if (
        grid[0][1] != "Full CI"
        or grid[0][4] != "HF"
        or grid[13][0] != ("Ionization Rates Γ")
    ):
        msg = f"unexpected layout of Table S3: {grid[0]!r}, {grid[13]!r}"
        raise ValueError(msg)
    if grid[1][2] != "aug-cc-pVQZ+3s3p3d3f" or grid[1][5] != "aug-cc-pVQZ+3s3p3d3f":
        msg = f"unexpected basis columns of Table S3: {grid[1]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[14:24]:
        field, direction = cells[0].split(" ")
        angle = {PARALLEL: "0", PERPENDICULAR: "90"}[direction]
        for level, column in (("fci", 2), ("hf", 5)):
            rate = _rate(cells[column])
            if rate is None:
                continue
            files.setdefault(
                f"H2/jagau_2018__tableS3_{level}_H2_{angle}deg.csv", []
            ).append(Point(number(field), number(field), rate, "", " | ".join(cells)))
    return files


def _co() -> dict[str, list[Point]]:
    """Return Table S6's HF and CCSD rates of CO at 13 angles and two fields."""
    grid = json_table(DIRECTORY, "tbl_aef780ef04ac")
    if (
        grid[0][1] != "F=0.06"
        or grid[0][3] != "F=0.09"
        or grid[16][0] != ("Ionization Rates \\Gamma")
    ):
        msg = f"unexpected layout of Table S6: {grid[0]!r}, {grid[16]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[17:30]:
        angle = cells[0].removesuffix(DEGREE)
        for level, columns in (("hf", (1, 3)), ("ccsd", (2, 4))):
            points = [
                Point(field, field, number(cells[column]), "", " | ".join(cells))
                for field, column in zip(("0.06", "0.09"), columns, strict=True)
            ]
            files[f"CO/jagau_2018__tableS6_{level}_CO_{angle}deg.csv"] = points
    return files


def _ch4() -> dict[str, list[Point]]:
    """Return Table S8's HF and CCSD rates of CH4 in orientations A, B and C."""
    grid = json_table(DIRECTORY, "tbl_b453fabca2b8")
    if grid[1] != ["", "HF", "CCSD", "HF", "CCSD"]:
        msg = f"unexpected header of Table S8: {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:17]:
        field, orientation = cells[0][:-3], cells[0][-2]
        if cells[0][-3:] != f"({orientation})" or orientation not in "ABC":
            msg = f"unexpected row label in Table S8: {cells[0]!r}"
            raise ValueError(msg)
        for level, column in (("hf", 3), ("ccsd", 4)):
            files.setdefault(
                f"CH4/jagau_2018__tableS8_{level}_CH4_{orientation}.csv", []
            ).append(
                Point(
                    number(field),
                    number(field),
                    number(cells[column]),
                    "",
                    " | ".join(cells),
                )
            )
    return files


def _o2() -> dict[str, list[Point]]:
    """Return Tables S9 and S10's HF and CCSD rates of O2."""
    files: dict[str, list[Point]] = {}
    grid = json_table(DIRECTORY, "tbl_13bba2b740f3")
    if grid[1] != ["", *["HF", "CCSD"] * 4]:
        msg = f"unexpected header of Table S9: {grid[:2]!r}"
        raise ValueError(msg)
    for cells in grid[2:9]:
        field = number(cells[0])
        for angle, level, column in (
            ("0", "hf", 5),
            ("0", "ccsd", 6),
            ("90", "hf", 7),
            ("90", "ccsd", 8),
        ):
            rate = _rate(cells[column])
            if rate is None:
                continue
            files.setdefault(
                f"O2/jagau_2018__tableS9_{level}_O2_{angle}deg.csv", []
            ).append(Point(field, field, rate, "", " | ".join(cells)))
    grid = json_table(DIRECTORY, "tbl_9352914b0c27")
    if grid[0][5] != "\\Gamma(F=0.06)" or grid[0][7] != "\\Gamma(F=0.10)":
        msg = f"unexpected header of Table S10: {grid[0]!r}"
        raise ValueError(msg)
    for cells in grid[2:11]:
        angle = cells[0]
        if angle in {"0", "90"}:
            continue  # Table S9's rows
        for level, columns in (("hf", (5, 7)), ("ccsd", (6, 8))):
            files[f"O2/jagau_2018__tableS10_{level}_O2_{angle}deg.csv"] = [
                Point(field, field, number(cells[column]), "", " | ".join(cells))
                for field, column in zip(("0.06", "0.10"), columns, strict=True)
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the rates of Tables S3, S6, S8, S9 and S10, one file per series."""
    return _h2() | _co() | _ch4() | _o2()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Danichkin et al. (1981): a survey of vibrational Q-branch Raman parameters.

Table 1 (p. 1058, ``tbl_fede8ab1f190``; continued p. 1059, ``tbl_e069277a54e8``, and
p. 1060, ``tbl_c9af12ff6f78``) lists per molecule and transition the frequency shift,
the depolarization ratios measured by others (with the reference of each), the value
the survey adopts, and the cross sections relative to the N2 Q branch, each with the
excitation wavelength of the experiment and its reference, all 'referred to a single
system of measurements, corresponding to lambda_0 = 488.0 nm, linearly polarized
exciting radiation, observation angle theta = pi/2, scattering cross section equal to
the sum of both polarization components for the Q branch, unless otherwise stated'
(p. 1061). Table 2 (the library's unlabelled ``tbl_d60a413f5ef7``, p. 1062; its caption
is printed above it on the page, image checked) gives the absolute N2 cross sections
behind the ratios; Table 3 (``tbl_a78a07d029aa``, p. 1063) the survey's recommended
values.

STRUCTURAL FIX (Table 1). The library has no CSV for Table 1, and its cells join the
stacked lines of each column without separators (``0,0280,020±0,0010,0220,020...``),
join the references of a braced group into one, and lose lines at the foot of the
hydrogen block: the last two wavelengths (``488,0``, ``488,0``), the ratio
``1,6 (Q, J=1)`` and its reference ``[12]``. The in-scope lines are therefore
transcribed below from the table images (read at twice their size), each with its
wavelength and reference, and the script checks that the transcribed lines of each
column, concatenated, are exactly the library cell (after removing the words
'entire band' and the hyphenation 'en-tire band'); lost lines are marked ``lost``.

What the images show and the cells lose:
- 'entire band' under a ratio means the whole vibrational band, not the Q branch;
- braces join lines to one reference: ``[5, 10]`` (Yoshino and Bernstein as reduced
  by Murphy et al.), ``[32]`` (ozone), ``[18]`` (the two hydrogen lines of Carlson and
  Fenner), ``[12]`` (hydrogen's 2,44 and 1,6);
- hydrogen's ``2,44+0,24`` is printed with a plain ``+`` (the underline of the ``±``
  not visible), read as ``±`` like every other entry;
- carbon monoxide's first depolarization ratio is printed ``0.055`` with a point.

Numbers are printed with decimal commas; the value column has them as points, the
printed column keeps the print. Molecules outside the scope (F2, NO, H2O, HCl, SO2,
HBr, HF, Cl2, DCl, DBr, BrCl, Br2, I2, ICl, IBr, H2S) are not transcribed.
"""

from dataclasses import dataclass

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Danichkin_1981_RamanScatteringParameters"
PART1, PART2, PART3 = "tbl_fede8ab1f190", "tbl_e069277a54e8", "tbl_c9af12ff6f78"
TABLE_2 = "tbl_d60a413f5ef7"
TABLE_3 = "tbl_a78a07d029aa"
SYSTEM = "488.0"  # nm: the survey's reference excitation, the abscissa of every file
# Words the cells keep beside the numbers (removed before comparing).
WORDS = ("en-tire band", "entire band")


@dataclass(frozen=True)
class Ratio:
    """One line of the cross-section columns: wavelength, ratio, reference."""

    wavelength: str
    printed: str
    ref: str
    note: str = ""
    lost: bool = False
    """The whole line is lost from the library cells."""
    wavelength_lost: bool = False
    """Only its wavelength is lost from the library cells."""


@dataclass(frozen=True)
class Block:
    """One molecule and transition of Table 1, as the image shows it."""

    part: str
    row: int
    label: str
    material: str
    tag: str
    shifts: tuple[tuple[str, str], ...]
    rho: tuple[tuple[str, str], ...]
    ratios: tuple[Ratio, ...]


def R(  # noqa: N802
    wavelength: str, printed: str, ref: str, note: str = "", lost: str = ""
) -> Ratio:
    """Return one transcribed ratio line; ``lost`` is "", "wavelength" or "all"."""
    return Ratio(wavelength, printed, ref, note, lost == "all", lost == "wavelength")


BLOCKS = (
    Block(PART1, 2, "N_20—1", "N2", "N2", (("v=0-1", "2330"),),
          (("0,028", "5"), ("0,020±0,001", "6"), ("0,022", "7"), ("0,02", "8"),
           ("0,0211±0,0043", "9")),
          ()),
    Block(PART1, 3, "O_20—1", "O2", "O2", (("v=0-1", "1556"),),
          (("0,061", "5"), ("0,047", "7"), ("0,051±0,004", "6"), ("0,043±0,003", "9")),
          (R("514,5", "1,23±0,12", "11"), R("514,5", "1,23±0,06", "7"),
           R("488,0", "1,33±0,13", "12"), R("337,1", "1,23±0,06", "13"),
           R("488,0", "1,22±0,07", "14"), R("337,1", "1,30±0,15", "15"),
           R("337,1", "1,4±0,4", "37"), R("337,1", "1,4±0,14", "38"),
           R("435,8", "1,2±0,1", "5, 10"), R("435,8", "1,6", "5, 10", "entire band"))),
    Block(PART1, 4, 'D_20—1J\'→J"0—02—2 (max)', "D2", "D2",
          (("v=0-1 Q(0)", "2991"), ("v=0-1 Q(2) (maximum)", "2986")),
          (("0,0094", "5"), ("0,0099", "16"), ("0,0125±0,001", "6")),
          (R("435,8", "2,2", "5, 10"), R("435,8", "2,46", "5, 10", "entire band"))),
    Block(PART1, 5, 'H_20—1J\'→J"0—01—1 (max)', "H2", "H2",
          (("v=0-1 Q(0)", "4162"), ("v=0-1 Q(1) (maximum)", "4156")),
          (("0,0094", "5"), ("0,0099±0,0013", "16"), ("0,0112", "17"),
           ("0,0125±0,001", "6")),
          (R("514,5", "2,26", "12"), R("488,0", "4,21±0,84", "18"),
           R("514,5", "5,88±1,18", "18"), R("337,1", "1,62±0,21", "15"),
           R("337,1", "2,6±0,8", "37"), R("337,1", "2,8±0,3", "38"),
           R("435,8", "2,56", "5, 10", "entire band"),
           R("488,0", "2,44+0,24", "12", "entire band", lost="wavelength"),
           R("488,0", "1,6", "12", "Q(1) (printed 'Q, J=1')", lost="all"))),
    Block(PART1, 11, "O_3v_1", "O3", "O3", (("nu1", "1103"),),
          (("0,026±0,005", "32"), ("0,015±0,005", "32"), ("0,020±0,005", "32"),
           ("0,053±0,005", "32")),
          (R("514,5", "3,8", "22"), R("476,5", "2,2±0,3", "32"),
           R("488,0", "2,4±0,3", "32"), R("496,5", "2,7±0,3", "32"),
           R("514,5", "2,1±0,3", "32"))),
    Block(PART2, 2, "CO0-1", "CO", "CO", (("v=0-1", "2143"),),
          (("0.055", "5"), ("0,038", "7")),
          (R("514,5", "0,91", "11"), R("514,5", "1,01±0,05", "7"),
           R("337,1", "0,98", "13"), R("488,0", "1,02±0,10", "12"),
           R("488,0", "0,75±0,09", "14"), R("337,1", "1,0±0,3", "37"),
           R("337,1", "0,9±0,09", "38"), R("435,8", "0,93±0,09", "5, 10"),
           R("435,8", "1,15±0,11", "5, 10", "entire band"))),
    Block(PART2, 5, "CO_2v_1", "CO2", "CO2_nu1", (("nu1", "1388"),),
          (("0,028", "5"), ("0,027", "7")),
          (R("435,8", "1,43±0,14", "5, 10"), R("514,5", "1,61", "11"),
           R("514,5", "1,50±0,10", "7"), R("337,1", "1,44±0,07", "13"),
           R("488,0", "1,41±0,14", "12"), R("488,0", "1,65±0,11", "14"),
           R("337,1", "1,4±0,4", "37"), R("337,1", "1,20±0,12", "38"))),
    Block(PART2, 6, "CO_22v_2", "CO2", "CO2_2nu2", (("2nu2", "1285"),),
          (("0,044", "5"),),
          (R("435,8", "1,00±0,10", "5, 10"), R("514,5", "1,07", "11"),
           R("337,1", "1,08±0,05", "13"), R("488,0", "0,91±0,09", "12"),
           R("488,0", "1,15±0,11", "14"), R("337,1", "0,9±0,3", "37"),
           R("337,1", "0,84±0,08", "38"))),
    Block(PART2, 7, "CO_20110-1110", "CO2", "CO2_1409", (("01(1)0-11(1)0", "1409"),),
          (("0,081", "5"),), (R("435,8", "1,86", "5, 10"),)),
    Block(PART2, 8, "CO_20110-0310", "CO2", "CO2_1265", (("01(1)0-03(1)0", "1265"),),
          (("0,127", "5"),), (R("435,8", "1,18", "5, 10"),)),
    Block(PART3, 3, "N_2Ov_1", "N2O", "N2O_nu1", (("nu1", "1285"),), (),
          (R("488,0", "2,50±0,25", "12"), R("514,5", "3,07", "22"))),
    Block(PART3, 4, "N_2Ov_3", "N2O", "N2O_nu3", (("nu3", "2224"),), (),
          (R("488,0", "0,62±0,06", "12"), R("514,5", "0,60±0,06", "22"))),
    Block(PART3, 16, "CH_4v_1", "CH4", "CH4", (("nu1", "2917"),),
          (("0", "5, 7, 10"),),
          (R("435,8", "7,4", "23"), R("435,8", "7,0", "5, 10"),
           R("514,5", "7,6±0,4", "7"), R("337,1", "4,3±1,3", "37"),
           R("337,1", "7,2", "43"), R("514,5", "8,0±0,8", "11"),
           R("488,0", "6,0±0,6", "12"))),
)  # fmt: skip


def _check(block: Block, cells: list[str]) -> None:
    """Check the transcription of one block against its library cells."""
    if cells[0] != block.label:
        msg = f"{block.part} row {block.row} is {cells[0]!r}, not {block.label!r}"
        raise ValueError(msg)
    kept = [ratio for ratio in block.ratios if not ratio.lost]
    ratio_cell = cells[6]
    for word in WORDS:
        ratio_cell = ratio_cell.replace(word, "")
    expected = {
        1: "".join(shift for _, shift in block.shifts),
        2: "".join(printed for printed, _ in block.rho),
        5: "".join(ratio.wavelength for ratio in kept if not ratio.wavelength_lost),
        # N2's ratio is the reference, printed 1 (not transcribed).
        6: "".join(ratio.printed for ratio in kept)
        or ("1" if block.tag == "N2" else ""),
    }
    actual = {1: cells[1], 2: cells[2], 5: cells[5], 6: ratio_cell}
    for column, text in expected.items():
        if actual[column] != text:
            msg = f"{block.label} column {column}: {actual[column]!r} is not {text!r}"
            raise ValueError(msg)


def _split(printed: str) -> tuple[str, str]:
    """Return the value and uncertainty of a printed entry, with decimal points."""
    text = printed.replace(",", ".")
    for sign in ("±", "+"):
        if sign in text:
            value, sigma = text.split(sign)
            return value, sigma
    return text, ""


def _table_1() -> dict[str, list[Point]]:
    """Return Table 1's in-scope lines as data files."""
    grids = {part: json_table(DIRECTORY, part) for part in (PART1, PART2, PART3)}
    files: dict[str, list[Point]] = {}
    for block in BLOCKS:
        _check(block, grids[block.part][block.row])
        stem = f"{block.material}/danichkin_1981__"
        files[f"{stem}shift_{block.tag}.csv"] = [
            Point(row, row, shift, "", f"{block.label} | {shift}")
            for row, shift in block.shifts
        ]
        if block.rho:
            points = []
            for index, (printed, ref) in enumerate(block.rho):
                row = f"ref. [{ref}]"
                if [r for _, r in block.rho].count(ref) > 1:  # ozone's braced four
                    row += f", line {index + 1}"
                value, sigma = _split(printed)
                points.append(Point(row, SYSTEM, value, sigma, f"{printed} | [{ref}]"))
            files[f"{stem}rho_{block.tag}.csv"] = points
        if block.ratios:
            points = []
            for ratio in block.ratios:
                wavelength = ratio.wavelength.replace(",", ".")
                row = f"ref. [{ratio.ref}], {wavelength} nm"
                if ratio.note:
                    row += f", {ratio.note.split(' (')[0]}"
                cell = f"{ratio.printed} {ratio.note}".strip()
                printed = f"{ratio.wavelength} | {cell} | [{ratio.ref}]"
                value, sigma = _split(ratio.printed)
                extra = {"excitation_nm": wavelength}
                points.append(Point(row, SYSTEM, value, sigma, printed, extra))
            files[f"{stem}ratio_{block.tag}.csv"] = points
    return files


def _table_2() -> list[Point]:
    """Return Table 2's absolute N2 cross sections (converted to 488 nm)."""
    points = []
    for index, cells in enumerate(table(DIRECTORY, TABLE_2)[1:]):
        wavelength, value, ref = cells[0], cells[2], cells[3]
        printed = " | ".join(cells[:4])
        extra = {"excitation_nm": wavelength}
        if value == "5.5 and 6.2 (entire band)":
            # One printed row for Murphy et al. [10] and Yoshino and Bernstein [5]:
            # I_N2/I_CH4 = 0.14 (Q branch) and 0.17 (entire band).
            for part, number in (("Q branch", "5.5"), ("entire band", "6.2")):
                row = f"refs. [10][5], {wavelength} nm, {part}"
                points.append(Point(row, SYSTEM, number, "", printed, extra))
            continue
        # Ref. [12] has two rows at 488.0 nm (two methods): the row number tells.
        row = f"row {index + 1}, ref. {ref}, {wavelength} nm"
        points.append(Point(row, SYSTEM, value, "", printed, extra))
    return points


def _table_3() -> dict[str, list[Point]]:
    """Return Table 3's recommended depolarization ratios and cross-section ratios."""
    gases = {"N_{2}": "N2", "O_{2}": "O2", "CO_{2}": "CO2", "CO": "CO", "H_{2}": "H2",
             "O_{3}": "O3", "CH_{4}": "CH4"}  # fmt: skip
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_3)[1:]:
        if cells[0] not in gases:
            continue
        gas = gases[cells[0]]
        transition = (
            cells[1].replace("—", "-").replace("v_{1}", "nu1").replace("2v_{2}", "2nu2")
        )
        printed = " | ".join(cells)
        stem = f"{gas}/danichkin_1981__recommended_"
        files.setdefault(f"{stem}rho_{gas}.csv", []).append(
            Point(transition, SYSTEM, *_split(cells[3]), printed)
        )
        if gas != "N2":  # N2's ratio is 1 by definition
            files.setdefault(f"{stem}ratio_{gas}.csv", []).append(
                Point(transition, SYSTEM, *_split(cells[4]), printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the in-scope values of Tables 1-3 as data files."""
    files = _table_1()
    files["N2/danichkin_1981__absolute_N2.csv"] = _table_2()
    files.update(_table_3())
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

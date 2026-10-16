# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Schrötter and Klöckner (1979): Raman scattering cross sections of gases.

Chapter 4 (pp. 123-166) of *Raman Spectroscopy of Gases and Liquids*, A. Weber, ed.,
Topics in Current Physics 11 (Springer, 1979); the library's copy is the whole book
(``Weber_1979_RamanSpectroscopyGases``, formerly ``Unverified_6d54463e97c7``; PDF
page n is book page n - 10 in this chapter). Every table was checked against its
image, and the continuation order of Table 4.2 against the printed page numbers (137,
138, 139, 140, 141).

* Table 4.1 (``tbl_d34b5e440027``, p. 130): the absolute differential cross section of
  the N2 Q branch at each exciting line, with the chapter's normalised value
  (dsigma/dOmega divided by (nu_0 - 2331 cm^-1)^4) and the reference. One file per
  reference, the absolute value as the value and the printed row (both columns) kept.
* Table 4.2 (``tbl_5220e2a41fe5``, ``tbl_3bb750412546``, ``tbl_ce7fb84505c6``,
  ``tbl_20201534676f``, ``tbl_b9faf9e8a520``, pp. 137-141): the relative normalised
  cross sections Sigma_j (Eq. 4.14) of the Q branches, one column per exciting line (the
  column heads, 287 ... 633 nm, are the abscissa as printed), each value with a
  footnote letter naming its reference. One file per molecule, band and reference.
  Only the curated species are taken (N2 is the reference, 1 by definition, and is not
  written); part 3 (p. 139) has none.
* Table 4.4 (``tbl_a794533c145e``, p. 149): CH4, per mode, the shift, the depolarization
  ratio, Sigma_j and the absolute normalised cross section, with references.

STRUCTURAL FIX (Tables 4.1, 4.2 and 4.4). The library's cells join the stacked lines
of a cell without separators ('44 ± 1743 ± 242 ± 243.2 ± 0.8', '1.2_5^e1.0_9^f',
'0.520.66'), and in Table 4.2 the N2O row puts the lines of the 2224 cm^-1 band into
the cells of the 1285 cm^-1 band (the image shows 0.5_0^d at 488 nm and 0.5_1^b at
515 nm on the 2224 cm^-1 line, and nothing else on it). The lines are therefore
transcribed below from the images, each with its reference, and the script checks
that the lines of each cell, concatenated, are exactly the library cell (after
removing spaces and the markup ``_ ^ { }``), and that no in-scope cell is left
unread.

Table 4.2 prints the last digit of many values lowered ('1.2₅'): the transcription
keeps it as ``1.2_5`` in the printed column and the value is 1.25. Every value carries
a raised footnote letter, the reference (the key under the table, p. 141); the value
and its letter are written to the printed column as '1.2_5^e'.
"""

import re
from dataclasses import dataclass

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Weber_1979_RamanSpectroscopyGases"
TABLE_4_1 = "tbl_d34b5e440027"
# Table 4.2 in page order (pp. 137-141, the printed page numbers of the images).
TABLE_4_2 = (
    "tbl_5220e2a41fe5",
    "tbl_3bb750412546",
    "tbl_ce7fb84505c6",
    "tbl_20201534676f",
    "tbl_b9faf9e8a520",
)
TABLE_4_4 = "tbl_a794533c145e"
STEM = "schrotter_1979__"
WAVELENGTHS = ("287", "337", "347", "351", "364", "436", "458", "488", "515", "633")
_MARKUP = re.compile(r"[\s_^{}$]")


def _flat(text: str) -> str:
    """Return a cell without spaces and markup, for comparing with a transcription."""
    return _MARKUP.sub("", text)


def _plus_minus(cell: str) -> tuple[str, str]:
    """Return the value and uncertainty of a '43.2 ± 0.8' line ('970': no sigma)."""
    value, _, sigma = cell.partition("±")
    return value.strip(), sigma.strip()


# ---------------------------------------------------------------------------------
# Table 4.1: per library row, the stacked lines (dsigma/dOmega, normalised, reference)
# as the image shows them.
TABLE_4_1_LINES = {
    1: ("632.8", (("21 ± 3", "6.4 ± 1", "4.52"),)),
    2: (
        "514.5",
        (
            ("44 ± 17", "5.1 ± 2", "4.37"),
            ("43 ± 2", "5.0 ± 0.3", "4.38"),
            ("42 ± 2", "4.9 ± 0.3", "4.22"),
            ("43.2 ± 0.8", "5.05 ± 0.11", "4.45"),
        ),
    ),
    3: (
        "488.0",
        (
            ("33 ± 11", "3.0 ± 1", "4.21"),
            ("43", "4.0 ± 1", "4.21"),
            ("54 ± 3", "5.0 ± 0.3", "4.22"),
            ("55.8 ± 2", "5.13 ± 0.2", "4.45"),
        ),
    ),
    4: (
        "457.9",
        (("76 ± 5", "5.2 ± 0.4", "4.22"), ("73.7 ± 3", "5.09 ± 0.25", "4.45")),
    ),
    5: ("435.8", (("92 ± 10", "5.1 ± 0.5", "4.10"),)),
    6: ("363.8", (("204 ± 25", "5.1 ± 0.6", "4.45"),)),
    7: ("351.1", (("243 ± 30", "5.2 ± 0.7", "4.45"),)),
    8: ("300", (("970", "(10.5)", "4.44"),)),
}
# Reference -> file tag. Fenner et al. [4.21] give two values at 488.0 nm: 'too low a
# value by the latter method [an absolute spectrometer calibration] and a somewhat
# better, but still low result by comparison with liquid benzene' (p. 127), in that
# order in the table; their rows are named so.
TABLE_4_1_FILES = {
    "4.52": "kamiyama",
    "4.37": "fouche",
    "4.38": "penney_1972",
    "4.22": "hyatt",
    "4.45": "kloeckner",
    "4.21": "fenner",
    "4.10": "murphy",
    "4.44": "penney_1974",
}
FENNER_ROWS = ("488.0 absolute calibration", "488.0 against benzene")


def _table_4_1() -> dict[str, list[Point]]:
    grid = json_table(DIRECTORY, TABLE_4_1)
    if not grid[0][0].startswith("Exciting line") or grid[0][3] != "Ref.":
        msg = f"unexpected Table 4.1 header {grid[0]!r}"
        raise ValueError(msg)
    if not grid[9][0].startswith("Weighted average") or grid[9][2] != "5.05 ± 0.08":
        msg = f"unexpected Table 4.1 last row {grid[9]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    fenner = iter(FENNER_ROWS)
    for row, (wavelength, lines) in TABLE_4_1_LINES.items():
        cells = grid[row]
        if cells[0] != wavelength:
            msg = f"Table 4.1 row {row} is {cells[0]!r}, not {wavelength}"
            raise ValueError(msg)
        for column, joined in ((1, 0), (2, 1)):
            if _flat("".join(line[joined] for line in lines)) != _flat(cells[column]):
                msg = f"Table 4.1 row {row} column {column}: {cells[column]!r}"
                raise ValueError(msg)
        if "".join(f"[{line[2]}]" for line in lines) != cells[3]:
            msg = f"Table 4.1 row {row} references: {cells[3]!r}"
            raise ValueError(msg)
        for absolute, normalised, ref in lines:
            value, sigma = _plus_minus(absolute)
            label = next(fenner) if ref == "4.21" else wavelength
            printed = f"{wavelength} | {absolute} | {normalised} | [{ref}]"
            name = f"N2/{STEM}table4_1_{TABLE_4_1_FILES[ref]}_N2.csv"
            files.setdefault(name, []).append(
                Point(label, wavelength, value, sigma, printed)
            )
    return files


# ---------------------------------------------------------------------------------
# Table 4.2.


@dataclass(frozen=True)
class Band:
    """One in-scope molecule and band of Table 4.2, transcribed from the image."""

    part: int
    """Index into TABLE_4_2."""
    row: int
    """The library row of the band (its label and wavenumber cells)."""
    material: str
    label: str
    """The molecule cell as the library has it (checked)."""
    shift: str
    """The band's nu_j as printed (cm^-1)."""
    tag: str
    """File tag: material and band."""
    lines: dict[str, tuple[str, ...]]
    """Per exciting wavelength, the lines 'value^letter' in the image's order."""
    borrowed: dict[str, tuple[str, ...]]
    """Lines of this band that the library puts into another row's cell (N2O)."""
    host_row: int = 0
    """The library row whose cells hold the borrowed lines."""


def _lines(**columns: str) -> dict[str, tuple[str, ...]]:
    """Return {wavelength: lines} from keyword arguments nm337="1.2_5^e 1.0_9^f"."""
    return {
        key.removeprefix("nm"): tuple(text.split()) for key, text in columns.items()
    }


BANDS = (
    Band(0, 3, "O2", "O_2", "1555", "O2",
         _lines(nm337="1.2_5^e 1.0_9^f", nm347="1.0_8^h", nm351="1.0^p",
                nm364="1.0^p", nm436="1.0_5^a", nm458="1.0^p",
                nm488="1.0_9^d 1.0^p 1.0_3^dd", nm515="1.0_2^b 1.0_2^c 1.0^p",
                nm633="0.9_1^o"),
         {}),
    Band(0, 4, "H2", "H_2", "4156", "H2",
         _lines(nm337="4.2^e", nm436="4.1_5^a", nm488="3.6_4^d", nm515="3.4_5^b"),
         {}),
    Band(0, 5, "D2", "D_2", "2986", "D2", _lines(nm436="2.6^a"), {}),
    Band(0, 14, "CO", "CO", "2143", "CO",
         _lines(nm337="0.8_7^e 0.9_6^f", nm347="0.7_9^h", nm351="0.9_8^p",
                nm364="0.9_6^p", nm436="0.9_0^a", nm458="0.9_2^p",
                nm488="0.9_6^d 0.9_3^p 0.7_2^dd", nm515="0.8_7^b 0.9_6^c 0.9_5^p",
                nm633="1.0_4^o"),
         {}),
    Band(1, 3, "CO2", "CO_2", "1285", "CO2_1285",
         _lines(nm337="0.7_2^e 0.8_0^f", nm436="0.7_3^a", nm488="0.7_1^d 0.9_2^dd",
                nm515="0.7_8^b", nm633="0.7_8^o"),
         {}),
    Band(1, 4, "CO2", "", "1388", "CO2_1388",
         _lines(nm337="1.0_4^e 1.1_9^f", nm347="1.0_3^h", nm436="1.0_7^a",
                nm488="1.1_4^d 1.3_5^dd", nm515="1.2_1^b 1.2_1^c", nm633="1.1_3^o"),
         {}),
    Band(1, 5, "N2O", "N_2O", "1285", "N2O_1285",
         _lines(nm347="2.5^h", nm351="2.0^p", nm364="2.0^p", nm458="1.8_5^p",
                nm488="1.7_6^d 1.8^p", nm515="2.1_2^b 1.8^p"),
         {}),
    # The 2224 cm^-1 band's two lines are in the library's 1285 cm^-1 cells (row 5),
    # after that band's lines; the image aligns them with the 2224 label.
    Band(1, 6, "N2O", "", "2224", "N2O_2224",
         {}, _lines(nm488="0.5_0^d", nm515="0.5_1^b"), host_row=5),
    Band(1, 11, "O3", "O_3", "1103", "O3", _lines(nm287="816^n", nm515="3^b"), {}),
    Band(3, 8, "CH4", "CH_4", "2917", "CH4",
         _lines(nm337="8.2^e", nm347="7.8^h", nm351="8.4^p", nm364="8.3^p",
                nm436="8.0^a 9.5^j", nm458="9.2^p", nm488="6.8^d 9.2^p 8.7^r",
                nm515="9.1^b 8.7^c 9.3^p"),
         {}),
    Band(4, 1, "SF6", "SF_6", "775", "SF6",
         _lines(nm337="3.3^e", nm436="3.4^a", nm488="3.7_6^aa", nm515="3.9_6^aa"),
         {}),
)  # fmt: skip
# The rows of Table 4.2 holding only out-of-scope molecules, per part (checked, so
# that an in-scope row cannot be missed).
OUT_OF_SCOPE = {
    0: {"N_2", "F_2", "Cl_2", "Br_2", "I_2", "HF", "HCl", "HBr", "HI"},
    1: {"NO", "ONF", "H_2O", "H_2S", "SO_2"},
    2: {"SO_2", "NO_2", "C_2H_2", "C_2N_2", "NH_3", "ND_3", "C_2H_4",
        "cis-C_2H_2Cl_2", "C_2H_3Cl", "CH_3OH", "C_2H_6"},
    3: {"C_2H_5I", "C_2H_5OH", "C_3H_6", "CD_4", "CCl_4", "CH_2Cl_2", "CHCl_3"},
    4: {"UF_6", "C_3H_8", "C_4H_6", "(CH_3)_2CO", "C_6H_6", "C_6D_6", "C_6H_{12}"},
}  # fmt: skip
_LINE = re.compile(r"^(?P<value>[\d.]+(?:_\d)?)\^(?P<letter>[a-z]+)$")
# The footnote key (the last library row of part 5, p. 141) for the letters taken.
KEY = (
    "^a[4.10]", "^b[4.36,37]", "^c[4.38]", "^d[4.21]", "^e[4.39]", "^f[4.40]",
    "^h[4.56]", "^j[4.11]", "^n[4.44] ", "^o[4.52]", "^p[4.54]", "^r[4.61]",
    "^aa[4.100]", "^dd[4.102]",
)  # fmt: skip


def _header_row(part: int) -> int:
    """Return the library row of the wavelength heads (part 5 lost its top row)."""
    return 0 if part == 4 else 1


def _check_part(part: int, grid: list[list[str]]) -> None:
    """Check the wavelength heads and that every molecule row is known."""
    head = _header_row(part)
    if tuple(grid[head][2:]) != WAVELENGTHS:
        msg = f"Table 4.2 part {part + 1} heads {grid[head]!r}"
        raise ValueError(msg)
    known = OUT_OF_SCOPE[part] | {band.label for band in BANDS if band.part == part}
    for row in grid[head + 1 :]:
        if row[0].startswith("^a[4.10]"):
            missing = [entry for entry in KEY if entry not in row[0]]
            if missing:
                msg = f"Table 4.2 key lacks {missing!r}"
                raise ValueError(msg)
        elif row[0] and row[0] not in known:
            msg = f"Table 4.2 part {part + 1}: unexpected molecule {row[0]!r}"
            raise ValueError(msg)


def _check_band(band: Band, grids: list[list[list[str]]]) -> None:
    """Check a band's transcription against its library cells."""
    grid = grids[band.part]
    cells = grid[band.row]
    if (cells[0], cells[1]) != (band.label, band.shift):
        msg = (
            f"Table 4.2 row {band.row} is {cells[:2]!r}, not {band.label}/{band.shift}"
        )
        raise ValueError(msg)
    for column, wavelength in enumerate(WAVELENGTHS, start=2):
        own = "".join(band.lines.get(wavelength, ()))
        if band.host_row:
            # The library cell of this row is empty: the lines are in the host's.
            if cells[column]:
                msg = f"Table 4.2 row {band.row} {wavelength}: {cells[column]!r}"
                raise ValueError(msg)
            continue
        guest = "".join(
            "".join(other.borrowed.get(wavelength, ()))
            for other in BANDS
            if other.part == band.part and other.host_row == band.row
        )
        if _flat(own + guest) != _flat(cells[column]):
            msg = (
                f"Table 4.2 {band.tag} {wavelength}: {cells[column]!r} != {own + guest}"
            )
            raise ValueError(msg)


def _table_4_2() -> dict[str, list[Point]]:
    grids = [json_table(DIRECTORY, table_id) for table_id in TABLE_4_2]
    for part, grid in enumerate(grids):
        _check_part(part, grid)
    # N2 is the reference: 1 at every exciting line (not written).
    if grids[0][2][:2] != ["N_2", "2331"] or set(grids[0][2][2:]) != {"1"}:
        msg = f"unexpected N2 row {grids[0][2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for band in BANDS:
        _check_band(band, grids)
        for wavelength, lines in (band.lines | band.borrowed).items():
            for line in lines:
                match = _LINE.match(line)
                if match is None:
                    msg = f"not a 'value^letter' line: {line!r}"
                    raise ValueError(msg)
                value = match.group("value").replace("_", "")
                letter = match.group("letter")
                name = f"{band.material}/{STEM}sigma_{band.tag}_{letter}.csv"
                files.setdefault(name, []).append(
                    Point(wavelength, wavelength, value, "", f"{band.shift} | {line}")
                )
    return files


# ---------------------------------------------------------------------------------
# Table 4.4, CH4 (library rows 1-4): per mode the lines (Sigma_j, normalised absolute,
# reference) as the image shows them; nu_1's reference cell is empty (its Sigma_j
# carries the note a, 'average from Table 4.2').
TABLE_4_4_MODES = {
    1: ("nu1(a1)", "CH_4v_1(a_1)", (("8.56a", "43.2", ""),)),
    2: ("nu2(e)", "v_2(e)", (("0.52", "2.6", "4.61"), ("0.66", "3.3", "4.11"))),
    3: ("nu3(f2)", "v_3(f_2)", (("4.82", "24.3", "4.61"), ("4.0", "20.3", "4.11"))),
    4: ("nu4(f2)", "v_4(f_2)", (("<0.01", "<0.04", "4.61"),)),
}
TABLE_4_4_FILES = {"": "average", "4.61": "bermejo", "4.11": "kagel"}


def _table_4_4() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_4_4)
    if grid[0][:4] != [
        "Molecule and vibration",
        "\\nu_{j}[cm^{-1}]",
        "^{\\rho}_{j}",
        "\\Sigma_{j}",
    ]:
        msg = f"unexpected Table 4.4 header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for row, (mode, label, lines) in TABLE_4_4_MODES.items():
        cells = grid[row]
        if cells[0] != label:
            msg = f"Table 4.4 row {row} is {cells[0]!r}, not {label}"
            raise ValueError(msg)
        joined = (
            "".join(line[0] for line in lines),
            "".join(line[1] for line in lines),
            "".join(f"[{line[2]}]" for line in lines if line[2]),
        )
        if joined != (cells[3], cells[4], cells[5]):
            msg = f"Table 4.4 row {row}: {cells!r} != {joined!r}"
            raise ValueError(msg)
        printed = " | ".join(cells)
        files.setdefault(f"CH4/{STEM}table4_4_shift_CH4.csv", []).append(
            Point(mode, mode, cells[1], "", printed)
        )
        files.setdefault(f"CH4/{STEM}table4_4_rho_CH4.csv", []).append(
            Point(mode, "not stated", cells[2], "", printed)
        )
        for sigma, absolute, ref in lines:
            if sigma.startswith("<"):
                continue  # an upper bound (nu4 < 0.01): in the record's notes
            line = f"{sigma} | {absolute} | [{ref}]" if ref else f"{sigma} | {absolute}"
            name = f"CH4/{STEM}table4_4_sigma_{TABLE_4_4_FILES[ref]}_CH4.csv"
            files.setdefault(name, []).append(
                Point(mode, mode, sigma.removesuffix("a"), "", f"{printed} || {line}")
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the files of Tables 4.1, 4.2 and 4.4."""
    return _table_4_1() | _table_4_2() | _table_4_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

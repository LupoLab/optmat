# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chan, Cooper, Guo & Brion (1992), part II: neon oscillator strengths.

Continuum (low-resolution dipole (e,e), 1 eV FWHM, TRK-normalised): Table I
(``tbl_66ad40730b41``) and its continuation (``tbl_b252935d2cb5``) list photon energy
(eV) and df/dE in ``10^-2 eV^-1`` in three column pairs per page, read column by
column (images checked: every cell agrees). One file; the rows are the energies (all
distinct), in reading order.

Discrete lines (high-resolution dipole (e,e), 0.048 eV FWHM, and literature values):
Tables II (``tbl_3545e6ef35ab``, the 3s and 3s' resonance lines f_1 and f_2), III
(``tbl_bbf7e53a7149``, 4s ... 4d') and IV (``tbl_f0fb0b7927fd``, 6s ... 5d' and the
"Total to ionization"), cells files only (images checked). Unlike part III (argon,
krypton, xenon: ``chan_1992``), the column headers name the transitions without their
energies. The text gives the energies of the two resonance lines only, 16.671 eV (f_1)
and 16.848 eV (f_2) (p. 2), so each source row becomes up to three files:

- ``<source>_Ne``: the resonance lines, rows and abscissa the text's energies (eV);
- ``<source>_higher_Ne``: the lines of Tables III and IV, whose energies the paper does
  not print: rows are the transition labels ("4s", "4s'", ...), the abscissa an index,
  the column order of the tables (1 = 4s ... 12 = 5d');
- ``<source>_sum_Ne``: summed values of a pair nl + nl' (note a, "Total oscillator
  strength (f_1 + f_2)" in Table II, "Summed oscillator strength as indicated" in
  Tables III and IV). The images print each sum centred between the two columns of its
  pair; the cells put it under one of them (Table II: under f_2; Tables III and IV:
  under the unprimed line). Rows are the pairs ("3s+3s'", ...), indexed 1 to 7.

Handled structurally:

- per-cell uncertainties are printed either inline, ``0.0129(0.0006)`` (Tables III,
  IV), or on the next row, in parentheses under the value (Table II, where that row's
  label is empty or the method in parentheses); both go to the sigma column. Spaces
  inside digit groups (``0.003 30``) are removed;
- rows holding only a method in parentheses ("(Total absorption)") are skipped; the
  method is in the source record;
- "(a) dipole length", "(b) dipole velocity", "(a) wave function", "(b)
  semiempirical", Klose's "(i) IC and HFS", "(ii) IC and CF" and Lawrence and Liszt's
  "(1)" and "(2)^c" are variants of the source named on the row above;
- Aleksandrov et al. [32] appear under Theory and Experiment: the section is part of
  the file name of the calculation;
- the "Total to ionization" column (Table IV) goes to ``<source>_total_Ne`` with an
  index abscissa (the discrete sum up to 21.616 eV, not a line);
- footnote marks (``^a``, ``^b``, ``^c``) are stripped; the record says what they mean.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Chan_1992_AbsoluteOpticalOscillator_ec241a"
CONTINUUM = ("tbl_66ad40730b41", "tbl_b252935d2cb5")
TABLE_II = "tbl_3545e6ef35ab"
TABLE_III = "tbl_bbf7e53a7149"
TABLE_IV = "tbl_f0fb0b7927fd"
PREFIX = "Ne/chan_1992_neon__"

# Resonance-line energies (eV), from the text (p. 2): f_1 and f_2.
RESONANCE = {"(2P_{3/2})3s (f_1)": "16.671", "(2P_{1/2})3s' (f_2)": "16.848"}
# Tables III and IV: header label -> (row id, index in column order).
HIGHER = {
    "(^{2}P_{3/2})4s": "4s",
    "(^{2}P_{1/2})4s'": "4s'",
    "(^{2}P_{3/2})3d": "3d",
    "(^{2}P_{1/2})3d'": "3d'",
    "(^{2}P_{3/2})5s": "5s",
    "(^{2}P_{1/2})5s'": "5s'",
    "(^{2}P_{3/2})4d": "4d",
    "(^{2}P_{1/2})4d'": "4d'",
    "(^{2}P_{3/2})6s": "6s",
    "(^{2}P_{1/2})6s'": "6s'",
    "(^{2}P_{3/2})5d": "5d",
    "(^{2}P_{1/2})5d'": "5d'",
}
ORDER = list(HIGHER.values())
PAIRS = ["3s+3s'", "4s+4s'", "3d+3d'", "5s+5s'", "4d+4d'", "6s+6s'", "5d+5d'"]

# Row label (without footnote marks) -> file stem. Aleksandrov et al. appear in both
# sections; the theory row gets "_theory".
SOURCES = {
    "Amus'ya (1990) [50]": "amusya_1990",
    "Amusy'ya (1990) [50]": "amusya_1990",  # sic, Table III
    "Kelly (1964) [49]": "kelly_1964",
    "Cooper (1962) [4]": "cooper_1962",
    "Aleksandrov et al. (1983) [32]": "aleksandrov_1983",
    "Aleksandrov et al. (1983) [32](Total absorption)": "aleksandrov_1983",
    "Stewart (1975) [55]": "stewart_1975",
    "Albat and Gruen (1974) [54]": "albat_1974",
    "Gruzdev and Loginov (1973) [53]": "gruzdev_loginov_1973",
    "Aymar, Feneuille, and Klapisch (1970) [52]": "aymar_1970",
    "Gruzdev (1967) [3]": "gruzdev_1967",
    "Gold and Knox (1959) [51]": "gold_1959",
    "Klose (1969) [34]": "klose_1969",
    "Klose (1969) [75]": "klose_1969b",
    "Present work [(HR dipole (e,e)]": "present",
    "Present work [HR dipole (e,e)]": "present",
    "Tsurubuchi, Watanabe, and Arikawa (1990) [44]": "tsurubuchi_1990",
    "Westerveld, Mulder, and Van Eck (1979) [43]": "westerveld_1979",
    "Bhaskar and Luiro (1976) [37]": "bhaskar_1976",
    "Knystautas and Drouin (1974) [35]": "knystautas_1974",
    "Irwin, Livingston, and Kernahan (1973) [36]": "irwin_1973",
    "Natali, Kuyatt, and Mielczarek (1973) [47]": "natali_unpublished",
    "Natali, Kuyatt, and Mielczarek (1973) [47](Electron impact)": (
        "natali_unpublished"
    ),
    "Jongh and Eck (1971) [42]": "de_jongh_1971",
    "Kazantsev and Chaika (1971) [38]": "kazantsev_1971",
    "Geiger (1970) [45,46]": "geiger_1970",
    "Lawrence and Liszt (1969) [33]": "lawrence_1969",
    "Lawrence and Liszt (1960) [33]": "lawrence_1969",  # sic, Table III
    "Lewis (1967) [41]": "lewis_1967",
    "Korolev, Odinstov, and Fursova (1964) [31]": "korolev_1964",
    "Ducloy (1973\N{EN DASH}1974) [40]": "ducloy_1973",
    "Klose (1970) [34]": "klose_1969_lifetime",  # Ref. [34] is Klose (1969)
    "Decomps and Dumont (1968) [39]": "decomps_1968",
}
VARIANTS = {
    "(a) Dipole length": "length",
    "(b) Dipole velocity": "velocity",
    "(a) Wave function": "wave_function",
    "(b) Semiempirical": "semiempirical",
    "(i) IC and HFS": "ic_hfs",
    "(ii) IC and CF": "ic_cf",
    "(Lifetime: delayed coincidence) (1)": "",  # Lawrence and Liszt's own values
    "(2)": "westerveld",  # note c: recalculated by Westerveld et al.
}
_MARK = re.compile(r"\^[a-c]")
_VALUE = re.compile(r"^(?P<value>\d*\.\d+)(?:\((?P<sigma>\d*\.\d+)\))?$")
_SIGMA = re.compile(r"^\((?P<sigma>\d*\.\d+)\)$")


def _continuum() -> list[Point]:
    points = []
    for table_id in CONTINUUM:
        rows = table(DIRECTORY, table_id)[1:]
        for column in range(0, 6, 2):
            for cells in rows:
                energy, strength = cells[column].strip(), cells[column + 1].strip()
                if not energy and not strength:
                    continue  # the short last columns of the continuation
                printed = f"{energy} | {strength}"
                points.append(
                    Point(energy, number(energy), number(strength), "", printed)
                )
    return points


def _cell(cell: str) -> tuple[str, str]:
    """Return the value and the inline uncertainty of a value cell."""
    match = _VALUE.match(_MARK.sub("", cell).replace(" ", ""))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    return match["value"], match["sigma"] or ""


def _columns(header: list[str]) -> dict[int, tuple[str, str, str]]:
    """Map a header row's columns to (kind, row id, abscissa)."""
    columns: dict[int, tuple[str, str, str]] = {}
    for column, label in enumerate(header):
        if label in RESONANCE:
            columns[column] = ("line", RESONANCE[label], RESONANCE[label])
        elif label in HIGHER:
            row = HIGHER[label]
            columns[column] = ("higher", row, str(ORDER.index(row) + 1))
    return columns


def _pair(row: str) -> Point:
    """Return the sum-dataset point of the pair whose unprimed member is ``row``."""
    name = f"{row}+{row}'"
    return Point(name, str(PAIRS.index(name) + 1), "", "", "")


def _attach_sigmas(
    cells: list[str], last: list[tuple[str, int]], files: dict[str, list[Point]]
) -> None:
    """Put a row of parenthesised uncertainties on the previous row's points."""
    sigmas = [v.strip() for v in cells[1:] if v.strip()]
    for (path, index), sigma in zip(last, sigmas, strict=True):
        point = files[path][index]
        found = _SIGMA.match(sigma)
        if found is None:
            msg = f"unreadable uncertainty {sigma!r}"
            raise ValueError(msg)
        printed = f"{point.printed} || {' | '.join(cells)}"
        files[path][index] = Point(
            point.row, point.x, point.value, found["sigma"], printed
        )


def _stem(plain: str, parent: str, section: str) -> str:
    """Return the file stem of a source row, given its parent and section."""
    if plain in VARIANTS:
        suffix = VARIANTS[plain]
        return f"{parent}_{suffix}" if suffix else parent
    stem = SOURCES[plain]
    if stem == "aleksandrov_1983" and section == "(a) Theory":
        return "aleksandrov_1983_theory"
    return stem


def _point(
    cells: list[str], column: int, stem: str, headers: list[list[str]], *, summed: bool
) -> tuple[str, Point]:
    """Return the file and point of one value cell."""
    value, sigma = _cell(cells[column])
    printed_row = " | ".join(cells)
    if headers[0][column] == "Total to ionization":
        return f"{PREFIX}{stem}_total_Ne.csv", Point(
            "total", "1", value, sigma, printed_row
        )
    kind, row, abscissa = _columns(headers[1])[column]
    if summed:
        pair = _pair("3s" if kind == "line" else row)
        return f"{PREFIX}{stem}_sum_Ne.csv", Point(
            pair.row, pair.x, value, sigma, printed_row
        )
    suffix = "" if kind == "line" else "_higher"
    printed = f"{headers[1][column]} | {printed_row}"
    return f"{PREFIX}{stem}{suffix}_Ne.csv", Point(row, abscissa, value, sigma, printed)


def _discrete_table(table_id: str, files: dict[str, list[Point]]) -> None:
    grid = json_table(DIRECTORY, table_id)
    section, parent = "", ""
    last: list[tuple[str, int]] = []  # (file, index) of the previous row's points
    for cells in grid[2:]:
        label = cells[0].strip()
        values = [c.strip() for c in cells[1:]]
        if label in {"(a) Theory", "(b) Experiment"}:
            section = label
        elif any(values) and all(_SIGMA.match(v) or not v for v in values):
            _attach_sigmas(cells, last, files)
            last = []
        elif not any(values):
            if not label.startswith("("):  # a parent of variant rows
                parent = SOURCES[_MARK.sub("", label)]
            # else a method line such as "(Total absorption)": skipped
        else:
            plain = _MARK.sub("", label)
            stem = _stem(plain, parent, section)
            if plain not in VARIANTS:
                parent = stem
            summed = label.endswith("^a") or (
                table_id == TABLE_II and "^a" in "".join(values)
            )
            last = []
            for column in range(1, len(cells)):
                if cells[column].strip():
                    path, point = _point(cells, column, stem, grid, summed=summed)
                    files.setdefault(path, []).append(point)
                    last.append((path, len(files[path]) - 1))


def extract() -> dict[str, list[Point]]:
    """Return the continuum df/dE and the discrete oscillator strengths of neon."""
    files: dict[str, list[Point]] = {f"{PREFIX}continuum_Ne.csv": _continuum()}
    for table_id in (TABLE_II, TABLE_III, TABLE_IV):
        _discrete_table(table_id, files)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

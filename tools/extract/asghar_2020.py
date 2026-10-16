# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Asghar, Piracha, Ali and Baig (2020): Ar I 4p-4s transition probabilities.

Appl. Opt. 59, 4240-4246 (2020) = PDF pages 1-7 (journal page = PDF page + 4239); a
vector PDF with a text layer. Table 1 has no CSV rendering and is read from its JSON
cells; Table 2 from its CSV.

- Table 1 (p. 4243, ``tbl_8164cd4b15b2``): the measured A (1e6 s-1, "Present Work") of
  the 27 Ar I 4p-4s lines with a value, with the printed absolute uncertainty (the
  738.40 nm cell prints ``8.23 0.68``, the ``±`` lost); the three lines printed ``-``
  (747.2, 857.8, 1148.81 nm) have none and are not written.
- Table 1, column "Lifetimes tau (ns)": the upper levels' lifetimes the A are
  normalised to, quoted from Chang and Setser (footnote b, eight levels) and Erman and
  Huldt (footnote a, 4p'[1/2]0 and 4p[1/2]0), one file each, with the printed
  absolute uncertainty (the percentage in parentheses is its ratio).

Each A row is the printed (air) wavelength, unique in the table; ``x`` is the
transition in the paper's Racah notation, upper-lower as printed
(``4p'[1/2]0-4s'[1/2]1``, the upper level spanning its group's rows). A lifetime
row and its ``x`` are the upper level.

**Checks** (:func:`checks`; the script fails if one is missed): every line's printed
wavelength identifies one NIST Ar I level pair (air wavelengths from NIST's levels,
``external/nist_asd/Ar_I_levels.tsv``), which is the printed label but for the two
misprinted lower levels of rows without a value (747.2 nm, 857.8 nm); the Paschen
column is the Racah label's; the printed percentage of each lifetime is its ratio; the
column "NIST [24]" is NIST ASD's A (``Ar_I_lines_nir.tsv``) to the printed digits; Eq.
(2) (sum of A times tau = 1 per level, the omitted lines being < 0.2e6 s-1) holds
within 2 % for six levels and fails for exactly 4p'[1/2]0, 4p'[1/2]1, 4p'[3/2]2 and
4p'[3/2]1; Table 2's line strengths and f are Table 1's A through Eqs. (3) and (4)
within 3 % but three f (``TABLE_2_OFF``).
"""

import csv
import re

from _json_table import json_table
from _library import DATA, Point, number, table, write
from wiese_1989 import _air, _levels

DIRECTORY = "Asghar_2020_MeasurementsExperimentalTransition"
TABLE_1 = "tbl_8164cd4b15b2"
TABLE_2 = "tbl_d60935dc2edf"
FILE_A = "Ar/asghar_2020__A_4p_4s.csv"
FILE_CHANG = "Ar/asghar_2020__lifetime_chang_1978.csv"
FILE_ERMAN = "Ar/asghar_2020__lifetime_erman_1978.csv"
NO_VALUE = {"-", "\N{EN DASH}", "\N{EM DASH}", "\N{MINUS SIGN}"}
_A = re.compile(r"^([\d.]+)\s*±?\s*([\d.]+)$")
_TAU = re.compile(r"^([ab])\s*([\d.]+)\s*±\s*([\d.]+)\s*\(([\d.]+)%\)$")
# Lower levels misprinted in rows without a value (not written): printed -> NIST's.
MISPRINTED = {"747.2": ("4s'[3/2]1", "4s[3/2]1"), "857.8": ("4s[1/2]1", "4s'[1/2]1")}
# Levels whose printed A do not satisfy Eq. (2) with the printed lifetime.
EQ2_OFF = {"4p'[1/2]0", "4p'[1/2]1", "4p'[3/2]2", "4p'[3/2]1"}
# Table 2 cells (S, f) not reproduced by Eqs. (3), (4) from Table 1's A within 3 %.
TABLE_2_OFF = {("714.71", "f"), ("747.2", "f"), ("751.47", "f")}
PASCHEN = {
    "4s[3/2]2": "1s5",
    "4s[3/2]1": "1s4",
    "4s'[1/2]0": "1s3",
    "4s'[1/2]1": "1s2",
    "4p'[1/2]0": "2p1",
    "4p'[1/2]1": "2p2",
    "4p'[3/2]2": "2p3",
    "4p'[3/2]1": "2p4",
    "4p[1/2]0": "2p5",
    "4p[3/2]2": "2p6",
    "4p[3/2]1": "2p7",
    "4p[5/2]2": "2p8",
    "4p[5/2]3": "2p9",
    "4p[1/2]1": "2p10",
}


def _rows() -> list[tuple[str, str, list[str]]]:
    """Return Table 1's rows as (upper, lower, cells), the upper level filled in."""
    rows = []
    upper = ""
    for cells in json_table(DIRECTORY, TABLE_1)[2:]:
        label = cells[0].replace(" ", "")
        if label.startswith("4p"):
            upper, lower = label.split("-")
        else:
            lower = label
        rows.append((upper, lower, cells))
    return rows


def _transitions() -> list[Point]:
    """Return the 27 lines with a measured A."""
    points = []
    for upper, lower, cells in _rows():
        if cells[4].strip() in NO_VALUE:
            continue
        match = _A.match(cells[4].strip())
        if not match:
            msg = f"{cells[2]}: not a value and uncertainty: {cells[4]!r}"
            raise ValueError(msg)
        value, sigma = (number(group) for group in match.groups())
        points.append(
            Point(cells[2], f"{upper}-{lower}", value, sigma, " | ".join(cells))
        )
    return points


def _lifetimes() -> dict[str, list[Point]]:
    """Return the quoted lifetimes by footnote (a Erman-Huldt, b Chang-Setser)."""
    files: dict[str, list[Point]] = {FILE_CHANG: [], FILE_ERMAN: []}
    for upper, _, cells in _rows():
        if not cells[3]:
            continue
        match = _TAU.match(cells[3].strip())
        if not match:
            msg = f"{upper}: not a footnoted lifetime: {cells[3]!r}"
            raise ValueError(msg)
        mark, value, sigma, _ = match.groups()
        name = FILE_ERMAN if mark == "a" else FILE_CHANG
        files[name].append(Point(upper, upper, value, sigma, cells[3]))
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table 1's A and the two sets of quoted lifetimes."""
    return {FILE_A: _transitions(), **_lifetimes()}


def _identify(levels: dict[str, float], wavelength: str) -> str:
    """Return the upper-lower 4p-4s pair whose air wavelength is the printed one."""
    printed = float(wavelength)
    tolerance = 0.5 * 10.0 ** -len(wavelength.split(".")[1]) + 0.001
    if wavelength == "747.2":  # NIST's 747.117 nm: printed 747.2, not 747.1
        tolerance = 0.1
    found = [
        f"{up}-{lo}"
        for up, e_up in levels.items()
        for lo, e_lo in levels.items()
        if up[:2] == "4p"
        and lo[:2] == "4s"
        and abs(_air(1e7 / (e_up - e_lo)) - printed) < tolerance
    ]
    if len(found) != 1:
        msg = f"{wavelength} nm: {found}"
        raise ValueError(msg)
    return found[0]


def _nist() -> dict[str, float]:
    """NIST ASD's A (1e6 s-1) of the Ar I 4p-4s lines, ``upper-lower``."""
    path = DATA / "external" / "nist_asd" / "Ar_I_lines_nir.tsv"
    lines = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in list(csv.reader(handle, delimiter="\t"))[1:]:
            cells = [cell.strip('"') for cell in row]
            if len(cells) < 14 or not cells[3]:
                continue
            labels = []
            for conf, term, j in (cells[8:11], cells[11:14]):
                core = re.search(r"\(2P\*<(\d)/2>\)\.(4[sp])$", conf)
                k = re.match(r"2\[(\d/2)\]", term)
                if core and k:
                    prime = "'" if core.group(1) == "1" else ""
                    labels.append(f"{core.group(2)}{prime}[{k.group(1)}]{j}")
            if len(labels) == 2:
                lines[f"{labels[1]}-{labels[0]}"] = float(cells[3]) / 1e6
    return lines


def checks() -> None:
    """Check the labels, the NIST column, Eq. (2) and Table 2 against Table 1."""
    levels = _levels()
    nist = _nist()
    found_misprints = {}
    for upper, lower, cells in _rows():
        pair = _identify(levels, cells[2])
        if pair != f"{upper}-{lower}":
            found_misprints[cells[2]] = (lower, pair.split("-")[1])
        paschen = cells[1].split("-")[-1]
        if paschen != PASCHEN[pair.split("-")[1]] or (
            "-" in cells[1] and cells[1].split("-")[0] != PASCHEN[upper]
        ):
            msg = f"{cells[2]}: Paschen {cells[1]!r} is not {pair}'s"
            raise ValueError(msg)
        printed = cells[5].strip()
        if printed in NO_VALUE:
            if pair in nist:
                msg = f"{cells[2]}: NIST column empty but NIST has {nist[pair]}"
                raise ValueError(msg)
            continue
        rounding = 0.5 * 10.0 ** -len([*printed.split("."), ""][1])
        if abs(float(printed) - nist[pair]) > rounding + 1e-9:
            msg = f"{cells[2]}: NIST column {printed} against {nist[pair]}"
            raise ValueError(msg)
    if found_misprints != MISPRINTED:
        msg = f"misprinted labels {found_misprints}, expected {MISPRINTED}"
        raise ValueError(msg)
    for points in _lifetimes().values():
        for point in points:
            percent = float(point.printed.split("(")[1].rstrip("%)"))
            ratio = 100 * float(point.sigma) / float(point.value)
            if abs(ratio - percent) > 0.06:
                msg = f"{point.row}: {point.printed} is not its own ratio {ratio:.2f}"
                raise ValueError(msg)
    _check_eq2()
    _check_table_2()


def _check_eq2() -> None:
    """Sum of A times the printed lifetime, per upper level."""
    files = extract()
    tau = {p.row: float(p.value) for p in files[FILE_CHANG] + files[FILE_ERMAN]}
    total: dict[str, float] = {}
    for point in files[FILE_A]:
        upper = point.x.split("-")[0]
        total[upper] = total.get(upper, 0.0) + float(point.value) * tau[upper] * 1e-3
    off = {upper for upper, value in total.items() if abs(value - 1) > 0.02}
    if off != EQ2_OFF:
        msg = f"Eq. (2) fails for {sorted(off)}, expected {sorted(EQ2_OFF)}: {total}"
        raise ValueError(msg)


def _check_table_2() -> None:
    """Table 2's S (Eq. 3) and f (Eq. 4) from Table 1's A, printed air wavelengths."""
    a = {p.row: float(p.value) * 1e6 for p in _transitions()}
    pairs = {cells[2]: (upper, lower) for upper, lower, cells in _rows()}
    off = set()
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        upper, lower = pairs[cells[2]]
        g_k, g_i = 2 * int(upper[-1]) + 1, 2 * int(lower[-1]) + 1
        angstrom = 10 * float(cells[2])
        value = a.get(cells[2], 1e6 * _nist().get(f"{upper}-{lower}", 0.0))
        expected = {
            "S": 4.935e-19 * angstrom**3 * g_k * value,
            "f": 1.499e-16 * g_k / g_i * angstrom**2 * value,
        }
        for column, cell in (("S", cells[3]), ("f", cells[5])):
            if cell.strip() in NO_VALUE:
                continue
            printed = float(cell)
            digits = 0.5 * 10.0 ** -len([*cell.split("."), ""][1].split("E")[0])
            scale = 10.0 ** int(cell.split("E")[1]) if "E" in cell else 1.0
            tolerance = max(0.03 * printed, digits * scale)
            if abs(expected[column] - printed) > tolerance:
                off.add((cells[2], column))
    if off != TABLE_2_OFF:
        msg = f"Table 2 rows off Eqs. (3), (4): {sorted(off)}"
        raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

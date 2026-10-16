# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Aymar and Coulombe (1978): calculated A values and lifetimes, Kr I and Xe I.

At. Data Nucl. Data Tables 21, 537-566 (1978); journal p. 537 is PDF page 1, so the
journal page is the PDF page plus 536. A scanned (raster) PDF: the cells come from the
library's ``.json`` (OCR), and every taken cell was checked against the rendered pages
(2026-10-12). Length and velocity forms (the paper's ``r`` and ``nabla``) are separate
files. Taken:

* **Table III** (Kr I, PDF p. 8, journal p. 544, ``tbl_f7938e7a2c11``): the 30
  transitions of the 4p5 5s - 4p5 5p, 5p' array, AKI (10^6 s^-1).
* **Table IV** (Xe I): the 30 transitions of the 5p5 6s - 6p, 6p' array (PDF p. 20,
  journal p. 556, ``tbl_d3ae3ba3383f``); and the 19 transitions from 6p' (and 6p[1/2]0)
  down to 5d (PDF p. 22, journal p. 558, ``tbl_a138bdfcb9ef``, the block of rows
  ``5D ... - 6P'``), AKI (10^6 s^-1).
* **Table VI** (Kr I, PDF p. 28, journal p. 564, ``tbl_b2d644c2e2b3``): the lifetimes
  (ns) of 5p, 5p', 6p and 6p'.
* **Table VIII** (Xe I, PDF p. 30, journal p. 566, ``tbl_fc287c6cfec5``): the lifetimes
  (ns) of 6p, 6p', 7p and 6d.

Rows are the paper's Jl labels, ``6p'[1/2]0`` (nl, a prime for the 2P1/2 core, [K], J);
a transition is ``lower-upper`` as Table III/IV prints it. The OCR drops or garbles the
primes and brackets of Table III/IV's labels (``5P [3/2]1`` for 5p'[3/2]1, ``6S*``), so
each level is identified by its printed energy (EI, EK) among the levels of Tables VI
and VIII and the lower levels of :data:`LOWER`, and the label is built from that.
Values are written as decimals (``.284+02`` -> ``28.4``); ``printed`` keeps the row.

**OCR readings corrected** (:data:`READINGS`; each extracted string is checked to be in
the row, and the page reading is what the rendered page prints): Table IV p. 556,
6s'[1/2]0 - 6p[3/2]1 length AKI ``.840+02`` printed ``.840-02``; 6s'[1/2]1 - 6p'[1/2]0
EK ``69861. 69861.`` printed ``89861.``; f of 6s[3/2]2 - 6p[5/2]3 ``.661+00`` printed
``.641+00``; p. 557, 7s[3/2]1 - 6p'[3/2]2 length AKI ``.843+04`` printed ``.843-04``; p.
558, three velocity FIK whose exponent sign the OCR lost; Table III p. 544, S of
5s'[1/2]0 - 5p[1/2]1 ``.428-00`` printed ``.428+00``. None is a misprint of the paper.

**Checks** (:func:`checks`; the script fails if one is missed): the counts; every
taken row's vacuum wavelength against its printed energies (1e8 / LAMBDA = EK - EI
within 1 cm^-1, the energies being rounded), its statistical weights against the
labels' J, and its AKI and FIK against the printed line strength S through Eqs. (1) and
(2) (A = 2.026e18 S / (g_k lambda^3), f = 303.7 S / (g_i lambda)) within 1 %, both
forms; each label's level among NIST's (``external/nist_asd/<X>_I_levels.tsv``:
configuration, term and J from the label, energy within 1 cm^-1 of the printed one);
and the lifetimes of Kr 5p, 5p' and Xe 6p, 6p' against 1 / sum A over every row of
Table III/IV the library holds with that upper level, within 1 % (:func:`_check_sums`;
for Kr 6p, 6p' and Xe 7p, 6d the library lost rows and the sums fall short).
"""

import csv
import functools
import re
from dataclasses import dataclass
from decimal import Decimal

from _json_table import json_table
from _library import DATA, Point, write

DIRECTORY = "Aymar_1978_TheoreticalTransitionProbabilities"
SOURCE = "aymar_1978"
TABLE_III = "tbl_f7938e7a2c11"  # Kr, PDF p. 8: 5s - 5p, 5p'
TABLE_IV_6S = "tbl_d3ae3ba3383f"  # Xe, PDF p. 20: 6s - 6p, 6p'
TABLE_IV_5D = "tbl_a138bdfcb9ef"  # Xe, PDF p. 22: 5d - 6p'
TABLE_VI = "tbl_b2d644c2e2b3"  # Kr lifetimes, PDF p. 28
TABLE_VIII = "tbl_fc287c6cfec5"  # Xe lifetimes, PDF p. 30
# Every library table of Table III (Kr) and Table IV (Xe), in page order, for the
# lifetime sums; tbl_381ebc0c445f (PDF p. 26) is the end of Table IV, though the
# library captions it Table V.
ALL_III = (
    "tbl_f7938e7a2c11",
    "tbl_d6b4e25cc275",
    "tbl_cdd39ba1cc6e",
    "tbl_4e509184057f",
    "tbl_3e1362a4ffb7",
    "tbl_f6143fc6821d",
    "tbl_5f045103902c",
    "tbl_37741e167e8e",
    "tbl_cb105b033686",
    "tbl_ff6a7a271710",
    "tbl_434c6fb74aec",
    "tbl_fb1ce226704d",
    "tbl_10dabab39091",
    "tbl_2d5c40e44559",
    "tbl_2d045db29324",
)
ALL_IV = (
    "tbl_d3ae3ba3383f",
    "tbl_bc8a50ae85db",
    "tbl_b6490f51ac03",
    "tbl_a138bdfcb9ef",
    "tbl_08bedb70f22f",
    "tbl_049fa2662c81",
    "tbl_351616df57a6",
    "tbl_740a3ec61948",
    "tbl_381ebc0c445f",
)
# Lower levels of the taken arrays that Tables VI and VIII do not list: (gas, rounded
# energy as Table III/IV prints it) -> label.
LOWER = {
    ("Kr", "79972"): "5s[3/2]2",
    ("Kr", "80917"): "5s[3/2]1",
    ("Kr", "85192"): "5s'[1/2]0",
    ("Kr", "85847"): "5s'[1/2]1",
    ("Xe", "67068"): "6s[3/2]2",
    ("Xe", "68046"): "6s[3/2]1",
    ("Xe", "76197"): "6s'[1/2]0",
    ("Xe", "77186"): "6s'[1/2]1",
}
CORE = {"Kr": "4s2.4p5", "Xe": "5p5"}
# OCR readings checked against the rendered page: (table, row) -> (extracted, page).
READINGS = {
    (TABLE_IV_6S, 21): (".840+02", ".840-02"),  # 6s'[1/2]0-6p[3/2]1, length AKI
    (TABLE_IV_6S, 32): ("69861. 69861.", "89861."),  # 6s'[1/2]1-6p'[1/2]0, EK
    (TABLE_IV_6S, 5): (".661+00", ".641+00"),  # 6s[3/2]2-6p[5/2]3, length FIK
    (TABLE_III, 21): (".428-00", ".428+00"),  # 5s'[1/2]0-5p[1/2]1, length S
    ("tbl_bc8a50ae85db", 12): (".843+04", ".843-04"),  # 7s[3/2]1-6p'[3/2]2, length AKI
    (TABLE_IV_5D, 37): (".485+03", ".485-03"),  # 5d[7/2]3-6p'[3/2]2, velocity FIK
    (TABLE_IV_5D, 38): (".563+02", ".563-02"),  # 5d[5/2]2-6p'[3/2]1, velocity FIK
    (TABLE_IV_5D, 40): (".153+04", ".153-04"),  # 5d[5/2]2-6p'[1/2]1, velocity FIK
}
_E = r"(-?[.,]\d{3}[+-]\d\d)"
_ROW = re.compile(
    r"(\d+[.,]\d)\s+(\d+)[.,]\s+(?:\d+[.,]\s+)?(\d+)[.,]\s+(\d)\s*(\d)\s+"
    + r"\s+".join([_E, _E, _E, r"(\S+)", _E, _E, _E, r"(\S+)"])
)
_LABEL = re.compile(r"^(\d[spd]) ?('?) ?\[(\d)/2\](\d)$")
FORMS = ("length", "velocity")


@dataclass(frozen=True)
class Row:
    """One transition of Table III/IV as printed (after the checked OCR readings)."""

    table: str
    index: int
    text: str
    wavelength: float  # vacuum, Angstrom
    lower: str  # printed energies, cm^-1, rounded
    upper: str
    g_lower: int
    g_upper: int
    rate: dict[str, str]  # AKI by form, as printed (.284+02)
    f: dict[str, str]
    s: dict[str, str]


def decimal(printed: str) -> str:
    """Return a Table III/IV number (``.284+02``) as a plain decimal (``28.4``)."""
    match = re.fullmatch(r"(-?)[.,](\d{3})([+-]\d\d)", printed)
    if match is None:
        msg = f"not a table number: {printed!r}"
        raise ValueError(msg)
    sign, digits, exponent = match.groups()
    value = Decimal(f"{sign}0.{digits}E{exponent}").normalize()
    return format(value, "f")


def _float(printed: str) -> float:
    return float(decimal(printed))


def rows(tables: tuple[str, ...]) -> list[Row]:
    """Return every parsed transition row of the given library tables."""
    found = []
    for table in tables:
        for index, cells in enumerate(json_table(DIRECTORY, table)):
            text = " ".join(cell for cell in cells if cell).replace("&#x27;", "'")
            if (table, index) in READINGS:
                extracted, page = READINGS[table, index]
                if extracted not in text:
                    msg = f"{table} row {index}: no {extracted!r} in {text!r}"
                    raise ValueError(msg)
                text = text.replace(extracted, page, 1)
            match = _ROW.search(text)
            if match is None:
                continue
            g = match.groups()
            found.append(
                Row(
                    table,
                    index,
                    text,
                    float(g[0].replace(",", ".")),
                    g[1],
                    g[2],
                    int(g[3]),
                    int(g[4]),
                    {"length": g[5], "velocity": g[9]},
                    {"length": g[6], "velocity": g[10]},
                    {"length": g[7], "velocity": g[11]},
                )
            )
    return found


def _jl(cell: str) -> str:
    """Return a Table VI/VIII Jl cell (``6p' [1/2]0``) as ``6p'[1/2]0``."""
    match = _LABEL.match(cell.replace("&#x27;", "'").strip())
    if match is None:
        msg = f"not a Jl label: {cell!r}"
        raise ValueError(msg)
    shell, prime, k, j = match.groups()
    return f"{shell}{prime}[{k}/2]{j}"


def lifetime_rows(table: str) -> list[tuple[str, str, str, str, str, str]]:
    """Return Table VI/VIII's (label, Paschen, energy, tau_r, tau_v, printed) rows."""
    grid = json_table(DIRECTORY, table)
    return [
        (_jl(cells[1]), cells[0], cells[2], cells[3], cells[4], " | ".join(cells[:5]))
        for cells in grid[2:]
    ]


@functools.cache
def levels(gas: str) -> dict[str, float]:
    """Return the gas's levels and printed energies: Table VI/VIII, and LOWER."""
    table = TABLE_VI if gas == "Kr" else TABLE_VIII
    found = {label: float(energy) for label, _, energy, *_ in lifetime_rows(table)}
    found |= {label: float(energy) for (g, energy), label in LOWER.items() if g == gas}
    return found


def label(gas: str, energy: str) -> str:
    """Return the level of a Table III/IV energy (rounded) within 0.51 cm^-1, or ""."""
    near = [k for k, e in levels(gas).items() if abs(e - float(energy)) < 0.51]
    if len(near) > 1:
        msg = f"{gas} {energy}: levels {near}"
        raise ValueError(msg)
    return near[0] if near else ""


def _taken(gas: str, table: str, lower: str, upper: tuple[str, ...]) -> list[Row]:
    """Return the rows of a table whose lower and upper labels start as given."""
    return [
        row
        for row in rows((table,))
        if label(gas, row.lower).startswith(lower)
        and label(gas, row.upper)[:2] in upper
    ]


ARRAYS = {
    # stem: (gas, table, lower shell prefix, upper shells)
    "5s5p": ("Kr", TABLE_III, "5s", ("5p",)),
    "6s6p": ("Xe", TABLE_IV_6S, "6s", ("6p",)),
    "5d6p": ("Xe", TABLE_IV_5D, "5d", ("6p",)),
}
LIFETIMES = {
    # gas: (table, the shells taken)
    "Kr": (TABLE_VI, ("5p", "6p")),
    "Xe": (TABLE_VIII, ("6p", "7p", "6d")),
}


def _transition(gas: str, row: Row) -> str:
    return f"{label(gas, row.lower)}-{label(gas, row.upper)}"


def extract() -> dict[str, list[Point]]:
    """Return the A values and the lifetimes, one file per gas, array and form."""
    files: dict[str, list[Point]] = {}
    for stem, (gas, table, lower, upper) in ARRAYS.items():
        for form in FORMS:
            key = f"{gas}/{SOURCE}__A_{form}_{gas}_{stem}.csv"
            files[key] = [
                Point(
                    _transition(gas, row),
                    _transition(gas, row),
                    decimal(row.rate[form]),
                    "",
                    row.text,
                )
                for row in _taken(gas, table, lower, upper)
            ]
    for gas, (table, shells) in LIFETIMES.items():
        taken = [r for r in lifetime_rows(table) if r[0][:2] in shells]
        for column, form in ((3, "length"), (4, "velocity")):
            key = f"{gas}/{SOURCE}__tau_{form}_{gas}.csv"
            files[key] = [Point(r[0], r[0], r[column], "", r[5]) for r in taken]
    return files


def _nist_levels(gas: str) -> list[dict[str, str]]:
    path = DATA / "external" / "nist_asd" / f"{gas}_I_levels.tsv"
    with path.open(encoding="utf-8") as handle:
        return [
            {k: (v or "").strip('"') for k, v in row.items() if k}
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def nist_level(gas: str, label: str) -> dict[str, str]:
    """Return NIST's level of a Jl label (``6p'[1/2]0``), by configuration, term, J."""
    match = _LABEL.match(label)
    if match is None:
        msg = f"not a label: {label!r}"
        raise ValueError(msg)
    shell, prime, k, j = match.groups()
    core = "1/2" if prime else "3/2"
    configuration = f"{CORE[gas]}.(2P*<{core}>).{shell}"
    term = f"2[{k}/2]" + ("*" if shell[1] in "sd" else "")
    found = [
        r
        for r in _nist_levels(gas)
        if r["Configuration"] == configuration and r["Term"] == term and r["J"] == j
    ]
    if len(found) != 1:
        msg = f"{gas} {label}: {len(found)} NIST levels"
        raise ValueError(msg)
    return found[0]


def _check_rows() -> None:
    """Check each taken row's wavelength, weights, A and f against S and the levels."""
    for stem, (gas, table, lower, upper) in ARRAYS.items():
        for row in _taken(gas, table, lower, upper):
            name = f"{stem} {_transition(gas, row)}"
            gap = float(row.upper) - float(row.lower)
            if abs(1e8 / row.wavelength - gap) > 1.0:
                msg = f"{name}: LAMBDA {row.wavelength} against EK - EI {gap}"
                raise ValueError(msg)
            for energy, g in ((row.lower, row.g_lower), (row.upper, row.g_upper)):
                level = nist_level(gas, label(gas, energy))
                if 2 * int(level["J"]) + 1 != g:
                    msg = f"{name}: g {g} against J {level['J']}"
                    raise ValueError(msg)
                if abs(float(level["Level (cm-1)"]) - float(energy)) > 1.0:
                    msg = f"{name}: E {energy} against NIST {level['Level (cm-1)']}"
                    raise ValueError(msg)
            for form in FORMS:
                s = _float(row.s[form])
                rate = 2.026e18 * s / (row.g_upper * row.wavelength**3) / 1e6
                f = 303.7 * s / (row.g_lower * row.wavelength)
                for what, computed, printed in (
                    ("A", rate, _float(row.rate[form])),
                    ("f", f, _float(row.f[form])),
                ):
                    if abs(computed / printed - 1) > 0.01:
                        msg = f"{name} {form} {what}: {printed} against {computed:.4g}"
                        raise ValueError(msg)


# The levels whose every decay the library's rows of Table III/IV hold, so that 1 / sum
# A reproduces the printed lifetime: Kr 5p, 5p' and Xe 6p, 6p'. For Kr 6p, 6p' and Xe
# 7p, 6d the library lost rows (truncated rows at the feet of pages, PDF p. 21's last
# rows) and the sums are not checked.
SUMMED = {"Kr": ("5p",), "Xe": ("6p",)}


def _check_sums() -> None:
    """Check each summable lifetime against 1 / sum A of the library's rows, 1 %."""
    for gas, (table, _) in LIFETIMES.items():
        found = rows(ALL_III if gas == "Kr" else ALL_IV)
        for name, _, energy, tau_r, tau_v, _ in lifetime_rows(table):
            if name[:2] not in SUMMED[gas]:
                continue
            mine = [r for r in found if abs(float(r.upper) - float(energy)) < 0.51]
            for form, printed in (("length", tau_r), ("velocity", tau_v)):
                tau = 1e3 / sum(_float(r.rate[form]) for r in mine)
                if abs(tau / float(printed) - 1) > 0.01:
                    msg = f"{gas} {name} {form}: 1/sum A {tau:.3g} against {printed}"
                    raise ValueError(msg)


def checks() -> None:
    """Check the counts, the rows, the NIST levels and the lifetime sums."""
    files = extract()
    counts = {key.split("__")[1]: len(points) for key, points in files.items()}
    expected = {
        "A_length_Kr_5s5p.csv": 30,
        "A_velocity_Kr_5s5p.csv": 30,
        "A_length_Xe_6s6p.csv": 30,
        "A_velocity_Xe_6s6p.csv": 30,
        "A_length_Xe_5d6p.csv": 19,
        "A_velocity_Xe_5d6p.csv": 19,
        "tau_length_Kr.csv": 20,
        "tau_velocity_Kr.csv": 20,
        "tau_length_Xe.csv": 24,
        "tau_velocity_Xe.csv": 24,
    }
    if counts != expected:
        msg = f"counts {counts}"
        raise ValueError(msg)
    _check_rows()
    for gas, (table, shells) in LIFETIMES.items():
        for label, _, energy, *_ in lifetime_rows(table):
            if label[:2] in shells:
                level = nist_level(gas, label)
                if abs(float(level["Level (cm-1)"]) - float(energy)) > 1.0:
                    msg = f"{gas} {label}: {energy} against NIST"
                    raise ValueError(msg)
    _check_sums()


def main() -> None:
    """Check the tables, and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

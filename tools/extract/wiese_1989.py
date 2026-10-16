# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wiese, Brault, Danzmann, Helbig and Kock (1989): Ar I transition rates.

Phys. Rev. A 39, 2461-2471 (1989) = PDF pages 1-11 (journal page = PDF page + 2460); a
scan, so every cell taken was checked against the table images (2026-10-09). Library
tables without a CSV rendering are read from their JSON cells.

- Table IV (p. 2467, ``tbl_c910a03b3fcc``): the recommended A (1e6 s-1) of the 30 Ar I
  4s-4p lines with their "Estimated uncertainty" (per row, relative); 29 values and the
  bound ``<0.001`` (857.806 nm), kept in the file and excluded from the dataset.
- Table VI (p. 2469, ``tbl_9316afb55907``): the recommended A of the 4s-5p lines, two
  standard deviations; the 28 rows with a value (397.972 nm has none, 395.7 nm was not
  observed). Footnote marks (``0.0383^{a}``) are dropped from the value.
- Table I (p. 2464, ``tbl_62e5805103d7``): the A of the 15 measured branches of
  5p[5/2]2 (the 16th, 6123.946 nm, has 100R < 0.06 and no A), with Delta R (%, two
  standard deviations) as each row's relative uncertainty.

Each row is the printed (air) wavelength, unique in each table; ``x`` is the transition
in the paper's Racah notation, lower-upper (``4s[3/2]2-4p[5/2]3``: nl[K]J, a prime for
the 2P1/2 core, the parity mark dropped), built from the printed cells. Labels the paper
misprints are kept as printed and corrected by errata in the source record.

**Checks** (:func:`checks`; the script fails if one is missed): every row's line is
identified by its printed wavelength among NIST's Ar I level pairs (air wavelengths
from NIST's levels, ``external/nist_asd/Ar_I_levels.tsv``), and the label after the
record's errata is that pair's, so the errata are exactly the misprinted labels; the
record's per-row relative uncertainties are the printed ones; Table I's A are 100R /
153 ns within the rounding (but 3045.369 nm, 0.148 printed against 0.1471), and its
100R sum to 100.01.
"""

import csv
import re
import tomllib

from _json_table import json_table
from _library import DATA, Point, number, write

DIRECTORY = "Wiese_1989_UnifiedSetAtomic"
SOURCE = "wiese_1989"
TABLE_I = "tbl_62e5805103d7"
TABLE_IV = "tbl_c910a03b3fcc"
TABLE_VI = "tbl_9316afb55907"
FILE_IV = "Ar/wiese_1989__A_4s_4p.csv"
FILE_VI = "Ar/wiese_1989__A_4s_5p.csv"
FILE_I = "Ar/wiese_1989__A_5p52_branches.csv"
UPPER_I = "5p[5/2]2"  # Table I's upper level (caption), 116 999.39 cm-1
LIFETIME_NS = 153.0  # the normalising lifetime of Table I (caption)
DATASETS = {FILE_IV: "A_4s_4p", FILE_VI: "A_4s_5p", FILE_I: "A_5p52_branches"}
_RACAH = re.compile(r"^(\d[spd])('?)\s*\[\s*\\frac\{(\d)\}\{(\d)\}\s*\]\s*_\{?(\d)\}?")
_MARK = re.compile(r"\^\{?[a-c]\}?$")


def racah(cell: str) -> str:
    """Return a printed level as ``nl'[K]J`` (``4p'[1/2]0``)."""
    match = _RACAH.match(cell.strip())
    if not match:
        msg = f"not a Racah label: {cell!r}"
        raise ValueError(msg)
    shell, prime, numerator, denominator, j = match.groups()
    return f"{shell}{prime}[{numerator}/{denominator}]{j}"


def _array(table_id: str) -> list[Point]:
    """Return the rows of Table IV or VI that have a value."""
    points = []
    lower = ""
    for cells in json_table(DIRECTORY, table_id)[1:]:
        lower = cells[0] or lower  # the lower level spans its group's rows
        value = _MARK.sub("", cells[3].strip())
        if not value:
            continue  # no A (Table VI: 397.972 nm, and 395.7 nm "not observed")
        value = value.removeprefix("<")  # the bound, excluded in the record
        label = f"{racah(lower)}-{racah(cells[1])}"
        points.append(Point(cells[2], label, number(value), "", " | ".join(cells)))
    return points


def _table_i() -> list[Point]:
    """Return Table I's branches with an A."""
    return [
        Point(
            cells[2],
            f"{racah(cells[0])}-{UPPER_I}",
            number(cells[5].replace(" ", "")),
            "",
            " | ".join(cells),
        )
        for cells in json_table(DIRECTORY, TABLE_I)[1:]
        if cells[5].strip()
    ]


def extract() -> dict[str, list[Point]]:
    """Return Tables IV, VI and I."""
    return {FILE_IV: _array(TABLE_IV), FILE_VI: _array(TABLE_VI), FILE_I: _table_i()}


def _air(vacuum_nm: float) -> float:
    """Return the air wavelength of a vacuum one (NIST ASD's dispersion formula)."""
    s2 = (1e3 / vacuum_nm) ** 2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return vacuum_nm / n


def _levels() -> dict[str, float]:
    """NIST's Ar I levels by Racah label, ``4p'[1/2]0`` -> cm-1."""
    path = DATA / "external" / "nist_asd" / "Ar_I_levels.tsv"
    levels = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in list(csv.reader(handle, delimiter="\t"))[1:]:
            conf, term, j, _, energy = (cell.strip('"') for cell in row[:5])
            core = re.search(r"\(2P\*<(\d)/2>\)\.(\d+[spd])$", conf)
            k = re.match(r"2\[(\d/2)\]", term)
            if core and k and energy:
                prime = "'" if core.group(1) == "1" else ""
                levels[f"{core.group(2)}{prime}[{k.group(1)}]{j}"] = float(energy)
    return levels


def _identify(levels: dict[str, float], wavelength: str, upper: str | None) -> str:
    """Return the lower-upper pair whose air wavelength is the printed one, alone."""
    printed = float(wavelength)
    tolerance = max(0.006, 0.5 * 10.0 ** -len(wavelength.split(".")[1]))
    if upper:  # Table I: the 5p[5/2]2 branches, 1929-4939 nm from older levels
        tolerance = 0.006 if printed < 1000 else 0.08
    found = []
    for lo, e_lo in levels.items():
        for up, e_up in levels.items():
            if (upper and up != upper) or e_up <= e_lo or up[1] != "p" or lo[1] == "p":
                continue
            if abs(_air(1e7 / (e_up - e_lo)) - printed) < tolerance:
                found.append(f"{lo}-{up}")
    if len(found) != 1:
        msg = f"{wavelength} nm: {found}"
        raise ValueError(msg)
    return found[0]


def checks() -> None:
    """Check the labels (after errata) against NIST, the uncertainties and Table I."""
    files = extract()
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    levels = _levels()
    for name, points in files.items():
        dataset = DATASETS[name]
        errata = {
            e["row"]: e for e in record.get("errata", []) if e["dataset"] == dataset
        }
        upper = UPPER_I if name == FILE_I else None
        for point in points:
            label = point.x
            if point.row in errata:
                if errata[point.row]["printed"] != label:
                    msg = f"{dataset} {point.row}: erratum of {label!r} misnamed"
                    raise ValueError(msg)
                label = errata.pop(point.row)["corrected"]
            if label != _identify(levels, point.row, upper):
                msg = f"{dataset} {point.row}: {label!r} is not the line's pair"
                raise ValueError(msg)
        if errata:
            msg = f"{dataset}: errata of no row {sorted(errata)}"
            raise ValueError(msg)
    _check_uncertainties(files, record)
    _check_table_i(files[FILE_I])


def _check_uncertainties(files: dict[str, list[Point]], record: dict) -> None:
    """Check that the record's relative_rows are the printed uncertainties."""
    entries = {d["id"]: d for d in record["datasets"]}
    for name, points in files.items():
        stated = entries[DATASETS[name]]["uncertainty"]["relative_rows"]
        printed = {}
        for point in points:
            cells = point.printed.split(" | ")
            cell = cells[4] if name == FILE_I else cells[-1]
            if cell:
                printed[point.row] = float(cell.strip("±%")) / 100
        if stated != printed:
            msg = f"{name}: relative_rows {stated} against the printed {printed}"
            raise ValueError(msg)


def _check_table_i(points: list[Point]) -> None:
    """Table I's A against 100R / 153 ns, and the sum of 100R."""
    rows = json_table(DIRECTORY, TABLE_I)[1:]
    total = sum(float(cells[3]) for cells in rows if not cells[3].startswith("<"))
    if f"{total:.2f}" != "100.01":
        msg = f"Table I: sum of 100R is {total}"
        raise ValueError(msg)
    for point in points:
        r100, value = point.printed.split(" | ")[3], point.value
        expected = float(r100) / 100 / (LIFETIME_NS * 1e-9) / 1e6
        rounding = 0.5 * 10.0 ** -len(value.split(".")[1])
        rounding += 0.5 * 10.0 ** -len(r100.split(".")[1]) / 100 / (LIFETIME_NS * 1e-3)
        off = abs(expected - float(value)) > rounding
        if off != (point.row == "3045.369"):
            msg = f"Table I {point.row}: A {value} against 100R/153 ns {expected:.5f}"
            raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

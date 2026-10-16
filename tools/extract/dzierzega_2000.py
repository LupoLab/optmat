# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Dzierżęga, Volz, Nave and Griesmann (2000): Kr I 5p-5s transition rates.

Phys. Rev. A 62, 022505 (2000), pages 022505-1 to -9 = PDF pages 1-9. Library tables
without a CSV rendering are read from their JSON cells; every cell taken is also found
verbatim in the PDF's text layer, whose font maps "±" to "6", "<" to "," and ">" to "."
(checked by :func:`checks`), and was checked against the table images (2026-10-09).

- Table IV (p. 5, ``tbl_31985c0969ea``): the absolute transition rates A_ki (1e6 s-1)
  of all 30 5p-5s lines with their uncertainties (``40.64±0.2``: value and sigma); the
  bound ``<0.04`` (557.31 nm) is kept in the file and excluded from the dataset. The
  branching fractions and the lifetime column are not taken (the lifetimes are Table
  II's and Schmitt et al.'s).
- Table II (p. 3, ``tbl_bf6340c2979f``): the semiempirical lifetimes of the ten 5p
  levels ("This work"), printed ``24.58±0.85``; the extraction reads the "±" as "6"
  (``24.5860.85``), and the value and its uncertainty both have two decimals.

Table IV's rows are the printed (air) wavelengths, unique; ``x`` is the transition as
printed, upper-lower in the paper's J_cK (Racah) notation (``5p'[1/2]0-5s[3/2]1``: a
prime for the 2P1/2 core), without the Paschen labels the cells add in parentheses
(they stay in ``printed``). Table II's rows are the Paschen labels, ``x`` the Racah
label decoded from the cell's font (``5p8@1/2#0`` is 5p'[1/2]0). A misprinted label is
kept as printed and corrected by an erratum in the source record.

**Checks** (:func:`checks`; the script fails if one is missed): the text layer holds
every taken cell; Table V's A column is Table IV's, row by row; Table IV's italic
(semiempirical) lifetimes of 2p1, 2p2, 2p5 and 2p10 are Table II's; the sum of A per
upper level is (1 - far-IR fraction) / tau within 0.2 %; every line is identified by
its printed wavelength among NIST's Kr I level pairs (air wavelengths from
``external/nist_asd/Kr_I_levels.tsv``), and the label after the record's errata is that
pair's, so the errata are exactly the misprinted labels.
"""

import csv
import re
import tomllib

import pypdfium2 as pdfium
from _figure import pdf
from _json_table import json_table
from _library import DATA, Point, number, write

DIRECTORY = "Dzierzega_2000_AccurateTransitionRates"
SOURCE = "dzierzega_2000"
TABLE_II = "tbl_bf6340c2979f"
TABLE_IV = "tbl_31985c0969ea"
TABLE_V = "tbl_7bdf8c4ef802"
FILE_A = "Kr/dzierzega_2000__A_5p_5s.csv"
FILE_TAU = "Kr/dzierzega_2000__lifetime_5p.csv"
DATASET_A = "A_5p_5s"
# Table II's far-IR (5p-4d) branches, % of the total decay rate (2p1-2p4).
FAR_IR = {"2p1": 0.246, "2p2": 0.034, "2p3": 0.106, "2p4": 0.0045}
SEMIEMPIRICAL = ("2p1", "2p2", "2p5", "2p10")  # Table IV's italic lifetimes
_LEVEL = re.compile(r"^(5[sp]'?\[\d/2\])_\{?(\d)\}?\s*\((\d[sp])_\{?(\d+)\}?\)$")
_FONT = re.compile(r"^5p(8?)@(\d/2)#(\d)$")
_TAU = re.compile(r"^(\d+\.\d\d)6(\d+\.\d\d)$")


def level(cell: str) -> tuple[str, str]:
    """Return a Table IV level as (Racah ``5p'[1/2]0``, Paschen ``2p1``)."""
    match = _LEVEL.match(cell.strip())
    if not match:
        msg = f"not a level: {cell!r}"
        raise ValueError(msg)
    racah, j, shell, index = match.groups()
    return f"{racah}{j}", f"{shell}{index}"


def _table_iv() -> list[tuple[str, list[str]]]:
    """Return Table IV's rows with the upper level's Paschen label carried down."""
    rows = []
    upper = ""
    for cells in json_table(DIRECTORY, TABLE_IV)[1:]:
        upper = cells[0] or upper
        rows.append((upper, cells))
    return rows


def _rates() -> list[Point]:
    points = []
    for upper, cells in _table_iv():
        value, _, sigma = cells[6].partition("±")
        label = f"{level(upper)[0]}-{level(cells[1])[0]}"
        points.append(
            Point(
                cells[2],
                label,
                number(value.removeprefix("<")),
                sigma.strip(),
                " | ".join(cells),
            )
        )
    return points


def _lifetimes() -> list[Point]:
    points = []
    for cells in json_table(DIRECTORY, TABLE_II)[1:]:
        font = _FONT.match(cells[1].split()[0])
        tau = _TAU.match(cells[8])
        if not (font and tau):
            msg = f"Table II row {cells!r}"
            raise ValueError(msg)
        prime, k, j = font.groups()
        racah = f"5p{chr(39) if prime else ''}[{k}]{j}"
        points.append(
            Point(cells[0], racah, tau.group(1), tau.group(2), " | ".join(cells))
        )
    return points


def extract() -> dict[str, list[Point]]:
    """Return Table IV's rates and Table II's semiempirical lifetimes."""
    return {FILE_A: _rates(), FILE_TAU: _lifetimes()}


def _air(vacuum_nm: float) -> float:
    """Return the air wavelength of a vacuum one (NIST ASD's dispersion formula)."""
    s2 = (1e3 / vacuum_nm) ** 2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return vacuum_nm / n


def _levels() -> dict[str, float]:
    """NIST's Kr I 5s and 5p levels by Racah label, ``5p'[1/2]0`` -> cm-1."""
    path = DATA / "external" / "nist_asd" / "Kr_I_levels.tsv"
    levels = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in list(csv.reader(handle, delimiter="\t"))[1:]:
            conf, term, j, _, energy = (cell.strip('"') for cell in row[:5])
            core = re.search(r"\(2P\*<(\d)/2>\)\.(5[sp])$", conf)
            k = re.match(r"2\[(\d/2)\]", term)
            if core and k and energy:
                prime = "'" if core.group(1) == "1" else ""
                levels[f"{core.group(2)}{prime}[{k.group(1)}]{j}"] = float(energy)
    return levels


def _identify(levels: dict[str, float], wavelength: str) -> str:
    """Return the upper-lower pair whose air wavelength is the printed one, alone."""
    found = [
        f"{up}-{lo}"
        for up, e_up in levels.items()
        for lo, e_lo in levels.items()
        if up.startswith("5p")
        and lo.startswith("5s")
        and abs(_air(1e7 / (e_up - e_lo)) - float(wavelength)) < 0.006
    ]
    if len(found) != 1:
        msg = f"{wavelength} nm: {found}"
        raise ValueError(msg)
    return found[0]


def _text(page: int) -> str:
    """Return the text layer of a page (1-based)."""
    document = pdfium.PdfDocument(pdf(DIRECTORY))
    return document[page - 1].get_textpage().get_text_range()


def checks() -> None:
    """Check the cells, Tables V and II, the sums and the labels against NIST."""
    files = extract()
    layer_iv = _text(5)
    layer_iv = layer_iv[layer_iv.index("TABLE IV") :]
    for _, cells in _table_iv():
        for cell in cells[2:]:
            token = cell.replace("±", "6").replace("<10-3", ",1023")
            token = token.replace("<", ",").replace(">", ".")
            if token not in layer_iv:
                msg = f"Table IV cell {cell!r} not in the text layer"
                raise ValueError(msg)
    layer_ii = _text(3)
    for point in files[FILE_TAU]:
        if f"{point.value}6{point.sigma}" not in layer_ii:
            msg = f"Table II {point.row} not in the text layer"
            raise ValueError(msg)
    _check_table_v(files[FILE_A])
    _check_sums(files)
    _check_labels(files[FILE_A])


def _check_table_v(points: list[Point]) -> None:
    """Table V's A column ("From Table IV") is Table IV's, uncertainty included."""
    for point, cells in zip(points, json_table(DIRECTORY, TABLE_V)[1:], strict=True):
        cell = cells[3]
        if cells[2] != point.row:
            msg = f"Table V row {cells[2]} against {point.row}"
            raise ValueError(msg)
        if "(" in cell:
            value, uncertainty = cell.rstrip(")").split("(")
            decimals = len(value.split(".")[1])
            same = (
                abs(float(value) - float(point.value)) < 1e-9
                and abs(int(uncertainty) * 10.0**-decimals - float(point.sigma)) < 1e-9
            )
        else:
            same = cell == "<" + point.value
        if not same:
            msg = f"Table V {point.row}: {cell!r} against Table IV {point.value}"
            raise ValueError(msg)


def _check_sums(files: dict[str, list[Point]]) -> None:
    """Sum of A per level against (1 - far IR) / tau; Table IV's italic lifetimes."""
    tau = {p.row: p for p in files[FILE_TAU]}
    sums: dict[str, float] = {}
    printed_tau: dict[str, str] = {}
    for (upper, cells), point in zip(_table_iv(), files[FILE_A], strict=True):
        paschen = level(upper)[1]
        printed_tau[paschen] = cells[3] or printed_tau.get(paschen, "")
        if not cells[6].startswith("<"):
            sums[paschen] = sums.get(paschen, 0.0) + float(point.value)
    for paschen, total in sums.items():
        value = float(printed_tau[paschen].partition("±")[0])
        expected = (1 - FAR_IR.get(paschen, 0.0) / 100) / (value * 1e-3)
        if abs(total / expected - 1) > 0.002:
            msg = f"{paschen}: sum of A {total:.3f} against {expected:.3f}"
            raise ValueError(msg)
    for paschen in SEMIEMPIRICAL:
        if printed_tau[paschen] != f"{tau[paschen].value}±{tau[paschen].sigma}":
            msg = f"{paschen}: Table IV's {printed_tau[paschen]!r} is not Table II's"
            raise ValueError(msg)


def _check_labels(points: list[Point]) -> None:
    """Every label, after the record's errata, is NIST's pair for the wavelength."""
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    errata = {
        e["row"]: e for e in record.get("errata", []) if e["dataset"] == DATASET_A
    }
    levels = _levels()
    for point in points:
        label = point.x
        if point.row in errata:
            if errata[point.row]["printed"] != label:
                msg = f"{point.row}: erratum of {label!r} misnamed"
                raise ValueError(msg)
            label = errata.pop(point.row)["corrected"]
        if label != _identify(levels, point.row):
            msg = f"{point.row}: {label!r} is not the line's pair"
            raise ValueError(msg)
    if errata:
        msg = f"errata of no row {sorted(errata)}"
        raise ValueError(msg)


def main() -> None:
    """Check the tables and write the data files."""
    checks()
    write(extract())


if __name__ == "__main__":
    main()

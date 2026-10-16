# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Sato et al. (2011): He two-photon ionisation cross sections at the SCSS.

The paper is born digital (text layer on every page; its three figures are embedded
images). Its numbers are in two tables and the text:

* **Table 2** (p. 3, ``tbl_bc8ebb788440``, image checked): the two-photon ionisation
  cross section of He at 61.4, 58.4, 56.0 and 53.4 nm at 1e13 W/cm2, in 10^-52 cm^4
  s: the measurement ('Obs.', this work) with its error in parentheses on the last
  digit ('2.6 (6)'; none for 56.0 nm), and two columns of other papers' theory: the
  R-matrix Floquet values of van der Hart and Bingham (2005), Ref. 11, and the LOPT
  values of Nikolopoulos and Lambropoulos (2001), Ref. 12 ('-' at 53.4 nm: none). The
  row of H2 one-photon cross sections (Ref. 14, Chung et al. 1993) is not taken. The
  library's row labels are garbled (sigma 'H(2e)' and 'Obs.a'); the image reads
  sigma_He^(2) (10^-52 cm4 s) Obs.^a, Theory^b, Theory^c, sigma_H2 (10^-18 cm2)^d.
* **Table 1** (p. 3, ``tbl_076bec49a577``): the exponents n1 (He+) and n2 (H2+) of the
  ion yields' intensity dependence, Eqs. 1-2; not data of any quantity, read for the
  checks.
* **the text** (p. 5): the cross sections 'at 1 x 10^10 W cm^-2' extrapolated along
  Fig. 2's lines, 5.1, 550 and 2.2 x 10^-52 cm^4 s at 61.4, 58.4 and 56.0 nm, recorded
  inline in the source record (``text_1e10``).

**Checks** (:func:`checks`): the record's inline values against the sentence of the
text (the PDF's text layer, p. 6 of the PDF); the text's LOPT values (8.0, 360, 2.7)
against Table 2's Theory^c column; the extrapolation from 1e13 to 1e10 W/cm2 by Eq. 5,
sigma proportional to I^(n1 - n2 - 1) with Table 1's exponents, against the text's
values at 61.4 and 56.0 nm within :data:`EXTRAPOLATION` (at 58.4 nm it is reported,
not checked: Table 1's exponents give about 2500 x 10^-52 against the printed 550, the
register's finding). :func:`intensities` reports the peak intensity a flat-top
estimate gives from the stated pulse energies at the focus, the 10 um spot and the
300 fs pulse, against the text's central intensities.
"""

import math
import re
import tomllib
from typing import Any

import pypdfium2 as pdfium
from _figure import pdf
from _library import DATA, Point, number, table, write

DIRECTORY = "Sato_2011_DeterminationAbsoluteTwo"
SOURCE = "sato_2011"
TABLE_1 = "tbl_076bec49a577"
TABLE_2 = "tbl_bc8ebb788440"
ROWS = {
    "He/sato_2011__tbl2_obs_He.csv": "Obs",
    "He/sato_2011__tbl2_rmf_He.csv": "Theoryb",
    "He/sato_2011__tbl2_lopt_He.csv": "Theoryc",
}
"""Each data file and the start of its Table 2 row label."""
DASHES = ("-", "\N{EN DASH}")
"""Table 2's '-' (no value; an en dash in the image, a hyphen in the library's table):
the LOPT column at 53.4 nm."""
_CELL = re.compile(r"^(?P<value>[\d.]+)(?:\s*\((?P<error>\d+)\))?$")
TABLE_INTENSITY = 1e13
"""Table 2's intensity, W/cm2 (its caption and p. 3)."""
TEXT_INTENSITY = 1e10
"""The intensity the text extrapolates to, W/cm2 (p. 5)."""
EXTRAPOLATION = 0.10
"""The text's extrapolated values are matched within this (relative): Table 2's 0.9 has
one figure (5.6 %), and 0.01 in an exponent is 7 % over three decades of intensity."""


def _cell(cell: str) -> tuple[str, str]:
    """Return a cell's value and its error, both as printed digits: '2.6 (6)' -> 0.6.

    The error in parentheses is on the value's last digit(s); '' where none is printed.
    """
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    value, error = match["value"], match["error"]
    if error is None:
        return number(value), ""
    decimals = len(value.split(".")[1]) if "." in value else 0
    return number(value), f"{int(error) * 10.0**-decimals:.{decimals}f}"


def _table(table_id: str) -> tuple[list[str], dict[str, list[str]]]:
    """Return a table's wavelength header (nm, as printed) and its rows by label."""
    rows = table(DIRECTORY, table_id)
    header = rows[0]
    if not header[0].startswith("λ (nm)"):
        msg = f"{table_id}: unexpected header {header}"
        raise ValueError(msg)
    return header[1:], {row[0]: row[1:] for row in rows[1:]}


def _row(label: str) -> list[str]:
    """Return the cells of the Table 2 row whose label starts with (or contains) one."""
    _, rows = _table(TABLE_2)
    found = [cells for name, cells in rows.items() if name.startswith(label)]
    found = found or [cells for name, cells in rows.items() if label in name]
    if len(found) != 1:
        msg = f"Table 2: {len(found)} rows labelled {label!r}"
        raise ValueError(msg)
    return found[0]


def extract() -> dict[str, list[Point]]:
    """Return Table 2's cross-section rows, after the checks."""
    checks()
    wavelengths, _ = _table(TABLE_2)
    out: dict[str, list[Point]] = {}
    for path, label in ROWS.items():
        points = []
        for wavelength, cell in zip(wavelengths, _row(label), strict=True):
            if cell.strip() in DASHES:
                continue
            value, error = _cell(cell)
            points.append(Point(wavelength, wavelength, value, error, cell))
        out[path] = points
    return out


def exponents() -> dict[str, tuple[float, float]]:
    """Return Table 1's exponents n1 (He+) and n2 (H2+) by wavelength (as printed)."""
    wavelengths, rows = _table(TABLE_1)
    (n1,) = [cells for name, cells in rows.items() if name.startswith("n_{1}")]
    (n2,) = [cells for name, cells in rows.items() if name.startswith("n_{2}")]
    return {
        w: (float(_cell(a)[0]), float(_cell(b)[0]))
        for w, a, b in zip(wavelengths, n1, n2, strict=True)
    }


def _text() -> str:
    """Return the text layer of the paper's journal p. 5 (PDF p. 6), one line."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    text = document[5].get_textpage().get_text_range()
    return re.sub(r"\s+", " ", text.replace("\N{MINUS SIGN}", "-"))


TIMES = "\N{MULTIPLICATION SIGN}"
"""The text layer's multiplication sign."""
_EXTRAPOLATED = re.compile(
    rf"1 {TIMES} 1010 W cm-2 are estimated to be (?P<a>[\d.]+) {TIMES} 10-52, "
    rf"(?P<b>[\d.]+) {TIMES} 10-52 and (?P<c>[\d.]+) {TIMES} 10-52 cm4s at 61\.4, "
    r"58\.4, and 56\.0 nm"
)
_LOPT = re.compile(
    rf"consistent with the cross section data, (?P<a>[\d.]+) {TIMES} 10-52, "
    rf"(?P<b>[\d.]+) {TIMES} 10-52 and (?P<c>[\d.]+) {TIMES} 10-52 cm4s, reported "
    r"theoretically by Nikolopoulos"
)


def _sentence(pattern: re.Pattern[str]) -> dict[str, str]:
    """Return the three values of a sentence of p. 5, keyed 61.4, 58.4, 56.0 nm."""
    match = pattern.search(_text())
    if match is None:
        msg = f"the sentence {pattern.pattern[:40]!r}... is not on p. 5"
        raise ValueError(msg)
    return dict(
        zip(("61.4", "58.4", "56.0"), (match["a"], match["b"], match["c"]), strict=True)
    )


def _record() -> dict[str, dict[str, Any]]:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        return {d["id"]: d for d in tomllib.load(handle)["datasets"]}


def extrapolated() -> dict[str, tuple[float, float]]:
    """Return, by wavelength: Eq. 5's value at 1e10 W/cm2 and the text's (1e-52 cm4 s).

    sigma(I) = sigma(1e13) (I / 1e13)^(n1 - n2 - 1), from Table 2's measurement and
    Table 1's exponents (Eq. 5 with Eqs. 1-2: sigma is proportional to Y1 / (Y2 I)).
    """
    observed = {
        w: float(_cell(c)[0])
        for w, c in zip(_table(TABLE_2)[0], _row("Obs"), strict=True)
    }
    printed = _sentence(_EXTRAPOLATED)
    ratio = TEXT_INTENSITY / TABLE_INTENSITY
    return {
        w: (observed[w] * ratio ** (n1 - n2 - 1), float(printed[w]))
        for w, (n1, n2) in exponents().items()
        if w in printed
    }


def checks() -> dict[str, object]:
    """Return the checks against the text; raise if one fails."""
    out: dict[str, object] = {}
    # The record's inline values are the text's.
    (point_set,) = [_record()["text_1e10"]["points"]]
    recorded = {row: value for row, _, value, *_ in point_set}
    printed = _sentence(_EXTRAPOLATED)
    if recorded != printed:
        msg = f"the record's text_1e10 {recorded} is not the text's {printed}"
        raise ValueError(msg)
    out["text at 1e10 W/cm2 (1e-52 cm4 s)"] = printed
    # The text's LOPT values are Table 2's Theory^c column.
    wavelengths, _ = _table(TABLE_2)
    column = dict(zip(wavelengths, _row("Theoryc"), strict=True))
    lopt = _sentence(_LOPT)
    for wavelength, value in lopt.items():
        if float(value) != float(column[wavelength]):
            msg = f"{wavelength} nm: text LOPT {value}, Table 2 {column[wavelength]}"
            raise ValueError(msg)
    out["text's LOPT values = Table 2 Theory^c"] = lopt
    # Eq. 5's extrapolation with Table 1's exponents against the text's values.
    for wavelength, (computed, text) in extrapolated().items():
        out[f"{wavelength} nm at 1e10: Eq. 5 with Table 1, text"] = (
            round(computed, 2),
            text,
            round(computed / text, 3),
        )
        if wavelength != "58.4" and abs(computed / text - 1) > EXTRAPOLATION:
            msg = f"{wavelength} nm: Eq. 5 gives {computed:.3g}, the text {text}"
            raise ValueError(msg)
    # The exponent the text's 58.4 nm value implies, against Table 1's n1 - n2 - 1.
    observed = float(_cell(_row("Obs")[wavelengths.index("58.4")])[0])
    implied = math.log10(550 / observed) / math.log10(TEXT_INTENSITY / TABLE_INTENSITY)
    n1, n2 = exponents()["58.4"]
    out["58.4 nm exponent: implied by the text, Table 1's n1 - n2 - 1"] = (
        round(implied, 3),
        round(n1 - n2 - 1, 3),
    )
    return out


# The experiment's stated pulse energies at the focus (p. 2), uJ, the focal spot's
# diameter (10 um, p. 2) and the assumed pulse length (300 fs, p. 3), and the central
# intensities of the shot-to-shot fluctuation (p. 4), W/cm2.
PULSE_UJ = {"53.4": 1.3, "56.0": 2.7, "58.4": 4.8, "61.4": 6.7}
SPOT_CM = 10e-4
PULSE_S = 300e-15
CENTRAL = {"61.4": 4e13, "58.4": 3e13, "56.0": 2e13, "53.4": 1.3e13}


def intensities() -> dict[str, tuple[float, float]]:
    """Return, by wavelength: a flat-top intensity from the pulse energy; the text's.

    E / (tau pi (d / 2)^2), W/cm2; the paper does not say how it obtains its intensity.
    """
    area = math.pi * (SPOT_CM / 2) ** 2
    return {
        w: (energy * 1e-6 / (PULSE_S * area), CENTRAL[w])
        for w, energy in PULSE_UJ.items()
    }


def main() -> None:
    """Write the data files, and print the checks and the intensity estimate."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for wavelength, (flat, central) in intensities().items():
        print(
            f"{wavelength} nm: flat-top {flat:.2e} W/cm2, text {central:.1e} "
            f"({central / flat:.2f})"
        )


if __name__ == "__main__":
    main()

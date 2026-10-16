# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract SCHOTT's N-BK7 data sheet (as of 01-Dec-2023): indices, transmittance, dn/dT.

The data sheet is in the library as a document without a DOI. Its text layer (the PDF
is made from a spreadsheet), rendered with ``pdftotext -layout`` and checked against the
page image, is attached to it as ``n-bk7-datasheet.txt``. The rows are read from their
fixed positions in that rendering:

- "Refractive Indices": label, wavelength (nm) and index, the first three fields of
  lines 7-35; the rows at 296.7, 280.4 and 248.3 nm print no index and are skipped;
- "Internal Transmittance": wavelength (nm), tau_i of 10 mm and of 25 mm, the fields
  from column 58 on; empty cells are skipped;
- "Temperature Coefficients of the Refractive Index": for each interval (-40/-20,
  +20/+40, +60/+80), dn_rel/dT then dn_abs/dT (1e-6/K) at 1060.0 nm, e and g. The
  lines e and g are placed at the data sheet's own wavelengths for them (546.1 and
  435.8 nm).
"""

import re

from _library import Point, attachment, number, write

PAPER = "Unverified_b3287ae9ef98"
TEXT = "n-bk7-datasheet.txt"
PREFIX = "N-BK7/schott_2023__"

_INDEX = re.compile(r"^\s{3}(n\S*)\s+(\d+\.\d)(?:\s+(1\.\d{5}))?\s{2,}(.*)$")
_TAU = re.compile(r"^\s*(\d{3,4})(?:\s+(0\.\d+))?(?:\s+(0\.\d+))?(?:\s{2,}|$)")
_DNDT = re.compile(r"^\s+([-+]\d+/[-+]\d+)((?:\s+\d\.\d){6})")
# The columns of the dn/dT table, as printed, and where they sit.
_LINES = (("1060.0", "1060.0"), ("e", "546.1"), ("g", "435.8"))
_INTERVALS = {"-40/-20": "m30C", "+20/+40": "30C", "+60/+80": "70C"}


def _transmittance_part(line: str) -> str:
    """Return the part of a line from the transmittance column on (column 58)."""
    return line[58:]


def extract() -> dict[str, list[Point]]:
    """Return the data files and their points."""
    lines = attachment(PAPER, TEXT).read_text(encoding="utf-8").splitlines()
    index: list[Point] = []
    tau = {"10mm": [], "25mm": []}
    dndt: dict[str, list[Point]] = {}
    for line in lines[6:35]:
        match = _INDEX.match(line)
        if match:
            label, wavelength, value, _ = match.groups()
            if value is not None:
                raw = " ".join((label, wavelength, value))
                index.append(Point(label, number(wavelength), number(value), "", raw))
    for line in lines[6:50]:
        match = _TAU.match(_transmittance_part(line))
        if match and match.group(1):
            wavelength, ten, twenty_five = match.groups()
            raw = " ".join(v for v in (wavelength, ten, twenty_five) if v)
            if ten:
                tau["10mm"].append(Point(wavelength, wavelength, number(ten), "", raw))
            if twenty_five:
                tau["25mm"].append(
                    Point(wavelength, wavelength, number(twenty_five), "", raw)
                )
    for line in lines[70:88]:
        match = _DNDT.match(line)
        if match:
            interval, cells = match.group(1), match.group(2).split()
            name = _INTERVALS[interval]
            for kind, values in (("rel", cells[:3]), ("abs", cells[3:])):
                dndt[f"{PREFIX}dndT_{kind}_{name}.csv"] = [
                    Point(label, x, number(v), "", f"{interval} {label} {v}")
                    for (label, x), v in zip(_LINES, values, strict=True)
                ]
    if len(index) != 20 or len(tau["10mm"]) != 26 or len(tau["25mm"]) != 25:
        msg = f"unexpected rows: {len(index)}, {len(tau['10mm'])}, {len(tau['25mm'])}"
        raise ValueError(msg)
    if len(dndt) != 6:
        msg = f"dn/dT: {len(dndt)} rows of 6"
        raise ValueError(msg)
    return {
        f"{PREFIX}n_catalogue.csv": index,
        f"{PREFIX}tau_10mm.csv": tau["10mm"],
        f"{PREFIX}tau_25mm.csv": tau["25mm"],
        **dndt,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

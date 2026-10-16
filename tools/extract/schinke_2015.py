# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Schinke et al. (2015): band-to-band alpha of silicon at 295 K (Table I).

Table I (Appendix A, PDF pp. 20-21; ``tbl_80266e561237`` and its continuation
``tbl_8c2f326870f4``) prints, in three column groups per page, the wavelength (nm),
alpha_bb (cm^-1) and its relative expanded uncertainty U(alpha_bb)/alpha_bb (%, coverage
factor k = 2, "rounded to two significant digits"): 250-1450 nm at 10 nm steps, 121
rows. Each alpha cell prints its own power of ten, ``m \times 10^{p}``; the library's
CSV loses the superscript ("1.804 x 106" for 1.804 x 10^6, "9.503 x 10-1", with the
multiplication sign and no spaces), and since every cell is the mantissa, the sign and
"10", the exponent is the rest of the cell. It is stored in e-notation (``1.804e6``,
scale 1); the cell stays in ``printed``.

Both parts were checked against their images and every cell against the PDF's text
layer (``pdftotext -layout``), on 2026-10-14: they agree. The column groups are read
group by group and the rows sorted by wavelength; the grid is checked to be complete.

The relative uncertainty has no column in the format (``sigma`` is absolute, in the
printed unit), so the record carries it as ``relative_rows`` (U/100, level "2 sigma").
:func:`extract` checks that the record's ``relative_rows`` are this table's printed
percentages divided by 100, row by row, so the record cannot drift from the table.

Table II (Appendix B, the normalised temperature coefficient c_T = (1/alpha)
d alpha/dT at 295 K, ``tbl_75b95580c5bd`` and ``tbl_b00ddc0c0d86``) has no quantity
and is not extracted.
"""

import re
import tomllib
from decimal import Decimal

from _library import DATA, MINUS, Point, number, table, write

DIRECTORY = "Schinke_2015_UncertaintyAnalysisCoefficient"
PARTS = ("tbl_80266e561237", "tbl_8c2f326870f4")
ALPHA = "\N{GREEK SMALL LETTER ALPHA}"
TIMES = "\N{MULTIPLICATION SIGN}"
GROUP = [
    "\N{GREEK SMALL LETTER LAMDA} [nm]",  # codespell:ignore lamda (the Unicode name)
    f"{ALPHA}bb [cm-1]",
    f"U{ALPHA}({ALPHA}bbbb) [%]",
]
FILE = "Si/schinke_2015__alpha_295K.csv"
ROWS = 121

_POWER = re.compile(rf"^(?P<mantissa>\d\.\d{{3}}){TIMES}10(?P<power>-?\d)$")


def scientific(cell: str) -> str:
    """Return a cell "m x 10p" (the superscript lost) as ``mep``."""
    match = _POWER.match(cell.strip().replace(MINUS, "-"))
    if match is None:
        msg = f"not a power-of-ten cell: {cell!r}"
        raise ValueError(msg)
    value = f"{match['mantissa']}e{match['power']}"
    float(value)
    return value


def rows() -> list[list[str]]:
    """Return Table I's (lambda, alpha, U) triples, sorted by wavelength."""
    found: list[list[str]] = []
    for table_id in PARTS:
        cells = table(DIRECTORY, table_id)
        if cells[0] != GROUP * 3:
            msg = f"{table_id}: unexpected header {cells[0]}"
            raise ValueError(msg)
        for group in range(3):
            for line in cells[1:]:
                triple = line[3 * group : 3 * group + 3]
                if any(triple):
                    found.append(triple)
    found.sort(key=lambda triple: int(triple[0]))
    if [int(t[0]) for t in found] != list(range(250, 1451, 10)) or len(found) != ROWS:
        msg = f"Table I: {len(found)} rows, not the 250-1450 nm grid"
        raise ValueError(msg)
    return found


def _check_record(found: list[list[str]]) -> None:
    """Check that the record's relative_rows are the table's U / 100."""
    with (DATA / "sources" / "schinke_2015.toml").open("rb") as handle:
        record = tomllib.load(handle)
    (dataset,) = (d for d in record["datasets"] if d.get("file") == FILE)
    recorded = dataset["uncertainty"]["relative_rows"]
    expected = {t[0]: float(Decimal(number(t[2])) / 100) for t in found}
    if recorded != expected:
        wrong = sorted(set(recorded.items()) ^ set(expected.items()))
        msg = f"schinke_2015: relative_rows differ from Table I: {wrong}"
        raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Return alpha_bb at 295 K, checking the record's relative uncertainties."""
    found = rows()
    _check_record(found)
    points = [Point(t[0], t[0], scientific(t[1]), "", " | ".join(t)) for t in found]
    return {FILE: points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()

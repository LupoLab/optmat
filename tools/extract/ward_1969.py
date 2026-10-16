# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ward & New (1969): third-harmonic coefficients chi at 6943 A (Table III).

Table III (``tbl_ffb0f106d17e``, published only as JSON cells) compares, in
10^-39 esu/atom, the two experimental sets (i) and (ii) of the third-harmonic
coefficient chi(-3w; w, w, w), "chi theoretical", "chi derived" (from the Kerr or dc
coefficients by Eq. (IV.4)), omega_n, chi^K and chi^0. The cells were compared with
the table image; they agree. What the grid does:

- helium's absolute experimental values print as ``7.0\times 3^{\pm 1}`` and
  ``24\times 3^{\pm 1}`` (a factor-of-3 uncertainty, not a sigma): the value is taken
  and the factor is left to the record's notes;
- helium's three theoretical values share one cell, ``4.0^{a}12.0^{b}0.9^{c}`` (stacked
  in the image): split by their footnote letters;
- the four "chi derived" values of helium sit in rows below it, each on the row of its
  input (chi^K or chi^0); Ne and Ar have a second derived row (from chi^0);
- other cells are ``value\pm sigma`` with an optional footnote ``^{x}``.

Every other gas is "relative to chi_He = 4.0" (the image's bracket), as the caption
says. Table II's relative X (arbitrary units, the same values with smaller
uncertainties), its S_H/S_0, X/Delta k and literature Delta k are not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Ward_1969_OpticalThirdHarmonic"
TABLE = "tbl_ffb0f106d17e"
RUBY = "694.3"  # nm, 6943 A as printed (abstract, p. 57)
HE_NE = "632.8"  # nm, Buckingham and Dunmur's 6328 A laser (their p. 1777)
GASES = {"He": "He", "Ne": "Ne", "Ar": "Ar", "Kr": "Kr", "Xe": "Xe"}
GASES |= {"H_2": "H2", "CO_2": "CO2", "N_2": "N2"}

_VALUE = re.compile(
    r"^(?P<value>[\d.]+)(?:\\pm\s*(?P<sigma>[\d.]+))?(?:\^\{(?P<mark>[a-f])\})?$"
)
_FACTOR = re.compile(r"^(?P<value>[\d.]+)\\times 3\^\{\\pm 1\}$")
_STACKED = re.compile(r"([\d.]+)\^\{([a-f])\}")

# chi^K and chi^0 footnotes -> dataset prefix.
KERR = {
    "a": "kerr_sitz_yaris",
    "e": "kerr_buckingham_dunmur",
    "f": "kerr_buckingham_orr",
}
STATIC = {"a": "static_sitz_yaris", "d": "static_hurst"}
THEORY = {"a": "sitz_yaris", "b": "appendix_b", "c": "dawes"}
# The input each "chi derived" row is derived from, by the footnote of its chi^K or
# chi^0 cell on the same row.
DERIVED_ROW = {
    ("K", "a"): "from_kerr_sitz_yaris",
    ("K", "e"): "from_kerr_buckingham_dunmur",
    ("K", "f"): "from_kerr_buckingham_orr",
    ("0", "a"): "from_static_sitz_yaris",
    ("0", "d"): "from_static_hurst",
}


def _cell(text: str) -> tuple[str, str, str]:
    """Return the value, sigma and footnote letter of a Table III cell."""
    match = _VALUE.match(text.replace(" ", ""))
    if match is None:
        msg = f"unexpected cell {text!r}"
        raise ValueError(msg)
    return (
        number(match["value"]),
        match["sigma"] or "",
        match["mark"] or "",
    )


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column."""
    files: dict[str, list[Point]] = {}

    def add(name: str, gas: str, point: Point) -> None:
        files.setdefault(f"{gas}/ward_1969__{name}_{gas}.csv", []).append(point)

    gas = ""
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        if cells[0]:
            gas = GASES[cells[0]]
        if cells[1].startswith("relative to"):
            continue
        printed = " | ".join(cells)
        if gas == "He" and cells[1]:
            for column, name in ((1, "chi_i_absolute"), (2, "chi_ii_absolute")):
                match = _FACTOR.match(cells[column])
                assert match is not None, cells[column]
                add(name, gas, Point("chi", RUBY, number(match["value"]), "", printed))
            for value, mark in _STACKED.findall(cells[3]):
                add(THEORY[mark], gas, Point("chi", RUBY, number(value), "", printed))
            continue
        if gas != "He" and cells[1]:
            for column, name in ((1, "chi_i"), (2, "chi_ii")):
                value, sigma, _ = _cell(cells[column])
                add(name, gas, Point("chi", RUBY, value, sigma, printed))
            if cells[3]:
                value, _, mark = _cell(cells[3])
                assert mark in {"b", "c"}, cells[3]  # Xe prints b: see the record
                add("dawes", gas, Point("chi", RUBY, value, "", printed))
        source = ""
        if cells[6]:
            value, sigma, mark = _cell(cells[6])
            add(
                KERR[mark],
                gas,
                Point("chi_K", HE_NE if mark == "e" else RUBY, value, sigma, printed),
            )
            source = DERIVED_ROW["K", mark]
        if cells[7]:
            value, sigma, mark = _cell(cells[7])
            add(STATIC[mark], gas, Point("chi_0", "inf", value, sigma, printed))
            source = DERIVED_ROW["0", mark]
        if cells[4]:
            value, _, _ = _cell(cells[4])
            add("derived", gas, Point(source, RUBY, value, "", printed))
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

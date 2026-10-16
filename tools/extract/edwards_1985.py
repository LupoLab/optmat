# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Edwards (1985): n and k of crystalline silicon, 2000-0.00372 eV (Table XII).

Table XII of the Handbook chapter runs over fourteen pages, each a library table (the
list ``PARTS``, in page order; pdf pp. 9-22, chapter pp. 555-568). Every row prints the
photon energy (eV), from 6.0 eV down the wavenumber (cm-1), the wavelength, and two
pairs of n and k columns. Each part was compared cell by cell with its table image; the
cells agree with the images except for the structural points below. The part on pdf
p. 15 has no CSV rendering (its "Silicon" title cell spans the columns) and every part
is read from its JSON cells, which are the same cells.

Structural points, handled here (each checked on the table images):

1. **Reference tags.** A tag in brackets after a value ("0.9999048 [30]", "1.010
   [20]") names the original of that value and of the values below it in the same
   column, until the next tag. The tags (text pp. 547-552): [30] Henke et al.
   (calculated by Edwards), [31] Gahwiller and Brown, [32] Brown et al., [29] Hunter
   (unpublished), [26] Philipp, [20] Aspnes and Theeten, [21] Hulthen, [23] McFarlane
   and Roberts, [1, 2, 3, 8] the selected infrared indices, [5, 6] Loewenstein et al.
   and Randall and Rawcliffe, [12] Hordvik and Skolnik, [15] Bendow et al., [14]
   Johnson, [19] Afsar and Hasted. One dataset per column kind (n or k) and tag; the
   two n and the two k columns share the mapping, since no tag appears in both.
2. **One untagged cell of another original.** At 630 A (19.68 eV) the k cell "9.3 x
   10^-2" is untagged under Philipp's "0.122 [26]", but it is Hunter's: its row's n is
   tagged [29], it lies on Hunter's 10-A grid between his 8.35e-2 (620 A) and 0.10
   (640 A), and Philipp's next value is the tagged 0.140 at 19.5 eV. ``TAG_OVERRIDES``
   gives it tag [29] without changing the column's running tag.
3. **Exponent inheritance.** A k cell printed as a bare mantissa (1.0 <= m < 10) after
   a cell "m x 10^-e" takes the same power of ten: "3.46" under "3.19 x 10^-5" is
   3.46e-5. The inheritance ends at the next cell with its own power, or at a plain
   decimal ("0.122"), and runs across empty cells, page breaks and tags, as the images
   print it: Hunter's first k, "1.0 [29]" at 125 A, is 1.0e-3 under the Henke column's
   1.54e-3; Brown et al.'s first, "1.08 [32]", is 1.08e-2 under Gahwiller and Brown's
   1.38e-2; Afsar and Hasted's first, "1.7 [19]", is 1.7e-4 under Johnson's 1.08e-4.
   The power is reset where the wavelength unit changes (the micrometre header of the
   p. 14 part): the second k column's next value, Aspnes and Theeten's "2.909 [20]" at
   6.0 eV, is a plain k of 2.909 (it continues smoothly to 0.945 at 3.26 eV), not a
   mantissa under Gahwiller and Brown's 4.9e-3 at 95 eV. Every inherited cell is
   written in E notation ("3.46e-5") with the raw line kept in ``printed``.
4. **Lines without an energy.** Six lines print neither energy nor wavelength, only a
   second value under the line above (Philipp's at 20.00 eV; Johnson's at 7.143,
   7.692, 8.333 and 9.091 um). Each takes the wavelength of the line above.
5. **Unit switch.** The wavelength column is in angstrom to 650 A and in micrometre
   from the p. 14 part (19.0 eV) on; every header row is read for its unit. Philipp's
   values span the switch (620 and 635.8 A, then 0.06526 um on), so they are two
   datasets per column kind, one per unit.

Printed misprints outside the values (the abscissa is the wavelength, right on each of
these rows, so none is an erratum): the 640 A row prints 10.37 eV for 19.37 eV; the
wavenumber is 1,492 for 1,429 cm-1 at 7.00 um, 608 for 606 at 16.50 um and 21,660 for
about 21,620 at 0.4626 um.

Row ids are the running index of the data lines of the whole table (1 = 2000 eV); the
abscissa is the wavelength as printed.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "EDWARDS_1985_SiliconSi"
SOURCE = "edwards_1985"

# Table XII, parts in page order (pdf pp. 9-22, chapter pp. 555-568).
PARTS = (
    "tbl_fe55a657b8ec",
    "tbl_646c7991193c",
    "tbl_c55280ac2760",
    "tbl_3ad03cca75e0",
    "tbl_5f02459e2946",
    "tbl_465e3c66d70e",
    "tbl_0f7731d732b9",
    "tbl_6a64c6e4156f",
    "tbl_9f49ec26fcd8",
    "tbl_66f11f14a901",
    "tbl_afc5f60ea2d9",
    "tbl_1eb73f7d4ec6",
    "tbl_f5e644f1b126",
    "tbl_f4f72cdf8301",
)

UNITS = {"Å": "angstrom", "μm": "um"}

# The value columns: (cell index, column kind). Columns 3-4 are the first n, k pair,
# 5-6 the second.
COLUMNS = ((3, "n"), (4, "k"), (5, "n"), (6, "k"))

# (column kind, reference tag, wavelength unit) -> dataset id. The unit only matters
# for Philipp's values, which are printed in both units.
DATASETS = {
    ("n", "30", "angstrom"): "henke_n",
    ("k", "30", "angstrom"): "henke_k",
    ("k", "31", "angstrom"): "gahwiller_brown_k",
    ("k", "32", "angstrom"): "brown_k",
    ("n", "29", "angstrom"): "hunter_n",
    ("k", "29", "angstrom"): "hunter_k",
    ("n", "26", "angstrom"): "philipp_n_angstrom",
    ("k", "26", "angstrom"): "philipp_k_angstrom",
    ("n", "26", "um"): "philipp_n_um",
    ("k", "26", "um"): "philipp_k_um",
    ("n", "20", "um"): "aspnes_theeten_n",
    ("k", "20", "um"): "aspnes_theeten_k",
    ("k", "21", "um"): "hulthen_k",
    ("k", "23", "um"): "mcfarlane_roberts_k",
    ("n", "1, 2, 3, 8", "um"): "ir_n",
    ("n", "5, 6", "um"): "far_ir_n",
    ("k", "12", "um"): "hordvik_skolnik_k",
    ("k", "15", "um"): "bendow_k",
    ("k", "14", "um"): "johnson_k",
    ("k", "19", "um"): "afsar_hasted_k",
}

# (cell index, printed wavelength) -> tag of an untagged cell that belongs to another
# original than the column's running tag (structural point 2).
TAG_OVERRIDES = {(4, "630"): "29"}

_TAG = re.compile(r"^(.*?)\s*\[([^\]]+)\]\s*$")
_POWER = re.compile(r"^(\d\.\d+)\s*\\times\s*10\^\{(-\d+)\}$")
_PLAIN = re.compile(r"^\d+\.?\d*$")


def _split_tag(cell: str) -> tuple[str, str | None]:
    """Return a cell's value text and its reference tag (None if untagged)."""
    text = cell.strip()
    match = _TAG.match(text)
    if match is None:
        return text, None
    return match.group(1).strip(), match.group(2).strip()


class _Column:
    """One value column of Table XII: its current reference tag and power of ten."""

    def __init__(self, index: int, kind: str) -> None:
        self.index = index
        self.kind = kind
        self.tag: str | None = None
        self.power: str | None = None

    def read(self, cell: str) -> tuple[str, str | None] | None:
        """Return the cell's value (E notation with a power of ten) and own tag."""
        text, tag = _split_tag(cell)
        if tag is not None:
            self.tag = tag
        if not text:
            return None
        match = _POWER.match(text)
        if match is not None:
            self.power = match.group(2)
            return f"{match.group(1)}e{self.power}", tag
        if not _PLAIN.match(text):
            msg = f"column {self.index}: unexpected cell {cell!r}"
            raise ValueError(msg)
        if text.startswith("0."):
            # A plain decimal ends the inheritance of a power of ten.
            self.power = None
            return text, tag
        if self.power is not None:
            # A bare mantissa under "m x 10^-e": the same power (see the docstring).
            if not 1.0 <= float(text) < 10.0:
                msg = f"column {self.index}: {cell!r} is not a mantissa"
                raise ValueError(msg)
            return f"{text}e{self.power}", tag
        return text, tag


def extract() -> dict[str, list[Point]]:
    """Return one data file per column kind, original and unit of Table XII."""
    columns = [_Column(index, kind) for index, kind in COLUMNS]
    files: dict[str, list[Point]] = {}
    unit = ""
    wavelength = ""
    index = 0
    for part in PARTS:
        for cells in json_table(DIRECTORY, part):
            first = cells[0].strip()
            if first == "Silicon":
                continue  # the title line of the p. 15 part
            if first == "eV":
                # A header row, at the top of each part.
                new = UNITS[cells[2].strip()]
                if unit and new != unit:
                    # The unit switch resets the powers (structural point 3).
                    for column in columns:
                        column.power = None
                unit = new
                continue
            index += 1
            if first or cells[2].strip():
                wavelength = cells[2].strip()
                float(wavelength)
            # Otherwise a line without an energy: the wavelength of the line above.
            printed = " | ".join(c.strip() for c in cells)
            for column in columns:
                read = column.read(cells[column.index])
                if read is None:
                    continue
                value, own = read
                tag = column.tag
                if own is None:
                    tag = TAG_OVERRIDES.get((column.index, wavelength), tag)
                if tag is None:
                    msg = f"row {index}: value before any reference tag"
                    raise ValueError(msg)
                dataset = DATASETS[column.kind, tag, unit]
                point = Point(str(index), wavelength, value, "", printed)
                files.setdefault(f"Si/{SOURCE}__{dataset}.csv", []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Philipp (1985): n, k and k_ABS of silica glass, 2000-0.00248 eV (Table X).

Table X of the Handbook chapter runs over eleven pages, each a library table (the list
``PARTS``, in page order). Every row prints the photon energy (eV), for the infrared
rows the wavenumber (cm-1), the wavelength, n and k, and from p. 12 on a column k_ABS
(k from transmission). Each part was compared cell by cell with its table image, twice;
the cells agree with the images except for the structural points below. Page numbers
"p. 5" to "p. 15" are pdf pages (chapter pp. 753-763).

Structural points, handled here (each checked on the table images):

1. **Reference tags.** A tag in brackets after a value ("0.99993 [17]", "0.325 [9, 10]")
   names the original of that value and of the values below it in the same column, until
   the next tag. The tags (pp. 750-752): [17] Henke et al. (computed by Philipp from
   atomic scattering factors), [16] Rife and Osantowski, [15] Hunter (unpublished),
   [9, 10] Philipp's smooth curve through Philipp's and Lamy's Kramers-Kronig results,
   [9] Philipp, [12] Weinberg et al., [12-14] a smooth curve through Weinberg's,
   Appleton's and Kaminow's transmission results, [7] Brixner's fit of Malitson's prism
   data, [2] Philipp 1979 (infrared), [1] Randall and Rawcliffe, and in the k_ABS column
   [4] Drummond and [2] Philipp 1979. One dataset per column and tag.
2. **Exponent inheritance.** A k or k_ABS cell printed as a bare mantissa (1.0 <= m <
   10) after a cell "m x 10^-e" takes the same power of ten: "1.636" under "1.503 x
   10^-5" is 1.636e-5. The inheritance ends at the next cell with its own power, or at a
   plain decimal ("0.107", "0.0339"), and runs across empty cells and page breaks (p. 7
   opens with "1.76" under p. 6's 10^-4, p. 8 with "6.19" under 10^-3, p. 9 with "8.2"
   under 10^-2). Every such cell is written in E notation ("1.636e-5") with the raw cell
   kept in ``printed``. The images show no power of ten that the extracted cells lose:
   each "x 10^-e" printed in the images is in the cells, and no bare mantissa in the
   images carries one. Two steps in k are absorption edges, not errors: at 529 eV (p. 6)
   the inherited "1.486" (1.486e-4) falls below 4.971e-4 at 542 eV, the oxygen K edge;
   at 1815 eV (p. 5) the printed 6.300e-6 falls below 1.936e-5 at 1860 eV, the silicon K
   edge.
3. **Unit switch.** Part ``tbl_19c745bbcfde`` (p. 9) repeats the header half way, with
   the wavelength column changing from angstrom to micrometre (from the 25 eV row).
   Every header row is read for its unit, and Hunter's values, which span the switch
   (300-490 A, then 0.0500-0.0600 um), are two datasets, one per unit. Part
   ``tbl_3ed3b519e90e`` (p. 12) likewise repeats the header to add the k_ABS column.

One printed misprint is an erratum in the record, not a fix here: the wavelength of the
22.96 eV row (p. 9) is printed "0.540" um for 0.0540.

Row ids are the running index of the data rows of the whole table (1 = 2000 eV); the
abscissa is the wavelength as printed.
"""

import re

from _library import Point, table, write

DIRECTORY = "PHILIPP_1985_SiliconDioxideSiO2"
SOURCE = "philipp_1985"

# Table X, parts in page order (pp. 5-15 of the chapter, pdf pages 5-15).
PARTS = (
    "tbl_252c849fca85",
    "tbl_163001deeaae",
    "tbl_dcfeca794199",
    "tbl_675d911da2e7",
    "tbl_19c745bbcfde",
    "tbl_f3221b4d7cc9",
    "tbl_7e7b88d780eb",
    "tbl_3ed3b519e90e",
    "tbl_b04abceaf3b9",
    "tbl_7604b3cb0747",
    "tbl_e7f8b49d5a9f",
)

UNITS = {"Å": "angstrom", "μm": "um"}

# (column, reference tag, wavelength unit) -> dataset id. The unit only matters for
# Hunter's values, which are printed in both units.
DATASETS = {
    ("n", "17", "angstrom"): "henke_n",
    ("k", "17", "angstrom"): "henke_k",
    ("n", "16", "angstrom"): "rife_osantowski_n",
    ("k", "16", "angstrom"): "rife_osantowski_k",
    ("n", "15", "angstrom"): "hunter_n_angstrom",
    ("k", "15", "angstrom"): "hunter_k_angstrom",
    ("n", "15", "um"): "hunter_n_um",
    ("k", "15", "um"): "hunter_k_um",
    ("n", "9, 10", "um"): "philipp_lamy_n",
    ("k", "9, 10", "um"): "philipp_lamy_k",
    ("n", "9", "um"): "philipp_vuv_n",
    ("k", "12", "um"): "weinberg_k",
    ("k", "12-14", "um"): "transmission_curve_k",
    ("n", "7", "um"): "brixner_malitson_n",
    ("n", "2", "um"): "philipp_ir_n",
    ("k", "2", "um"): "philipp_ir_k",
    ("k_abs", "4", "um"): "drummond_k_abs",
    ("k_abs", "2", "um"): "philipp_k_abs",
    ("n", "1", "um"): "randall_rawcliffe_n",
    ("k", "1", "um"): "randall_rawcliffe_k",
}

# Rife and Osantowski's sample is "a high-silica glass with 7.5% TiO2" (p. 751), not
# silica: its values are filed under that glass; every other dataset is SiO2.
MATERIALS = {"rife_osantowski_n": "SiO2-TiO2", "rife_osantowski_k": "SiO2-TiO2"}

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
    """One value column of Table X: its current reference tag and power of ten."""

    def __init__(self, name: str) -> None:
        self.name = name
        self.tag: str | None = None
        self.power: str | None = None

    def read(self, cell: str) -> str | None:
        """Return the cell's value (E notation when it has a power of ten), or None."""
        text, tag = _split_tag(cell)
        if tag is not None:
            self.tag = tag
        if not text:
            return None
        match = _POWER.match(text)
        if match is not None:
            self.power = match.group(2)
            return f"{match.group(1)}e{self.power}"
        if not _PLAIN.match(text):
            msg = f"column {self.name}: unexpected cell {cell!r}"
            raise ValueError(msg)
        if text.startswith("0."):
            # A plain decimal ends the inheritance of a power of ten.
            self.power = None
            return text
        if self.power is not None:
            # A bare mantissa under "m x 10^-e": the same power (see the docstring).
            if not 1.0 <= float(text) < 10.0:
                msg = f"column {self.name}: {cell!r} is not a mantissa"
                raise ValueError(msg)
            return f"{text}e{self.power}"
        return text


def extract() -> dict[str, list[Point]]:
    """Return one data file per value column and original of Table X."""
    columns = {name: _Column(name) for name in ("n", "k", "k_abs")}
    files: dict[str, list[Point]] = {}
    unit = ""
    index = 0
    for part in PARTS:
        for cells in table(DIRECTORY, part):
            if cells[0].strip() == "eV":
                # A header row, at the top of each part and where the layout changes.
                unit = UNITS[cells[2].strip()]
                continue
            index += 1
            wavelength = cells[2].strip()
            float(wavelength)
            printed = " | ".join(c.strip() for c in cells)
            values = {"n": cells[3], "k": cells[4], "k_abs": ""}
            if len(cells) > 5:
                values["k_abs"] = cells[5]
            for name, column in columns.items():
                value = column.read(values[name])
                if value is None:
                    continue
                if column.tag is None:
                    msg = f"row {index}: {name} value before any reference tag"
                    raise ValueError(msg)
                dataset = DATASETS[name, column.tag, unit]
                material = MATERIALS.get(dataset, "SiO2")
                point = Point(str(index), wavelength, value, "", printed)
                relative = f"{material}/{SOURCE}__{dataset}.csv"
                files.setdefault(relative, []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

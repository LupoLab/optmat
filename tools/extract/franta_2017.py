# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Franta et al. (2017): float-zone silicon's band gap, newAD2 table, n and k.

Table 2 (p. 409) lists the critical-point energies of the model; its first row is the
band gap E_g (the minimum energy of the indirect transitions) at 0 K and at 300 K, with
the average phonon temperature Theta of Eq. (54). Only E_g has a quantity; the other
energies are the record's formula coefficients. The table exists only as JSON cells
(no CSV); the cells lose the primes of some row labels (E'_0 against E_0), which does
not affect the E_g row (checked against the image).

The paper prints no n or k. Two tables of numbers go with it:

- The newAD2 site's unpublished universal dispersion model of float-zone c-Si,
  ``c-Si.Enk`` (attached to the paper in the library, with its URL): 4001 rows of photon
  energy (eV), n and k, 0.01-100 eV on a logarithmic grid, computed from ``c-Si.par``
  (the record's formula ``udm_c_si``). The rows up to 8.5 eV, the file's stated range,
  are extracted (check datasets); the rest is extrapolation.
- The refractiveindex.info pages ``main/Si/nk/Franta-<T>.yml`` ("Numerical data kindly
  provided by Daniel Franta"), kept in ``data/external/refractiveindex_info/``: the
  model of this paper as revised in Franta et al., J. Appl. Phys. 123, 185707 (2018).
  Each has one ``tabulated nk`` block of 4001 rows of wavelength (um), n and k,
  0.0309963-309.963 um. The pages at 20 and 25 degrees C and 300 K are extracted, at
  wavelengths up to 1e4/70 um (70 cm-1, the lowest wavenumber measured); the longer
  wavelengths are extrapolation. The script checks each page's REFERENCES and COMMENTS
  (the paper, the note, the 2018 update and the temperature) and keeps the numbers as
  the page prints them; the row id is the wavelength, or the energy for c-Si.Enk.
"""

from _json_table import json_table
from _library import DATA, Point, attachment, number, write

DIRECTORY = "Franta_2017_TemperatureDependentDispersion"
TABLE = "tbl_3d660fe93b4a"
PREFIX = "Si/franta_2017__"

ENK_ROWS = 4001
ENK_FITTED_EV = 8.5

DATABASE = DATA / "external" / "refractiveindex_info" / "main" / "Si" / "nk"
PAGE_ROWS = 4001
LONGEST_UM = 1e4 / 70  # 70 cm-1
# Page label -> the temperature its COMMENTS state.
PAGES = {
    "20C": "Temperature 20 °C.",
    "25C": "Temperature 25 °C.",
    "300K": "Temperature 300 K.",
}
PROVIDED = "(Numerical data kindly provided by Daniel Franta)"
PAPER = "Temperature-dependent dispersion model of float zone crystalline silicon"
UPDATE = "https://doi.org/10.1063/1.5026195"


def _band_gaps() -> dict[str, list[Point]]:
    """Return E_g at 0 K and at 300 K as two datasets."""
    rows = {cells[0]: cells for cells in json_table(DIRECTORY, TABLE)}
    cells = rows["E_{g}"]
    printed = " | ".join(cells)
    return {
        f"{PREFIX}band_gap_0K.csv": [Point("1", "", number(cells[1]), "", printed)],
        f"{PREFIX}band_gap_300K.csv": [Point("1", "", number(cells[2]), "", printed)],
    }


def _nk(rows: list[list[str]], n_file: str, k_file: str) -> dict[str, list[Point]]:
    """Return the n and k files of rows of abscissa, n and k as printed."""
    return {
        f"{PREFIX}{name}.csv": [
            Point(
                number(cells[0]),
                number(cells[0]),
                number(cells[column]),
                "",
                " ".join(cells),
            )
            for cells in rows
        ]
        for name, column in ((n_file, 1), (k_file, 2))
    }


def _newad2() -> dict[str, list[Point]]:
    """Return the newAD2 table's n and k up to 8.5 eV."""
    text = attachment(DIRECTORY, "c-Si.Enk").read_text(encoding="utf-8")
    rows = [line.split() for line in text.splitlines() if line.strip()]
    if len(rows) != ENK_ROWS or any(len(cells) != 3 for cells in rows):
        msg = f"c-Si.Enk: {len(rows)} rows, not {ENK_ROWS} of three numbers"
        raise ValueError(msg)
    fitted = [cells for cells in rows if float(cells[0]) <= ENK_FITTED_EV]
    return _nk(fitted, "newad2_n", "newad2_k")


def _page(label: str) -> list[list[str]]:
    """Return a page's rows of wavelength, n and k, after checking its preamble."""
    name = f"Franta-{label}.yml"
    lines = (DATABASE / name).read_text(encoding="utf-8").splitlines()
    text = "\n".join(lines)
    if (
        any(part not in text for part in (PROVIDED, PAPER, UPDATE, PAGES[label]))
        or text.count("type:") != 1
    ):
        msg = f"{name}: not the expected single Franta nk block"
        raise ValueError(msg)
    start = lines.index("  - type: tabulated nk") + 1
    if lines[start].strip() != "data: |":
        msg = f"{name}: unexpected block {lines[start - 1 : start + 1]}"
        raise ValueError(msg)
    rows = [line.split() for line in lines[start + 1 :] if line.strip()]
    if len(rows) != PAGE_ROWS or any(len(cells) != 3 for cells in rows):
        msg = f"{name}: {len(rows)} rows, not {PAGE_ROWS} of three numbers"
        raise ValueError(msg)
    return [cells for cells in rows if float(cells[0]) <= LONGEST_UM]


def extract() -> dict[str, list[Point]]:
    """Return the band gaps, the newAD2 table's n and k, and the pages' n and k."""
    files = _band_gaps() | _newad2()
    for label in PAGES:
        files |= _nk(_page(label), f"rii_{label}_n", f"rii_{label}_k")
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

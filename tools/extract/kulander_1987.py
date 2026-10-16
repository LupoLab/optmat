# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kulander (1987): time-dependent Hartree-Fock ionisation of helium.

The library copy is a scan of Phys. Rev. A 36, 2726 (1987): every page one 1-bit image
(5169 x 7129 pixels) with an OCR text layer; PDF page n is journal page 2725 + n. The
only numbers are those of **Table I** (p. 2732, library table ``tbl_8ef576c11381``,
JSON only; every cell checked against the table image):

* **(a)** 'Helium ionization rates (s^-1)' against the peak intensity (W/cm^2), at
  1.064 and 0.532 um (5e14, 1e15, 2e15 W/cm^2) and at a photon energy of 9.0 eV (1e15
  W/cm^2), each in two columns headed '1s^2' and "1s1s'": the time-dependent
  Hartree-Fock calculation (one doubly occupied orbital), its rate extrapolated to
  the static Hartree-Fock orbital energy, and the frozen one-electron model (one
  electron frozen in its ground-state orbital), as the record reads the headings.
  One file per wavelength and column.
* **(b)** '(generalized) cross sections': the two-photon sigma in cm^4/W^2 s at
  photon energies 13.6 and 17.7 eV and the one-photon sigma in cm^2 at 30.0 eV, in the
  columns TDHF (this paper's) and PT ('Refs. 25-27', 'Refs. 28': other papers'
  perturbation theory, not taken). One file for the two-photon TDHF values, one for
  the one-photon one.

Every cell is printed with its own power of ten ('4.7\times 10^{11}'): the value (and
the intensity) is written as the mantissa times that power, ``4.7e11``, and
``printed`` keeps the row's cells. The two-photon unit, cm^4 W^-2 s^-1, is held as
printed (``s-1 cm4/W2``): W / I^2, W the rate and I the intensity; the script also
reports sigma_2 = (W / I^2) (hbar omega)^2 in cm^4 s for comparison with the other
records (:func:`checks`), which the data files do not hold.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Kulander_1987_TimeDependentHartree"
TABLE_I = "tbl_8ef576c11381"
HEADER = ["I (W/cm2) \\ wavelength (μm)", "(a)", "", "", ""]
WAVELENGTHS = ["", "1.064", "", "0.532", ""]
COLUMNS = ["", "1s^2", "1s1s'", "1s^2", "1s1s'"]
LABELS = {"1s^2": "1s2", "1s1s'": "1s1sp"}
"""The column headings and the file labels: TDHF (1s^2) and frozen (1s1s')."""
_CELL = re.compile(r"^\s*([0-9.]+)\s*\\times\s*10\^\{(-?\d+)\}\s*$")
EV_J = 1.602176634e-19
"""The electronvolt, J (exact, SI 2019)."""


def _number(cell: str) -> str:
    """Return a cell's number as mantissa e power, as printed ('4.7e11')."""
    match = _CELL.match(cell)
    if match is None:
        msg = f"Table I: not a number times a power of ten: {cell!r}"
        raise ValueError(msg)
    return f"{match[1]}e{match[2]}"


def table() -> list[list[str]]:
    """Return Table I's rows, raising if its layout moved."""
    rows = json_table(DIRECTORY, TABLE_I)
    expected = {
        0: HEADER,
        1: WAVELENGTHS,
        2: COLUMNS,
        6: ["", "9.0 eV", "", "", ""],
        8: ["", "(b)", "", "", ""],
        9: ["", "TDHF", "PT", "", "\\hbar\\omega (eV)"],
    }
    if len(rows) != 13 or any(rows[i] != cells for i, cells in expected.items()):
        msg = "Table I: the layout moved"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return Table I(a)'s rates and Table I(b)'s TDHF cross sections."""
    rows = table()
    files: dict[str, list[Point]] = {}
    for start, column in ((1, "1064"), (3, "532")):
        for offset in (0, 1):
            label = LABELS[rows[2][start + offset]]
            files[f"He/kulander_1987__tbl1a_rate_{column}_{label}.csv"] = [
                Point(
                    _number(cells[0]),
                    _number(cells[0]),
                    _number(cells[start + offset]),
                    "",
                    " | ".join(cells),
                )
                for cells in rows[3:6]
            ]
    for offset in (0, 1):
        label = LABELS[rows[2][1 + offset]]
        cells = rows[7]
        files[f"He/kulander_1987__tbl1a_rate_9eV_{label}.csv"] = [
            Point(
                _number(cells[0]),
                _number(cells[0]),
                _number(cells[1 + offset]),
                "",
                " | ".join(cells),
            )
        ]
    two_photon = []
    for cells in rows[10:12]:
        if cells[0] != "\\sigma (cm4/W2s)" or cells[3] != "(Refs. 25-27)":
            msg = f"Table I(b): an unexpected row {cells}"
            raise ValueError(msg)
        two_photon.append(
            Point(cells[4], cells[4], _number(cells[1]), "", " | ".join(cells))
        )
    files["He/kulander_1987__tbl1b_sigma2_tdhf.csv"] = two_photon
    cells = rows[12]
    if cells[0] != "\\sigma (cm2)" or cells[3] != "(Refs. 28)":
        msg = f"Table I(b): an unexpected row {cells}"
        raise ValueError(msg)
    files["He/kulander_1987__tbl1b_sigma1_tdhf.csv"] = [
        Point(cells[4], cells[4], _number(cells[1]), "", " | ".join(cells))
    ]
    return files


def checks() -> dict[str, str]:
    """Return Table I(b)'s two-photon values as sigma_2 = (W / I^2) (hbar omega)^2.

    Both columns, TDHF and PT, in cm^4 s, for W = sigma_2 Phi^2 with Phi = I / (hbar
    omega) the photon flux; reported, not written.
    """
    out = {}
    for cells in table()[10:12]:
        photon = float(cells[4]) * EV_J
        for name, cell in (("TDHF", cells[1]), ("PT", cells[2])):
            sigma = float(_number(cell)) * photon**2
            out[f"{cells[4]} eV, {name}: sigma_2 (cm4 s)"] = f"{sigma:.3e}"
    return out


def main() -> None:
    """Write the data files and print the conversions."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

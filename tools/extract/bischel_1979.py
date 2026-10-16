# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bischel, Bokor, Kligler & Rhodes (1979): Xe two-photon ionisation estimate.

Table I (p. 381, library table ``tbl_246cdac4e41d``, JSON cells only), 'Estimates of
two-photon excitation cross sections', prints the coefficient alpha in cm^4/watt,
defined by sigma(nu_2) = alpha I_1 (Eq. 1) and the transition rate W_fg = alpha I^2 /
hbar omega (p. 382), for seven systems at 193 and 248 nm. Only the xenon row at 193
nm, whose excited product is Xe+ (two-photon ionisation into the autoionising region,
p. 386), is a quantity of the vocabulary: one file. The other rows are bound
two-photon excitations (no quantity holds a linewidth-dependent alpha; the record's
notes). The row's system cell is empty under the merged 'Xe' cell of the 248 nm row;
the script requires both. The cells were compared with the table image (2026-10-09):
no misread cell.

**Checks** (:func:`checks`): the row's alpha must be the value the text prints for
the same process, 'sigma / I ~ 10^-34 cm^4/W (see Table I)' (p. 388), and the value
mccown_1982 reproduces as Bischel et al.'s (its ``Xe/mccown_1982__bischel_Xe.csv``);
the bound rows' linewidth, footnote (a) '3 x 10^6 MHz', must be the text's '100
cm^-1' (p. 381) within 0.1 %, while the Xe+ row's is footnote (g)'s '~10^9 MHz',
'estimated for photoionization continuum'. Reported: sigma_2 = alpha hbar omega at
193 nm, and the ion density the text estimates from alpha (p. 388).
"""

import csv
import math
import re

from _json_table import json_table
from _library import DATA, Point, number, write

from optmat.constants import HBAR, C

DIRECTORY = "Bischel_1979_NonlinearOpticalProcesses"
TABLE_I = "tbl_246cdac4e41d"
FILE = "Xe/bischel_1979__alpha_Xe.csv"
ALPHA_HEADER = r"\alpha (cm^{4}/watt)"
WIDTH_HEADER = "Line Width(MHz)"
POWER = -34
"""The power of ten of the Xe+ cell; the record's scale."""
_POWER = re.compile(r"\\times\s*10\^\{(-?\d+)\}")

TEXT_ALPHA = 1e-34
"""p. 388: 'The characteristic cross sections and rates are sigma / I ~ 10^-34 cm^4/W
(see Table I)', for reaction (18), the two-photon ionisation of Xe at 193 nm."""
LASER_WIDTH_CM1 = 100.0
"""p. 381: Table I's estimates 'assuming a laser linewidth of 100 cm^-1'."""
MHZ_PER_CM1 = C * 100 / 1e6
"""One cm^-1 of linewidth in MHz."""
WAVELENGTH_M = 193e-9
# p. 388: 'Typical experimental conditions with an Xe density ~3 x 10^17 cm^-3 and an
# optical intensity of 10^9 W/cm^2 result in an estimated ion density of ~3 x 10^11
# cm^-3', with 'the laser pulse (20 ns)'.
DENSITY_CM3 = 3e17
INTENSITY_W_CM2 = 1e9
PULSE_S = 20e-9
PRINTED_IONS_CM3 = 3e11


def _rows() -> tuple[list[str], list[str], list[str]]:
    """Return Table I's header, its 248 nm xenon row and its 193 nm Xe+ row."""
    grid = json_table(DIRECTORY, TABLE_I)
    header = grid[0]
    if header[7] != ALPHA_HEADER or header[4] != WIDTH_HEADER:
        msg = f"Table I's header has changed: {header}"
        raise ValueError(msg)
    xenon = next(k for k, row in enumerate(grid) if row[0] == "Xe")
    upper, lower = grid[xenon], grid[xenon + 1]
    if upper[1] != "248" or lower[:3] != ["", "193", "Xe^{+}"]:
        msg = f"Table I's xenon rows have changed: {upper[:3]}, {lower[:3]}"
        raise ValueError(msg)
    return header, upper, lower


def extract() -> dict[str, list[Point]]:
    """Return the Xe+ row's alpha at 193 nm, the power of ten checked."""
    _, _, row = _rows()
    cell = row[7]
    match = _POWER.search(cell)
    if match is None or int(match.group(1)) != POWER:
        msg = f"Table I, Xe+ alpha: unexpected cell {cell!r}"
        raise ValueError(msg)
    printed = " | ".join(("Xe", *row[1:8]))
    return {FILE: [Point("Xe+ 193", "193", number(cell), "", printed)]}


def _mccown() -> float:
    """Return the value mccown_1982 reproduces as Bischel et al.'s, in cm^4 W^-1."""
    path = DATA / "values" / "ionisation" / "Xe" / "mccown_1982__bischel_Xe.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        (row,) = list(csv.DictReader(handle))
    return float(row["value"]) * 1e-34


def checks() -> dict[str, str]:
    """Raise unless the table, the text and mccown_1982's copy agree; report more."""
    _, bound, ionisation = _rows()
    (point,) = extract()[FILE]
    alpha = float(point.value) * 10.0**POWER
    failures = []
    if alpha != TEXT_ALPHA:
        failures.append(f"Table I's {alpha:.1e} is not the text's {TEXT_ALPHA:.0e}")
    if alpha != _mccown():
        failures.append(f"mccown_1982 reproduces {_mccown():.1e}, not {alpha:.1e}")
    width = float(number(bound[4].removesuffix("(a)"))) * 10.0**6
    if abs(width / MHZ_PER_CM1 / LASER_WIDTH_CM1 - 1) > 1e-3:
        failures.append(f"footnote (a)'s {width:.3g} MHz is not 100 cm-1")
    if ionisation[4] != r"\sim 10^{9}(g)":
        failures.append(f"the Xe+ row's linewidth is {ionisation[4]!r}")
    if failures:
        raise ValueError("; ".join(failures))
    photon_j = HBAR * 2 * math.pi * C / WAVELENGTH_M
    ions = alpha * INTENSITY_W_CM2**2 / photon_j * PULSE_S * DENSITY_CM3
    return {
        "alpha (Table I, p. 388 text, mccown_1982)": f"{alpha:.0e} cm4 W-1",
        "footnote (a) width": f"{width / MHZ_PER_CM1:.2f} cm-1 (text: 100)",
        "sigma_2 = alpha hbar omega at 193 nm": f"{alpha * photon_j:.3e} cm4 s",
        "ion density, p. 388 conditions": (
            f"{ions:.2e} cm-3 (printed ~{PRINTED_IONS_CM3:.0e}; "
            f"ratio {ions / PRINTED_IONS_CM3:.2f})"
        ),
    }


def main() -> None:
    """Write the data file and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

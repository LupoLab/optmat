# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Kroll and Bischel (1990): two-photon and excited-state cross sections of Xe.

The PDF is a scan (Lund University's repository copy, a cover and a licence page before
the article: journal p. 1340 is PDF p. 3). Both tables have merged header cells and no
CSV rendering, so they are read from their JSON cells; every cell was checked against
the table images and the rendered pages.

* **Table I** (``tbl_3fbaf4afd365``, p. 1346): the two-photon cross sections to the
  2p5 (J = 0), 2p6 (J = 2) and 2p9 (J = 2) levels of the 5p5 6p configuration (Paschen
  notation), in 10^-35 cm^4. Its upper part is the measured product sigma_tot G2(0)
  ('This work', row 2), with other authors' values (refs. 10, 8, 9; not extracted);
  its lower part sigma_tot, 'assuming G2(0) = 1.4' ('This work', row 7), with a
  theoretical range (ref. 12; not extracted). The uncertainty column: '30%' with note b,
  '50% for the 2p9 state'.
* **Table II** (``tbl_0a9507c3ebd9``, p. 1347): the photoionisation cross sections of
  the 2p5 and 2p6 levels at 248 and 193 nm, 10^-18 cm^2, ``value+-sigma``. The theory
  rows (refs. 13, 14) are printed centred across both states (one value per wavelength,
  a range for ref. 13) and are not extracted.

:func:`checks` reproduces what the paper derives from Table I: sigma_tot = sigma_tot
G2(0) / 1.4 to the printed figures, and Sec. V's resonant cross section 4.0e-45 cm^4 s
as the 2p5 sigma_tot G2(0) times the peak of a Gaussian of the Doppler width it states
(0.086 cm^-1), normalised in the two-photon angular frequency.
"""

import math

from _json_table import json_table
from _library import Point, number, write

from optmat.constants import ATOMIC_MASS, K_B, C
from optmat.gas.eos import equation_of_state

DIRECTORY = "Kroll_1990_TwoPhotonAbsorption"
TABLE_I = "tbl_3fbaf4afd365"
TABLE_II = "tbl_0a9507c3ebd9"
STATES = ("2p5", "2p6", "2p9")
"""Table I's columns 1-3, in order (headers ``2p_{5} state(J=0)`` ...)."""
MEASURED_ROW, CORRECTED_ROW = 2, 7
"""Table I's 'This work' rows under sigma_tot G2(0) and under sigma_tot."""
G2 = 1.4
"""G2(0), 'assuming this value holds also for the first anti-Stokes light' (p. 1346)."""
WAVELENGTHS = {
    2: ("2p5", "248"),
    3: ("2p6", "248"),
    4: ("2p5", "193"),
    5: ("2p6", "193"),
}
"""Table II's columns: state and ionising wavelength (nm) of each."""

DOPPLER_CM = 0.086
"""Sec. V (p. 1347): 'Doppler width 0.086 cm^-1' (FWHM, of the two-photon line)."""
SECTION_V_SIGMA = 4.0e-45
"""Sec. V: 'two-photon absorption cross section sigma^(2) = 4.0 x 10^-45 cm4 s'."""
SECTION_V_TOLERANCE = 0.03
"""The two-figure rounding of 7.0, 0.086 and 4.0 (0.7 %, 0.6 % and 1.3 %): the corrected
sigma_tot (5.0) would miss by 28 %."""
LEVEL_2P5_CM = 80119.0
"""The 2p5 level, '80 119 cm-1' (p. 1340), for the Doppler width check."""
ROOM_K = 295.0
"""Room temperature for the Doppler width check (the paper says 'room temperature')."""


def tables() -> tuple[list[list[str]], list[list[str]]]:
    """Return Tables I and II as cell grids."""
    return json_table(DIRECTORY, TABLE_I), json_table(DIRECTORY, TABLE_II)


def checks() -> dict[str, object]:
    """Return the checks of Table I's two parts and of Sec. V's resonant value.

    Raises if one fails.
    """
    out: dict[str, object] = {}
    failures = []
    grid, _ = tables()
    measured = grid[MEASURED_ROW]
    corrected = grid[CORRECTED_ROW]
    if measured[0] != "This work" or corrected[0] != "This work":
        failures.append("Table I's 'This work' rows moved")
    for column, state in enumerate(STATES, start=1):
        quotient = float(measured[column]) / G2
        out[f"{state}: sigma_tot G2(0) / 1.4"] = (quotient, corrected[column])
        if f"{quotient:.1f}" != corrected[column]:
            failures.append(
                f"{state}: {measured[column]} / 1.4 is not {corrected[column]}"
            )
    # Sec. V: sigma^(2)(0) = sigma_tot G2(0) x the Gaussian peak 2 sqrt(ln 2 / pi) /
    # Delta omega, Delta omega = 2 pi c x 0.086 cm^-1 (two-photon angular frequency).
    peak = 2 * math.sqrt(math.log(2) / math.pi) / (2 * math.pi * C * DOPPLER_CM * 100)
    product = float(measured[1]) * 1e-35 * peak
    out["Sec. V: 2p5 sigma_tot G2(0) x Doppler peak (cm4 s)"] = (
        product,
        SECTION_V_SIGMA,
    )
    if abs(product / SECTION_V_SIGMA - 1) > SECTION_V_TOLERANCE:
        failures.append(f"Sec. V's 4.0e-45 is not {product:.2e}")
    out["Sec. V: 2p5 sigma_tot x Doppler peak (cm4 s), not printed"] = (
        float(corrected[1]) * 1e-35 * peak
    )
    mass = equation_of_state("Xe").molar_mass / 1e-3 * ATOMIC_MASS
    width = LEVEL_2P5_CM * math.sqrt(8 * math.log(2) * K_B * ROOM_K / (mass * C**2))
    out["Doppler FWHM of the 2p5 line at 295 K (cm-1)"] = (width, DOPPLER_CM)
    if abs(width - DOPPLER_CM) > 0.0015:
        failures.append(f"Doppler width {width:.4f} cm-1, printed 0.086")
    ratio = float(measured[1]) / float(measured[2])
    out["2p5 / 2p6 ('~3', p. 1347)"] = ratio
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's two 'This work' rows and Table II's measured cross sections."""
    checks()
    table_i, table_ii = tables()
    files: dict[str, list[Point]] = {}
    for row, name in ((MEASURED_ROW, "sigma_tot_g2"), (CORRECTED_ROW, "sigma_tot")):
        cells = table_i[row]
        files[f"Xe/kroll_1990__{name}.csv"] = [
            Point(state, state, number(cells[column]), "", " | ".join(cells))
            for column, state in enumerate(STATES, start=1)
        ]
    cells = table_ii[2]
    if cells[0] != "Experiment":
        msg = "Table II's experiment row moved"
        raise ValueError(msg)
    for column, (state, wavelength) in WAVELENGTHS.items():
        value, _, sigma = cells[column].partition("±")
        files.setdefault(f"Xe/kroll_1990__photoionization_{state}.csv", []).append(
            Point(
                wavelength,
                wavelength,
                number(value),
                number(sigma),
                " | ".join(cells),
            )
        )
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

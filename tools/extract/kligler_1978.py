# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kligler, Pritchard, Bischel & Rhodes (1978): Xe two-photon ionisation.

Table I (p. 2222, library table ``tbl_51477d804fd7``, JSON cells only), 'Estimates of
two photon coupling of uv radiation to atomic and molecular systems', prints sigma/I in
cm^4/W, the cross section of Eq. (1), R = sigma (I / hbar omega) n, per unit intensity
(so the rate per atom is (sigma/I) I^2 / hbar omega), estimated from Eqs. (2)-(3) with
one intermediate state, mu_1^2 = mu_2^2 = 1 D^2 and an estimated linewidth. Only the
xenon row at 1930 Angstrom, whose excited product is Xe+ (two-photon ionisation), is a
quantity of the vocabulary: one file. The other rows are bound two-photon excitations
at line centre for an assumed linewidth (no quantity holds them; the record's notes).
The row's system cell is empty under the 'Xe' cell of the 2484 Angstrom row; the script
requires both. The cells were compared with the table image (2026-10-09): no misread
cell; the extraction drops the superscripts ('1 [times] 10-34' is 1 x 10^-34).

**Checks** (:func:`checks`): the 2484 Angstrom xenon row's sigma/I times the text's
intensity, 3 x 10^9 W/cm^2, must round to the text's 'sigma ~ 10^-21 cm2' (p. 2222);
the Xe+ value must equal the one bischel_1979 prints for the same estimate (its
``Xe/bischel_1979__alpha_Xe.csv``); the lower bound of p. 2220 must follow from the
printed counts, volume, duration, intensity and density (8 x 10^8 photons in 2 x 10^-4
cm^3 over 8 ns is 5 x 10^20 cm^-3 s^-1; sigma = R / ((I / hbar omega) n) >= 5 x 10^-27
cm^2). Reported: Eqs. (2)-(3) evaluated with Table I's inputs for each row and four
readings of the unstated linewidth factor g (Lorentzian or Gaussian, normalised in the
ordinary or the angular frequency), and the Xe 6p estimate in the per-level convention
of TC-952.
"""

import csv
import math
import re

from _json_table import json_table
from _library import DATA, Point, write

from optmat.constants import AU_ANGULAR_FREQUENCY, HBAR, C
from optmat.gas.two_photon_bound import SIGMA0_PER_DELTA

DIRECTORY = "Kligler_1978_XeOProductionThrough"
TABLE_I = "tbl_51477d804fd7"
FILE = "Xe/kligler_1978__sigma_over_I_Xe.csv"
HEADER = [
    "System",
    "Wave-length λ (Å)",
    "Excited product",
    "Energy denominator Δω (cm-1)",
    "Linewidth (MHz)",
    r"\mu_1^2(D)2",
    r"\mu_2^2(D)2",
    "\u03c3/I (cm4/W)",
    "Remarks",
]
POWER = -34
"""The power of ten of the Xe+ cell; the record's scale."""
_CELL = re.compile("^~?(\\d+(?:\\.\\d+)?)\u00d710(-\\d+)$")
"""A sigma/I cell as the extraction gives it, '3 [times] 10-31' for 3 x 10^-31."""
_DENOMINATOR = re.compile("^(\\d+)\u00d7103")
"""An energy denominator cell, '28 [times] 103' for 28 x 10^3 cm^-1."""

TEXT_INTENSITY_W_CM2 = 3e9
"""p. 2222: 'for a laser intensity of 3 x 10^9 W/cm^2'."""
TEXT_SIGMA_CM2 = 1e-21
"""p. 2222: 'the two-photon absorption in Xe will have sigma ~ 10^-21 cm^2'."""
# p. 2220: 'about 1.6 x 10^6 counts', 'overall detection efficiency of about 0.2 %',
# '8 x 10^8 photons', 'duration of peak laser power to be 8 nsec', 'focal volume of 2
# x 10^-4 cm^3', 'R_Xe = 5 x 10^20 cm^-3 sec^-1', '2.5 x 10^19 cm^-3 Xe', 'sigma_Xe >=
# 5 x 10^-27 cm^2 at our intensity', the intensity '3 x 10^9 W/cm^2' (p. 2219).
COUNTS = 1.6e6
EFFICIENCY = 0.002
PHOTONS = 8e8
DURATION_S = 8e-9
VOLUME_CM3 = 2e-4
RATE_CM3_S = 5e20
XE_DENSITY_CM3 = 2.5e19
SIGMA_BOUND_CM2 = 5e-27
KRF_M = 248.4e-9
"""The KrF* wavelength the paper states (248.4 nm, p. 2219 and Fig. 3)."""
DEBYE_ESU_CM = 1e-18
DEBYE_AU = 0.393430
"""One debye in atomic units (e a0), to six figures."""
CM1_PER_HARTREE = 219474.6313632
EQUATION_2_ROWS = frozenset({"H_2", "N_2", "Kr", "I_2", "Xe", ""})
"""The rows with one intermediate state and plain moments ('' is the Xe+ row); CO's
two intermediate states and N2O's approximate moments are not recomputed."""


def _grid() -> list[list[str]]:
    """Return Table I's rows, the header checked."""
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0] != HEADER:
        msg = f"Table I's header has changed: {grid[0]}"
        raise ValueError(msg)
    return grid


def _xenon() -> tuple[list[str], list[str]]:
    """Return Table I's 2484 Angstrom xenon row and its 1930 Angstrom Xe+ row."""
    grid = _grid()
    xenon = next(k for k, row in enumerate(grid) if row[0] == "Xe")
    upper, lower = grid[xenon], grid[xenon + 1]
    if upper[1:3] != ["2484", "Xe*(5p56p)"] or lower[:3] != ["", "1930", "Xe+"]:
        msg = f"Table I's xenon rows have changed: {upper[:3]}, {lower[:3]}"
        raise ValueError(msg)
    return upper, lower


def _sigma_over_i(cell: str) -> tuple[str, int]:
    """Return a sigma/I cell's mantissa as printed and its power of ten."""
    match = _CELL.match(cell)
    if match is None:
        msg = f"Table I: unexpected sigma/I cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), int(match.group(2))


def extract() -> dict[str, list[Point]]:
    """Return the Xe+ row's sigma/I at 1930 Angstrom, the power of ten checked."""
    _, row = _xenon()
    mantissa, power = _sigma_over_i(row[7])
    if power != POWER:
        msg = f"Table I, Xe+ sigma/I: power {power}, expected {POWER}"
        raise ValueError(msg)
    printed = " | ".join(("Xe", *row[1:8]))
    return {FILE: [Point("Xe+ 1930", "1930", mantissa, "", printed)]}


def _bischel() -> float:
    """Return bischel_1979's alpha for the same Xe+ estimate, in cm^4 W^-1."""
    path = DATA / "values" / "ionisation" / "Xe" / "bischel_1979__alpha_Xe.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        (row,) = list(csv.DictReader(handle))
    return float(row["value"]) * 1e-34


def _photon_j(wavelength_m: float) -> float:
    """Return the photon energy at a wavelength, J."""
    return 2 * math.pi * HBAR * C / wavelength_m


def _equation_2(row: list[str], g_s: float) -> float:
    """Return Eqs. (2)-(3) per unit intensity with a row's inputs, cm^4/W (CGS).

    sigma / I = (2 pi)^2 / (hbar c^2) omega |M|^2 g, M = 2 mu_1 mu_2 / (hc Delta), one
    intermediate state; g in s (its normalisation is the caller's reading).
    """
    angstrom = float(row[1].removesuffix("e"))
    match = _DENOMINATOR.match(row[3])
    if match is None:
        msg = f"Table I: unexpected energy denominator {row[3]!r}"
        raise ValueError(msg)
    delta_cm1 = float(match.group(1)) * 1e3
    mu1, mu2 = (float(cell.removeprefix("~").removesuffix("^c")) for cell in row[5:7])
    hbar = HBAR * 1e7
    c = C * 100
    omega = 2 * math.pi * c / (angstrom * 1e-8)
    energy = 2 * math.pi * hbar * c * delta_cm1
    m = 2 * math.sqrt(mu1 * mu2) * DEBYE_ESU_CM**2 / energy
    return (2 * math.pi) ** 2 / (hbar * c**2) * omega * m**2 * g_s * 1e7


def _widths(row: list[str]) -> float:
    """Return a row's linewidth in Hz from its cell ('10^6b', '~109g')."""
    cell = row[4].removeprefix("~")
    exponent = {"10^6b": 6, "109g": 9}[cell]
    return 10.0**exponent * 1e6


def checks() -> dict[str, str]:
    """Raise unless the table, the text and bischel_1979 agree; report the rest."""
    bound, _ = _xenon()
    mantissa, power = _sigma_over_i(bound[7])
    sigma_over_i = float(mantissa) * 10.0**power
    (point,) = extract()[FILE]
    alpha = float(point.value) * 10.0**POWER
    failures = []
    sigma = sigma_over_i * TEXT_INTENSITY_W_CM2
    if round(math.log10(sigma)) != round(math.log10(TEXT_SIGMA_CM2)):
        failures.append(f"Table I x 3e9 W/cm2 is {sigma:.1e}, the text ~1e-21")
    if alpha != _bischel():
        failures.append(f"bischel_1979 prints {_bischel():.1e}, not {alpha:.1e}")
    photons = COUNTS / EFFICIENCY
    rate = PHOTONS / (DURATION_S * VOLUME_CM3)
    flux = TEXT_INTENSITY_W_CM2 / _photon_j(KRF_M)
    bound_sigma = RATE_CM3_S / (flux * XE_DENSITY_CM3)
    if abs(photons / PHOTONS - 1) > 1e-9 or abs(rate / RATE_CM3_S - 1) > 1e-9:
        failures.append(f"photons {photons:.2e}, rate {rate:.2e}")
    if not SIGMA_BOUND_CM2 <= bound_sigma < 1.1 * SIGMA_BOUND_CM2:
        failures.append(f"the lower bound recomputes to {bound_sigma:.2e} cm2")
    if failures:
        raise ValueError("; ".join(failures))
    report = {
        "Xe 6p sigma/I x 3e9 W/cm2 (text: ~1e-21 cm2)": f"{sigma:.1e} cm2",
        "Xe+ sigma/I (bischel_1979: 1e-34)": f"{alpha:.0e} cm4 W-1",
        "Xe+ sigma_2 = (sigma/I) hbar omega at 193.0 nm": (
            f"{alpha * _photon_j(193.0e-9):.3e} cm4 s"
        ),
        "lower bound, p. 2220": (
            f"R = {rate:.1e} cm-3 s-1, sigma >= {bound_sigma:.2e} cm2, sigma/I >= "
            f"{bound_sigma / TEXT_INTENSITY_W_CM2:.2e} cm4/W, sigma_2 >= "
            f"{bound_sigma / flux:.2e} cm4 s"
        ),
    }
    root = math.sqrt(math.log(2) / math.pi)
    for row in _grid()[1:]:
        if row[0] not in EQUATION_2_ROWS:
            continue
        width = _widths(row)
        readings = {
            "Lorentz/Hz": 2 / (math.pi * width),
            "Gauss/Hz": 2 * root / width,
            "Lorentz/(rad/s)": 2 / (math.pi * 2 * math.pi * width),
            "Gauss/(rad/s)": 2 * root / (2 * math.pi * width),
        }
        name = row[0] or "Xe+"
        values = ", ".join(
            f"{k} {_equation_2(row, g):.2g}" for k, g in readings.items()
        )
        report[f"Eq. (2), {name} {row[1]} (printed {row[7]})"] = values
    mu2 = DEBYE_AU**2
    delta = (2 * mu2 / (28e3 / CM1_PER_HARTREE)) ** 2
    omega_res = 2 * math.pi * C * 80118.962e2 / 2
    sigma0_inputs = SIGMA0_PER_DELTA * (omega_res / AU_ANGULAR_FREQUENCY) ** 2 * delta
    lorentz = 2 / (math.pi * 2 * math.pi * 1e12)
    sigma0_printed = sigma_over_i * _photon_j(KRF_M) / lorentz
    delta_printed = (
        sigma0_printed
        * 1e-8
        / (SIGMA0_PER_DELTA * (omega_res / AU_ANGULAR_FREQUENCY) ** 2)
    )
    report["Xe 6p, TC-952 convention from Table I's inputs"] = (
        f"delta = |M|^2 = {delta:.2f} au, sigma_0 = {sigma0_inputs * 1e8:.2e} cm4"
    )
    report["Xe 6p, TC-952 convention from the printed sigma/I, 1 THz Lorentzian"] = (
        f"sigma_0 = {sigma0_printed:.2e} cm4, delta = {delta_printed:.1f} au"
    )
    return report


def main() -> None:
    """Write the data file and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

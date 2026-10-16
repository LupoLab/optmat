# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Check Sharma et al. (2018)'s record against the paper: Table 1 and its constants.

The source's one value, the O2 eight-photon cross section sigma_8 = (3.3 +- 0.3) x
10^-130 W^-8 m^16 s^-1, is printed in the text and recorded inline in the source
record; this script writes no data file. It checks the record's beam parameters
against Table 1 (p. 5, library table ``tbl_99fbcc35e0aa``) and reproduces the numbers
the record quotes from the paper's own constants:

* tau = FWHM / (2 (ln 2)^1/2) from the 164 fs FWHM, and that FWHM as 232 fs / 2^1/2
  (p. 4);
* the assumed electron-gas collision frequency nu_eg = N sigma_eg v_mean, with N the
  air density at 1 atm and 300 K and v_mean = (8 k T_e / (pi m_e))^1/2 for T_e = 0.4 eV
  and sigma_eg = 5 x 10^-16 cm^2 (p. 8: 5.18 x 10^11 s^-1);
* the O2 density n0 ~ 5.13 x 10^18 cm^-3 (p. 6) as the O2 fraction of that N.

The air's temperature ('about 300 K', p. 3), its pressure (atmospheric, not stated as
a number) and dry air's O2 fraction are the script's assumptions, named below.
"""

import math
import tomllib

from _library import DATA, Point, table

from optmat.constants import E_CHARGE, K_B, M_E, STD_ATMOSPHERE

DIRECTORY = "Sharma_2018_CountingElectronsMultiphoton"
SOURCE = "sharma_2018"
TABLE_1 = "tbl_99fbcc35e0aa"
# Table 1's rows: symbol, description, 'value +- standard error unit'; the record's
# parameters that hold each.
TABLE_ROWS = {
    "w0": ("waist_um", "waist_um_sigma", "μm"),
    "zR": ("rayleigh_length_mm", "rayleigh_length_mm_sigma", "mm"),
    "τ": ("tau_fs", "tau_fs_sigma", "fs"),
}
AUTOCORRELATION_FWHM_FS = 232.0  # p. 4
SIGMA_EG_M2 = 5e-16 * 1e-4  # p. 8, printed '5 x 10^-16 cm^-2' (cm^2)
T_E_EV = 0.4  # p. 8
NU_EG_PRINTED = 5.18e11  # s^-1, p. 8
N0_PRINTED_CM3 = 5.13e18  # p. 6
# The script's assumptions: atmospheric pressure, 300 K ('about 300 K', p. 3), and the
# O2 mole fraction of dry air, 0.2095.
AIR_PRESSURE_PA = STD_ATMOSPHERE
AIR_TEMPERATURE_K = 300.0
O2_FRACTION = 0.2095
# How closely the reproduced numbers must agree: the printed digits (0.5 in the last).
DIGITS = {"fwhm": 0.5, "tau": 0.15, "nu": 0.005e11, "n0": 0.005e18}


def _parameters() -> dict[str, object]:
    with (DATA / "sources" / f"{SOURCE}.toml").open("rb") as handle:
        record = tomllib.load(handle)
    (dataset,) = record["datasets"]
    return dataset["parameters"]


def check_table() -> None:
    """Raise unless the record's beam parameters are Table 1's, as printed."""
    parameters = _parameters()
    rows = table(DIRECTORY, TABLE_1)
    seen = set()
    for symbol, _, cell in rows:
        value, error_unit = (part.strip() for part in cell.split("±"))
        name, sigma_name, unit = TABLE_ROWS[symbol]
        error = error_unit.removesuffix(unit)
        if (str(parameters[name]), str(parameters[sigma_name])) != (value, error):
            msg = f"Table 1, {symbol}: {cell!r} is not the record's {name}"
            raise ValueError(msg)
        seen.add(symbol)
    if seen != set(TABLE_ROWS):
        msg = f"Table 1: rows {sorted(seen)}"
        raise ValueError(msg)


def reproduced() -> dict[str, float]:
    """Return the printed numbers the record quotes, computed from the paper's own."""
    parameters = _parameters()
    fwhm = AUTOCORRELATION_FWHM_FS / math.sqrt(2)
    tau = float(str(parameters["pulse_fs"])) / (2 * math.sqrt(math.log(2)))
    density = AIR_PRESSURE_PA / (K_B * AIR_TEMPERATURE_K)
    speed = math.sqrt(8 * T_E_EV * E_CHARGE / (math.pi * M_E))
    return {
        "fwhm": fwhm,
        "tau": tau,
        "nu": density * SIGMA_EG_M2 * speed,
        "n0": O2_FRACTION * density * 1e-6,
    }


def extract() -> dict[str, list[Point]]:
    """Check the record against the paper; there are no data files to write."""
    check_table()
    parameters = _parameters()
    printed = {
        "fwhm": float(str(parameters["pulse_fs"])),
        "tau": float(str(parameters["tau_fs"])),
        "nu": NU_EG_PRINTED,
        "n0": N0_PRINTED_CM3,
    }
    for name, value in reproduced().items():
        # tau: the Table's 98.6 is a mean of measurements; FWHM / (2 (ln 2)^1/2) of
        # the printed 164 fs gives 98.5, within the printed digit.
        if abs(value - printed[name]) > DIGITS[name]:
            msg = f"{name}: computed {value:.4g}, printed {printed[name]:.4g}"
            raise ValueError(msg)
    return {}


def main() -> None:
    """Run the checks and print the reproduced numbers."""
    extract()
    for name, value in reproduced().items():
        print(f"{name}: {value:.4g}")


if __name__ == "__main__":
    main()

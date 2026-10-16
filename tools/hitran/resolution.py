# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The fits' representation error inside their optical ranges, at coarser resolutions.

    uv run python tools/hitran/resolution.py [CO CO2 CH4 N2O]

OEP 3's representation bound is measured against the line-by-line reference smoothed
by a Gaussian of sigma = 1 cm^-1 (``windows.py``). Inside a fit's optical range the
bound is reported, not enforced (TC-749): near the weak overtone bands the reference's
line structure at that resolution is far beyond any pole model. This tool measures, for
each stored molecular fit and for sigma = 1, 10 and 100 cm^-1, the largest
|D_m[fit] - D_m[reference]| over the fit's qualified intervals inside its optical
range, in units of the fit's statistical one-sigma of D_m (m = 0, 2, 3), and the share
of the range where it exceeds that sigma; and the shortest wavelength where D_2, the
group-velocity dispersion's driver, misses OEP 3's tolerance, 1 % of the larger of
|D_2| and the electronic terms' |D_2| (Q18). A pulse whose bandwidth is much wider
than sigma sees the smoothed structure. The results are printed; nothing is written.
"""

import sys
import warnings

import numpy as np

from optmat.gas.models.compact import compact_intervals, statistical_sigma_au
from optmat.gas.models.fitted import fitted_models

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import windows

SIGMAS_CM = (1.0, 10.0, 100.0)


def main(gases: list[str]) -> None:
    """Print the representation error at each resolution."""
    for gas in gases:
        fit = fitted_models()[f"fit/{gas}"]
        intervals = compact_intervals(fit)
        lines = windows.load_lines(gas)
        parent = windows.load_parent(gas)
        for sigma_cm in SIGMAS_CM:
            settings = windows.Settings(sigma_cm=sigma_cm)
            reference = windows.build_reference(lines, parent, settings)
            nu = reference.grid
            wavelength = 1e-2 / nu
            inside = np.zeros(nu.size, dtype=bool)
            for low, high in intervals:
                inside |= (wavelength >= low) & (wavelength <= high)
            nu = nu[inside]
            sample = reference.sample(nu)
            omega = nu / windows.HCM
            sigma = statistical_sigma_au(fit, omega)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = {
                    0: fit.polarizability_au(omega),
                    2: fit.dispersion_au(omega, 2),
                    3: fit.dispersion_au(omega, 3),
                }
            parts = []
            for m in (0, 2, 3):
                ratio = np.abs(model[m] - sample[f"d{m}"]) / sigma[m]
                worst = float(np.max(ratio))
                share = float(np.mean(ratio > 1))
                where = 1e7 / nu[int(np.argmax(ratio))]
                parts.append(
                    f"D{m} {worst:9.3g} ({share:6.1%} > 1, worst {where:.0f} nm)"
                )
            # OEP 3's GVD tolerance: 1 % of the larger of |D_2| and the electronic
            # |D_2|, which stays finite where D_2 crosses zero.
            electronic = parent.electronic(omega)[2]
            scale = np.maximum(np.abs(sample["d2"]), np.abs(electronic))
            fails = np.abs(model[2] - sample["d2"]) > 0.01 * scale
            first = f"{1e7 / nu[fails].max():.0f} nm" if fails.any() else "nowhere"
            parts.append(
                f"D2 beyond 1 % of max(|D2|, |D2 el|) from {first} "
                f"({float(np.mean(fails)):.1%} of the points)"
            )
            span = ", ".join(f"{a * 1e9:.0f}-{b * 1e9:.0f}" for a, b in intervals)
            print(
                f"{gas:4s} sigma {sigma_cm:5.0f} cm-1 over {span} nm: "
                + "; ".join(parts)
            )


if __name__ == "__main__":
    main(sys.argv[1:] or list(windows.GASES))

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
CO2 from the ultraviolet to its infrared windows
================================================

The refitted CO2 model (OEP 3) from 0.25 to 12 µm at 1 atm and 296 K. Top: n - 1.
Bottom: the group-velocity dispersion's driver,
:math:`|D_2| = |d^2(\\omega\\alpha)/d\\omega^2|`, proportional to :math:`\\beta_2` in
the dilute limit, with its statistical one-sigma uncertainty and the representation
bound of the infrared terms. The shaded regions are the gaps between the qualified
intervals, where the gas absorbs or the bands' Gaussian rules do not meet the
tolerance. The model is not to be used there, and warns.

Inside the optical range (to 1.4 µm) the bound is reported, not enforced (TC-749): near
1 µm it exceeds the dispersion itself, from the weak combination bands' line structure,
which no pole model holds. Beyond the range each window is qualified because the bound
is within the tolerance: a tenth of the fit's statistical sigma in the polarizability,
and a hundredth of the larger of the full and the electronic dispersion. There the
fixed infrared terms carry no statistical uncertainty, and the bound, not the sigma,
sets how well the dispersion is known.
"""

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

from optmat.conditions import State
from optmat.constants import AU_POLARIZABILITY
from optmat.exceptions import ModelRangeWarning
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.models.fitted import fitted_models

co2 = fitted_models()["fit/CO2"]
state = State(pressure=101325.0, temperature=296.0)
fig, (top, bottom) = plt.subplots(
    2, 1, figsize=(7, 6), sharex=True, layout="constrained"
)
for low, high in co2.intervals:
    wavelength = np.geomspace(low * 1.0005, high * 0.9995, 120)
    omega = HARTREE_NM * 1e-9 / wavelength  # hartree
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        top.loglog(wavelength * 1e6, co2.refractivity(wavelength, state), "C0")
        d2 = np.abs(co2.dispersion_au(omega, 2))
        sigma = co2.dispersion_sigma_relative(omega, 2) * d2
        bound = co2.representation_bound_au(omega, 2)
    bottom.loglog(wavelength * 1e6, d2, "C0", label="|D2|" if low < 1e-6 else None)
    bottom.loglog(
        wavelength * 1e6, sigma, "C1", label="statistical sigma" if low < 1e-6 else None
    )
    bottom.loglog(
        wavelength * 1e6,
        bound,
        "C2",
        label="representation bound" if low < 1e-6 else None,
    )
edges = [b for _, b in co2.intervals[:-1]], [a for a, _ in co2.intervals[1:]]
for axis in (top, bottom):
    for start, stop in zip(*edges, strict=True):
        axis.axvspan(start * 1e6, stop * 1e6, color="0.9", lw=0)
    axis.axvspan(co2.intervals[-1][1] * 1e6, 12.5, color="0.9", lw=0)
    axis.grid(which="both", alpha=0.3)
top.set_ylabel("n - 1 (1 atm, 296 K)")
bottom.set_ylabel(f"a.u. (alpha in units of {AU_POLARIZABILITY:.3g} C m²/V)")
bottom.set_xlabel("vacuum wavelength (µm)")
bottom.set_xlim(0.25, 12.5)
bottom.legend(fontsize="small")
plt.show()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
The dispersion of a helium-argon fill against pressure
======================================================

A 50/50 helium-argon fill at 294 K, from 1 to 60 bar, at 800 nm. The reference is
optmat's mixture law (OEP 4): the mixture's equation of state for the densities and
Lorentz-Lorenz of the whole mixture, with its pair terms. Each curve leaves one piece
out and shows the relative change in :math:`n - 1` (left) and in the group-velocity
dispersion (right):

- **dilute law**: :math:`\\chi = \\sum_s N_s\\alpha_s/\\varepsilon_0`, no local field;
- **no pair terms**: Lorentz-Lorenz without the refractivity virial terms;
- **partial-pressure density**: each gas's own density at its partial pressure, in place
  of the mixture's equation of state;
- **multi-fluid density**: CoolProp's multi-fluid model in place of the default virial
  fallback. Its cross virials miss the measured ones (TC-779), so the gap between the
  two methods is a measure of how well the density is known.

The shaded band is the declared uncertainty of the unlike-pair estimate (OEP 4,
section 7.2).
"""

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np

from optmat.conditions import State
from optmat.constants import EPS0, C
from optmat.exceptions import ModelRangeWarning
from optmat.gas.eos import equation_of_state
from optmat.gas.mixture import mixture_optics
from optmat.gas.mixture_eos import mixture_equation_of_state
from optmat.gas.models.registry import default_model

WAVELENGTH = 800e-9
TEMPERATURE = 294.0
X = {"He": 0.5, "Ar": 0.5}
optics = mixture_optics(list(X))
with warnings.catch_warnings():
    warnings.simplefilter("ignore", ModelRangeWarning)
    default = mixture_equation_of_state(X)
    multi = mixture_equation_of_state(X, method="multi-fluid")


def refractivity(kind, wavelength, densities):
    """Return n - 1 by one variant of the law."""
    if kind in ("full", "pp density", "multi-fluid density"):
        return optics.refractivity(wavelength, densities)
    alpha = {s: default_model(s).polarizability(wavelength) for s in densities}
    if kind == "dilute":
        chi = sum(densities[s] * alpha[s] / EPS0 for s in densities)
    else:  # no pair terms
        y = sum(densities[s] * alpha[s] / (3 * EPS0) for s in densities)
        chi = 3 * y / (1 - y)
    return chi / (1 + np.sqrt(1 + chi))


def densities_for(kind, pressure):
    """Return the species' densities by the variant's density rule."""
    state = State(pressure=pressure, temperature=TEMPERATURE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        if kind == "pp density":
            return {
                s: equation_of_state(s).number_density(
                    State(pressure=x * pressure, temperature=TEMPERATURE)
                )
                for s, x in X.items()
            }
        if kind == "multi-fluid density":
            return multi.partial_densities(state)
        return default.partial_densities(state)


def gvd(kind, densities):
    """Return beta_2 (fs^2/m) by a centred second difference in frequency."""
    omega = 2 * np.pi * C / WAVELENGTH
    step = omega * 1e-3
    beta = [
        (omega + j * step)
        / C
        * refractivity(kind, 2 * np.pi * C / (omega + j * step), densities)
        for j in (-1, 0, 1)
    ]
    return (beta[0] - 2 * beta[1] + beta[2]) / step**2 * 1e30


pressures = np.geomspace(1e5, 6e6, 30)  # Pa
kinds = ("dilute", "no pair terms", "pp density", "multi-fluid density")
fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
reference_n, reference_gvd, band = [], [], []
for p in pressures:
    d = densities_for("full", p)
    reference_n.append(float(refractivity("full", WAVELENGTH, d)))
    reference_gvd.append(float(gvd("full", d)))
    band.append(float(optics.refractivity_uncertainty(WAVELENGTH, d)))
reference_n = np.array(reference_n)
reference_gvd = np.array(reference_gvd)
for kind in kinds:
    n, b2 = [], []
    for p in pressures:
        d = densities_for(kind, p)
        n.append(float(refractivity(kind, WAVELENGTH, d)))
        b2.append(float(gvd(kind, d)))
    left.semilogx(pressures / 1e5, np.array(n) / reference_n - 1, label=kind)
    right.semilogx(pressures / 1e5, np.array(b2) / reference_gvd - 1, label=kind)
left.fill_between(
    pressures / 1e5,
    -np.array(band) / reference_n,
    np.array(band) / reference_n,
    color="0.85",
    label="unlike-pair estimate",
)
for axis, title in ((left, "n - 1"), (right, "GVD")):
    axis.axhline(0, color="k", lw=0.5)
    axis.set_xlabel("pressure (bar), 294 K")
    axis.set_ylabel(f"relative change in {title}")
    axis.grid(which="both", alpha=0.3)
left.legend(fontsize="small")
plt.show()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""
Ionisation: the noble gases against all-electron rates at 390 nm
================================================================

Van der Hart's R-matrix calculations of the whole atom, with one state of the ion
(:func:`optmat.gas.ionisation_validation.rmatrix_rates`): R-matrix Floquet rates of He,
Ne and Ar (2006, read from the image of its Fig. 6), and time-dependent R-matrix rates
of Ar (2007, its Table I and the vector drawing of its Fig. 4). Values read from
figures follow data rule 8 (TC-699 to TC-701).

**Left:** the references (markers) and the default's cycle-averaged rate (lines). The
references resolve the atoms' resonances, which the model does not have.

**Right:** the default over each reference against the Keldysh parameter, with
Parker et al.'s two-electron rate of helium over the same Floquet rate (grey). Helium's
two all-electron methods agree away from the resonance at 1.5e14 W/cm². The current
default includes the frequency-dependent Coulomb factor adopted in TC-848; the
larger excesses reported in TC-701 describe the earlier static-factor model.
See :doc:`/user-guide/ionisation-laser` for the current ratios and the distinction
between resonant structure and the smooth rate model.
"""

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np
from _helpers import COLOURS, MARKERS

from optmat.constants import AU_ANGULAR_FREQUENCY
from optmat.exceptions import ModelRangeWarning
from optmat.gas import ionisation_validation as validation
from optmat.gas.ionisation_rates import ionisation_rate_model

GASES = {"He": COLOURS[0], "Ne": COLOURS[1], "Ar": COLOURS[2]}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
left, right = axes
with warnings.catch_warnings():
    warnings.simplefilter("ignore", ModelRangeWarning)
    for gas, colour in GASES.items():
        model = ionisation_rate_model(gas, 390e-9)
        for k, comparison in enumerate(validation.rmatrix_rates(gas)):
            intensity = comparison.x / 1e4
            left.semilogy(
                intensity,
                comparison.expected * AU_ANGULAR_FREQUENCY,
                MARKERS[k],
                color=colour,
                mfc="none" if k else colour,
                ms=4,
                label=f"{gas}: {comparison.label}",
            )
            gamma = model.keldysh(validation.amplitude(comparison.x))
            right.semilogy(
                gamma,
                comparison.ratio,
                MARKERS[k],
                color=colour,
                mfc="none" if k else colour,
                ms=4,
            )
        dense = np.geomspace(0.2e18, 3e18, 200)
        left.semilogy(
            dense / 1e4,
            model.cycle_averaged_rate(validation.amplitude(dense)),
            color=colour,
            lw=1.2,
        )
    (helium,) = validation.rmatrix_rates("He")
    gamma = ionisation_rate_model("He", 390e-9).keldysh(validation.amplitude(helium.x))
    right.semilogy(
        gamma,
        validation.parker_rate_au(helium.x, 390) / helium.expected,
        "x",
        color="#86858a",
        label="He: Parker et al. over the Floquet rate",
    )
left.set_xlabel("intensity (W/cm$^2$)")
left.set_ylabel("cycle-averaged rate (1/s)")
left.set_title("390 nm: references (markers), the default (lines)")
left.legend(fontsize=7)
right.axhline(1.0, color="#0b0b0b", lw=0.8)
right.set_xlabel("Keldysh parameter")
right.set_ylabel("default / reference")
right.set_title("the default over the all-electron rates")
right.legend(fontsize=8)
plt.show()

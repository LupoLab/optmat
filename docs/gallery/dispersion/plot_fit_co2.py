# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Carbon dioxide: the fitted model
================================

The CO2 model (``fit/CO2``): two effective electronic oscillators fitted to Old,
Gentili & Peck (480-1600 nm, level free), the visible data, the 633 and 1542 nm
metrology, Hohm's 325-633 nm dispersion, Shelton's N2/CO2 ratios with Marchetti's
CO2-laser values as validators, with HITRAN's vibrational bands below 7000 cm^-1
fixed (reduced
chi-squared 4.01; TC-797).
Shown over the model's main range, 250 nm to 1.4 um; seven additional qualified
infrared intervals extend to 8.941 um, with gaps listed in :doc:`/gas_models`.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("CO2", xlim=(250e-9, 1.4e-6))

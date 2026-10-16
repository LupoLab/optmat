# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Methane: the fitted model
=========================

The methane model (``fit/CH4``): two effective electronic oscillators fitted to the
Cuthbertsons' (1920) and Watson & Ramaswamy's visible dispersion (level free),
Achtermann's 633 nm and Burns et al.'s 458-633 nm values, Zhao et al.'s 1566 nm value
and Shelton's N2/CH4 ratios, with HITRAN's vibrational bands below 4500 cm^-1 fixed
(reduced chi-squared 0.68; TC-797). Shown from 240 nm to 1.4 um,
across the absorption gaps; the qualified intervals (240 to 1123 nm and the infrared
windows) are listed in :doc:`/gas_models` (TC-797).
"""

# %%
from _helpers import fit_overview

fig = fit_overview("CH4", xlim=(240e-9, 1.4e-6))

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Carbon monoxide: the fitted model
=================================

The CO model (``fit/CO``): two effective electronic oscillators fitted to Smith et al.'s
167-288 nm interferometry, the Cuthbertsons' 1920 visible dispersion (level free),
Watson & Ramaswamy's green line and Zhao et al.'s 2026 refractometry at 1566 nm, with
HITRAN's fundamental and first overtone fixed (reduced chi-squared 0.37).
Shown over the optical range. The main qualified interval extends from 167 to
2288 nm, with two further infrared windows; see :doc:`/gas_models` (TC-797).
"""

# %%
from _helpers import fit_overview

fig = fit_overview("CO", xlim=(167e-9, 2e-6))

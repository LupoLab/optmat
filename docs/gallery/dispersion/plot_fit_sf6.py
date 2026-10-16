# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Sulphur hexafluoride: the fitted model
======================================

The SF6 model (``fit/SF6``): two effective electronic oscillators fitted to Swift's,
Vukovic's, Achtermann's and Burns et al.'s values and Shelton's N2/SF6 ratios, with
HITRAN's measured infrared cross-sections as fixed bands (reduced chi-squared 4.1; the
two Shelton sets prefer the bands times 1.036, TC-146). The bands give the gas a
zero-dispersion wavelength near 1.2 um. Shown from 230 nm to 5 um.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("SF6", xlim=(230e-9, 5e-6))

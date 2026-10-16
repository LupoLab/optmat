# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Nitrogen: the fitted model
==========================

The nitrogen model (``fit/N2``): effective oscillators fitted to Griesmann & Burnett's
vacuum ultraviolet (level free), Peck & Khanna's visible and near infrared, and the 633
nm, microwave, dielectric and 10.57 um values agree to their uncertainties (reduced
chi-squared 1.1).
"""

# %%
from _helpers import fit_overview

fig = fit_overview("N2", xlim=(140e-9, 20e-6))

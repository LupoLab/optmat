# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Carbon dioxide: the curated data
================================

Every direct carbon dioxide dataset (a refractivity, polarizability or molar
refractivity) turned into the refractivity of the real gas at 0 °C and 1 atm, with the
published formulas (solid) and refractiveindex.info's entries (dotted), each over the
range its source supports. The lower panel shows the deviations from Old, Gentili &
Peck's formula.
"""

# %%
from _helpers import data_overview

fig = data_overview("CO2", "old_1971/sellmeier", xlim=(1.5e-07, 2e-05))

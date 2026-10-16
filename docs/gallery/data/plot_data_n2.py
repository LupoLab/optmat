# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Nitrogen: the curated data
==========================

Every direct nitrogen dataset (a refractivity, polarizability or molar refractivity)
turned into the refractivity of the real gas at 0 °C and 1 atm, with the published
formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from Peck & Khanna's formula at 0 °C (Eq. 4)."""

# %%
from _helpers import data_overview

fig = data_overview("N2", "peck_1966/eq4", xlim=(1e-07, 2e-05))

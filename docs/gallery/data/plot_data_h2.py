# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Hydrogen: the curated data
==========================

Every direct hydrogen dataset (a refractivity, polarizability or molar refractivity)
turned into the refractivity of the real gas at 0 °C and 1 atm, with the published
formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from Peck & Huang's two-term formula (Eq. 7)."""

# %%
from _helpers import data_overview

fig = data_overview("H2", "peck_1977/eq7", xlim=(9e-08, 2e-05))

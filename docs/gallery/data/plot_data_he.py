# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Helium: the curated data
========================

Every direct helium dataset (a refractivity, polarizability or molar refractivity)
turned into the refractivity of the real gas at 0 °C and 1 atm, with the published
formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from Mansfield & Peck's two-pole fit to their
own theoretical values (Eq. 8, corrected), which is smooth down to the first
resonance."""

# %%
from _helpers import data_overview

fig = data_overview("He", "mansfield_peck_1969/eq8", xlim=(4e-08, 2e-05))

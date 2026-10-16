# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Xenon: the curated data
=======================

Every direct xenon dataset (a refractivity, polarizability or molar refractivity)
turned into the refractivity of the real gas at 0 °C and 1 atm, with the published
formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from Bideau-Mehu et al.'s formula (corrected),
valid on both sides of the 146.9 nm line."""

# %%
from _helpers import data_overview

fig = data_overview("Xe", "bideau_mehu_1981/sellmeier_Xe", xlim=(1.2e-07, 2e-05))

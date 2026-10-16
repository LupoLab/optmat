# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Krypton: the curated data
=========================

Every direct krypton dataset (a refractivity, polarizability or molar refractivity)
turned into the refractivity of the real gas at 0 °C and 1 atm, with the published
formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from Bideau-Mehu et al.'s formula (known to lie
0.1-0.7 % below the paper's own column)."""

# %%
from _helpers import data_overview

fig = data_overview("Kr", "bideau_mehu_1981/sellmeier_Kr", xlim=(1e-07, 2e-05))

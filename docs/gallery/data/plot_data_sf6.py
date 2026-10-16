# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Sulphur hexafluoride: the curated data
======================================

Every direct sulphur hexafluoride dataset (a refractivity, polarizability or molar
refractivity) turned into the refractivity of the real gas at 0 °C and 1 atm, with
the published formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from He et al.'s formula (Eq. 8)."""

# %%
from _helpers import data_overview

fig = data_overview("SF6", "he_2021/SF6", xlim=(1.5e-07, 2e-05))

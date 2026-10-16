# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Nitrous oxide: the curated data
===============================

Every direct nitrous oxide dataset (a refractivity, polarizability or molar
refractivity) turned into the refractivity of the real gas at 0 °C and 1 atm, with
the published formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports.
The lower panel shows the deviations from He et al.'s formula (Eq. 9)."""

# %%
from _helpers import data_overview

fig = data_overview("N2O", "he_2021/N2O", xlim=(1.5e-07, 2e-05))

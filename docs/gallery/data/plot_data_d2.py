# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Deuterium: the curated data
===========================

Every direct deuterium dataset (a refractivity, polarizability or molar
refractivity) turned into the refractivity of the real gas at 0 °C and 1 atm, with
the published formulas (solid), the recorded ``luna_physdata`` formulas (dashed) and
refractiveindex.info's entries (dotted), each over the range its source supports. No
published formula for this gas is implemented yet, so there is no deviation panel."""

# %%
from _helpers import data_overview

fig = data_overview("D2", None, xlim=(9e-08, 2e-05))

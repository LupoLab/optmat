# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Deuterium: the fitted model
===========================

The deuterium model (``fit/D2``) uses Raj et al.'s (2018) ab initio average for
normal deuterium and Shelton's D2/H2 and H2/D2 ratios. Published refractometry is
shown as independent evidence where identified as a validator. See
:doc:`/gas_models` for the current data selection and :doc:`/accuracy` for limitations.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("D2", xlim=(105e-9, 20e-6))

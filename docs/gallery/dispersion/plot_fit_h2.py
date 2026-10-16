# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Hydrogen: the fitted model
==========================

The hydrogen model (``fit/H2``) is based on Raj et al.'s (2018) ab initio
average for normal hydrogen, with a 0.1 % allowance, the dispersion shapes of Peck
and Hung and of Smith et al., and Shelton's He/H2 phase-matching ratio. The current
model's statistical uncertainty and the remaining absolute-scale evidence are
summarised in :doc:`/gas_models` and :doc:`/accuracy`.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("H2", xlim=(105e-9, 20e-6))

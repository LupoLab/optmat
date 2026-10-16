# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Oxygen: the fitted model
========================

The oxygen model (``fit/O2``): effective oscillators fitted to Smith, Huber &
Parkinson's ultraviolet above the Schumann-Runge bands, Egan et al.'s 2026 values at 633
and 1542 nm, Birch at 633 nm, Hohm's 325-633 nm dispersion, Swift and Marchetti; Zhang's
800 nm value (A band) and Kumar & Thakkar's sums are validators.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("O2", xlim=(180e-9, 20e-6))

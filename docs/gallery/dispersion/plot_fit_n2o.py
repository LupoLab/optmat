# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Nitrous oxide: the fitted model
===============================

The N2O model (``fit/N2O``): two effective electronic oscillators fitted to Hohm's
325-633 nm dispersion (1e-5 relative) and the Cuthbertsons' (1914) visible dispersion
(both level free) and the 633 nm metrology of Birch and Egan et al., with HITRAN's
vibrational bands below 4500 cm^-1 fixed and the Thomas-Reiche-Kuhn sum as a prior
(reduced chi-squared 3.01: the two dispersions differ in shape by about 1e-4;
TC-797). The main qualified interval is 320-1516 nm, with three additional
infrared windows; see :doc:`/gas_models`.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("N2O", xlim=(320e-9, 2e-6))

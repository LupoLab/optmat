# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Neon: the fitted model
======================

The neon model (``fit/Ne``): the 73.6 nm and 74.4 nm resonance lines at their NIST
energies and two effective oscillators, with the magnetic susceptibility :math:`\\chi_0
= -8.48\\times10^{-5}` a.u. The level comes from Lesiuk, Przybytek & Jeziorski's ab
initio Cauchy moments and Gaiser & Fellmuth's 2018 dielectric static polarizability (via
Lesiuk; their 2010 value is 16 ppm higher); the shape from the Cuthbertsons' visible
data (level free: their standard state is ambiguous, and it comes out 0.3 % high) and
Bideau-Mehu's vacuum ultraviolet.

The fit's reduced chi-squared is 2.5: the experiments' dispersion is about 0.6 %
stronger, relative to the static term, than the ab initio S(-4) allows (TC-92),
as the refractivity-constrained distributions of Kumar & Thakkar and of Langhoff &
Karplus also show; the band is scaled by the reduced chi-squared. Birch (+0.14 %),
Orcutt & Cole (+0.33 %) and Leonard's compilation (+0.27 %) lie high, like the classical
helium values.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("Ne", xlim=(70e-9, 20e-6))

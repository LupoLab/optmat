# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Helium: the fitted model
========================

The helium model (``fit/He``): the 1s2p (58.4 nm) and 1s3p (53.7 nm) lines at their NIST
energies and four effective oscillators, with the magnetic susceptibility :math:`\\chi_0
= -2.12\\times10^{-5}` a.u. that makes the electric polarizability the refractive one.
The backbone is Puchalski et al.'s ab initio dispersion coefficients of 4He, with every
relativistic and QED correction (the static term to 1e-7), and Schmidt et al.'s
microwave static polarizability: the model reproduces both to their uncertainties, so
its one-sigma band is below 1 ppm from the infrared to 110 nm and 7 ppm at 76 nm.

Against it, the classical refractometry scatters by +-0.3 %: Mansfield & Peck +0.23 %
(fitted for its shape, which departs from theory by up to 1.1e-4), Birch +0.32 %,
Cuthbertson 1910 +0.24 %, Newell's 47.7 GHz value -0.30 %, Smith 1976 within its 0.2 %.
Liu's infinite-mass, nonrelativistic ab initio values lie 4.1e-4 below at long
wavelengths, as the mass and relativistic corrections require, the difference growing
towards the resonance as its shifted pole implies. Ermolov's formula, read at the 1 bar
its paper states, is 1.3 % high; read at 1 atm it is within 0.1 %.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("He", xlim=(55e-9, 20e-6))

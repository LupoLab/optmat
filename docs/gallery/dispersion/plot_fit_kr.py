# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Krypton: the fitted model
=========================

The krypton model (``fit/Kr``): the 116.5 nm and 123.6 nm resonance lines at their NIST
energies and two effective oscillators. The only ab initio dispersion in the library
(Pawlowski et al.'s CC3 moments) is not converged, so the fit rests on data: the
dielectric static polarizabilities of Orcutt & Cole and of Huot & Bose (via Kumar &
Thakkar; 3.6 sigma from the others), Birch's and Hohm & Kerl's 633 nm values, Burns et
al.'s 458-633 nm molar refractivities, Smith, Huber & Parkinson's 28 ultraviolet points
(0.06 %), Bideau-Mehu's revised (1976) and vacuum-ultraviolet (1981) values and the
Cuthbertsons' visible shape (level free). Leonard's absolute values lie 0.2 % low
and his recommended one 0.04 % high. The magnetic term (about -25 ppm) is neglected.
Chopra & Heddle's Lyman-alpha value lies between the two resonance lines, where the sign
of n - 1 is not measured.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("Kr", xlim=(110e-9, 20e-6))

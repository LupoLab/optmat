# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""
Xenon: the fitted model
=======================

The xenon model (``fit/Xe``): the 146.96 nm and 129.56 nm resonance lines at their NIST
energies and two effective oscillators, one of which settles at 121.6 nm, between the
strong 119.2 nm and the 125.0 nm lines. The level comes from Huot & Bose's dielectric
static polarizability (via Kumar & Thakkar), Burns et al.'s 458-633 nm molar
refractivities and Hohm & Kerl's 633 nm value; the shape from Abjean et al.'s and
Bideau-Mehu et al.'s ultraviolet and vacuum ultraviolet, which crosses the 146.96 nm
line (140.4 nm is between the two lines), and the Cuthbertsons' visible data (level
free). The fitted strength of the 146.96 nm line, 0.2688, is Bideau-Mehu's measured
0.268(8). Leonard's recommended value is 0.13 % high. The magnetic term (about -25 ppm)
is neglected.
"""

# %%
from _helpers import fit_overview

fig = fit_overview("Xe", xlim=(125e-9, 20e-6))

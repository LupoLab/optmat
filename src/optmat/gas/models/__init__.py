# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Dispersion models: published formulas and the models fitted to the curated data.

Every model gives the per-particle polarizability as a function of vacuum wavelength
(:class:`~optmat.gas.models.base.Model`); the refractive index at any state
follows by Lorentz-Lorenz. Published formulas are defined, with their printed
coefficients and any corrections, in the source records' ``[[formulas]]`` entries and
built by :mod:`~optmat.gas.models.formulas`.
"""

__all__: list[str] = []

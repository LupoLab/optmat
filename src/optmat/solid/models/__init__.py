# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Models of solids' dielectric functions: the interface, formulas and the registry.

- :mod:`optmat.solid.models.base`: the interface
  (:class:`~optmat.solid.models.base.Model`);
- :mod:`optmat.solid.models.formulas`: the published formulas and Franta et al.'s
  universal dispersion model, built from the curated records;
- :mod:`optmat.solid.models.tabulated`: curated tables of optical constants joined by
  range (a metal's complex index);
- :mod:`optmat.solid.models.registry`: every model by name, and the recommended one of
  each material with the reason.
"""

__all__: list[str] = []

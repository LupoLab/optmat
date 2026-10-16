# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The composition of gas mixtures, from the curated records: dry air (TC-548).

A mixture is not a material: its make-up is a ``mole_fraction`` dataset per constituent
(family ``composition``, TC-543), each naming the mixture as its source defines it. This
module reads them. The recommended composition of dry air is the CIPM-2007 one (Picard,
Davis, Gläser & Fujii, Metrologia 45, 149 (2008), Table 1): N2, O2, Ar and CO2 from the
CIPM's own evaluation, and eight trace constituents that the CIPM adopted from the US
Standard Atmosphere 1976.

The fractions are those of the printed table and sum to 1 - 2.9e-6 (the table's
rounding). A propagation code that models only some constituents renormalises them
itself and says so; this module does not.
"""

import math
from dataclasses import dataclass
from typing import Final

from optmat.datasets import datasets
from optmat.exceptions import OptmatError

__all__ = ["MIXTURES", "REDUCED", "Constituent", "composition"]

MIXTURES: Final = {"dry air": "dry air (CIPM-2007)"}
"""The recommended composition of each mixture: its name in the curated datasets."""

REDUCED: Final = {
    "dry air (N2, O2, Ar)": ("dry air (CIPM-2007)", ("N2", "O2", "Ar")),
    "dry air (N2, O2)": ("dry air (CIPM-2007)", ("N2", "O2")),
}
"""Named compositions that keep some constituents of a curated one, renormalised.

``"dry air (N2, O2, Ar)"`` keeps 99.957 % of CIPM-2007's molecules: the air of
:func:`optmat.gas.raman_response.raman_response` and the default air model
(OEP 1, section 4.4). ``"dry air (N2, O2)"`` is the argon-free reduction. Both leave
out CIPM-2007's CO2 (400 umol/mol) and its trace constituents."""


@dataclass(frozen=True)
class Constituent:
    """One constituent of a mixture.

    Attributes
    ----------
    material : str
        The constituent gas's material id.
    mole_fraction : float
        Its amount-of-substance fraction, mol/mol.
    sigma : float
        The standard uncertainty of the fraction where the source states one, else NaN.
    source : str
        The curated dataset, ``<source>/<id>``.
    primary : bool
        False where the source reprints the fraction from another (the trace
        constituents of CIPM-2007 dry air, from the US Standard Atmosphere 1976).
    """

    material: str
    mole_fraction: float
    sigma: float
    source: str
    primary: bool


def composition(mixture: str = "dry air") -> tuple[Constituent, ...]:
    """Return the constituents of a mixture, the most abundant first.

    Parameters
    ----------
    mixture : str, optional
        A key of :data:`MIXTURES` (the recommended composition, as printed), a
        mixture as the curated datasets name it (``"dry air (CIPM-2007)"``), or a
        reduced composition of :data:`REDUCED`, whose kept constituents' fractions
        (and their uncertainties) are divided by their sum.

    Returns
    -------
    tuple of Constituent
        The constituents.

    Raises
    ------
    OptmatError
        If no curated dataset describes the mixture.

    Examples
    --------
    >>> air = composition("dry air")
    >>> [(c.material, c.mole_fraction) for c in air[:4]]
    [('N2', 0.780848), ('O2', 0.20939), ('Ar', 0.009332), ('CO2', 0.0004)]
    """
    if mixture in REDUCED:
        base, kept = REDUCED[mixture]
        chosen = [c for c in composition(base) if c.material in kept]
        total = sum(c.mole_fraction for c in chosen)
        return tuple(
            Constituent(
                material=c.material,
                mole_fraction=c.mole_fraction / total,
                sigma=c.sigma / total,
                source=c.source,
                primary=c.primary,
            )
            for c in chosen
        )
    name = MIXTURES.get(mixture, mixture)
    found = []
    for dataset in datasets():
        if dataset.quantity != "mole_fraction":
            continue
        if dataset.parameters.get("mixture") != name:
            continue
        (value,) = dataset.scaled_values()
        (sigma,) = dataset.standard_uncertainty()
        found.append(
            Constituent(
                material=dataset.material,
                mole_fraction=float(value),
                sigma=float(sigma) if math.isfinite(sigma) else math.nan,
                source=dataset.key,
                primary=dataset.use == "primary",
            )
        )
    if not found:
        msg = f"no curated composition of {mixture!r}"
        raise OptmatError(msg)
    return tuple(sorted(found, key=lambda c: -c.mole_fraction))

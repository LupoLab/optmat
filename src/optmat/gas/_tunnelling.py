# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Tunnelling ionisation rates, for validating the ionisation parameters only.

optmat supplies the parameters of the ionisation rates (TC-439): the ionisation energy,
the structure coefficients of the ionising orbital and the polarizabilities. The rate
formulas belong to the propagation code. This private module evaluates the standard
static-field tunnelling rate only so that the curated parameters can be compared with
ab initio static-field rates; it is not public API.

The rate is that of the molecular ADK theory of X. M. Tong, Z. X. Zhao and C. D. Lin,
Phys. Rev. A 66, 033402 (2002), which reduces to the atomic ADK rate for one partial
wave. The asymptotic wave function of the ionising orbital is written (their eqs. 1-2)

.. math::

   \Psi^m(\mathbf r) = \sum_l C_l F_l(r) Y_{lm}(\hat{\mathbf r}), \qquad
   F_l(r \to \infty) \approx r^{Z_c/\kappa - 1} e^{-\kappa r},

with :math:`\kappa = \sqrt{2 I_p}` and :math:`Z_c` the charge of the ion, and the static
rate with the field along the molecular axis is their eq. (7),

.. math::

   w(F) = \frac{B^2(m)}{2^{|m|} |m|!} \frac{1}{\kappa^{2 Z_c/\kappa - 1}}
          \left(\frac{2\kappa^3}{F}\right)^{2 Z_c/\kappa - |m| - 1}
          e^{-2\kappa^3 / 3F},
   \qquad B(m) = \sum_l C_l Q(l, m),

:math:`Q(l, m) = (-1)^m \sqrt{(2l + 1)(l + |m|)! / (2 (l - |m|)!)}` (eq. 4). For a
molecule at the Euler angles :math:`\mathbf R` to the field, :math:`B(m')` takes the
rotation matrix :math:`D^l_{m'm}(\mathbf R)` (eqs. 8-9); averaged over orientations the
cross terms between partial waves vanish, and
:math:`\langle B^2(m') \rangle = \sum_l C_l^2 Q^2(l, m') / (2l + 1)`.

Everything is in atomic units.
"""

import math
from collections.abc import Mapping

__all__ = [
    "aligned_structure_factors",
    "isotropic_structure_factors",
    "q_factor",
    "static_rate_au",
]


def q_factor(l: int, m: int) -> float:  # noqa: E741 - l is the partial wave's symbol
    """Return Tong et al.'s normalisation factor Q(l, m) (their eq. 4).

    Parameters
    ----------
    l : int
        Orbital angular momentum, >= 0.
    m : int
        Magnetic quantum number along the field, ``|m| <= l``.

    Returns
    -------
    float
        ``(-1)^m sqrt((2l + 1) (l + |m|)! / (2 (l - |m|)!))``.

    Raises
    ------
    ValueError
        If ``|m| > l``.

    Examples
    --------
    >>> f"{q_factor(0, 0):.6f}"
    '0.707107'
    """
    if abs(m) > l:
        msg = f"|m| = {abs(m)} exceeds l = {l}"
        raise ValueError(msg)
    return (-1) ** m * math.sqrt(
        (2 * l + 1) * math.factorial(l + abs(m)) / (2 * math.factorial(l - abs(m)))
    )


def aligned_structure_factors(
    coefficients: Mapping[int, float], m: int
) -> dict[int, float]:
    """Return B^2(m) for an orbital aligned with the field (Tong et al., eq. 6).

    Parameters
    ----------
    coefficients : mapping of int to float
        ``C_l`` by partial wave ``l``, in Tong et al.'s normalisation (eq. 2).
    m : int
        The orbital's magnetic quantum number along the molecular axis.

    Returns
    -------
    dict of int to float
        ``{m: B(m)^2}``.
    """
    b = sum(c * q_factor(ell, m) for ell, c in coefficients.items() if ell >= abs(m))
    return {m: b * b}


def isotropic_structure_factors(
    coefficients: Mapping[int, float],
) -> dict[int, float]:
    """Return <B^2(m')> averaged over the orientations of the molecule.

    The average of ``|sum_l C_l D^l_{m'm}(R) Q(l, m')|^2`` over the Euler angles is
    ``sum_l C_l^2 Q(l, m')^2 / (2l + 1)``, by the orthogonality of the rotation
    matrices; it does not depend on the orbital's own ``m``.

    Parameters
    ----------
    coefficients : mapping of int to float
        ``C_l`` by partial wave ``l``, in Tong et al.'s normalisation (eq. 2).

    Returns
    -------
    dict of int to float
        ``<B^2(m')>`` for every ``m'`` from ``-l_max`` to ``l_max``.
    """
    l_max = max(coefficients)
    return {
        mp: sum(
            c * c * q_factor(ell, mp) ** 2 / (2 * ell + 1)
            for ell, c in coefficients.items()
            if ell >= abs(mp)
        )
        for mp in range(-l_max, l_max + 1)
    }


def static_rate_au(
    field_au: float,
    ionization_energy_au: float,
    structure_factors: Mapping[int, float],
    charge: float = 1.0,
) -> float:
    """Return the static-field tunnelling rate (Tong et al., eqs. 7 and 9).

    Parameters
    ----------
    field_au : float
        Field strength, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy of the orbital, hartree, > 0.
    structure_factors : mapping of int to float
        ``B^2(m')`` by ``m'``, from :func:`aligned_structure_factors` or
        :func:`isotropic_structure_factors`.
    charge : float, optional
        Charge ``Z_c`` of the ion left behind; 1 for a neutral.

    Returns
    -------
    float
        Rate per molecule in the orbital, atomic units (per atomic unit of time).

    Examples
    --------
    Hydrogen 1s (``C_0 = 2``) gives the exact asymptotic static rate
    ``4 / F exp(-2 / 3F)``:

    >>> w = static_rate_au(0.05, 0.5, aligned_structure_factors({0: 2.0}, 0))
    >>> f"{w / (4 / 0.05 * math.exp(-2 / 0.15)):.12f}"
    '1.000000000000'
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    power = 2 * charge / kappa
    rate = 0.0
    for m, b2 in structure_factors.items():
        rate += (
            b2
            / (2 ** abs(m) * math.factorial(abs(m)))
            * kappa ** (1 - power)
            * (2 * kappa**3 / field_au) ** (power - abs(m) - 1)
        )
    return rate * math.exp(-2 * kappa**3 / (3 * field_au))

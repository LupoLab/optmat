# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Convert between the per-particle polarizability and the refractive index.

The package stores gas dispersion as the dynamic dipole polarizability ``alpha(omega)``
of one atom or molecule, a property of the particle that does not depend on the state
of the gas. The refractive index at number density ``N`` follows from the
Lorentz-Lorenz (Clausius-Mossotti) relation

.. math::

   \frac{n^2 - 1}{n^2 + 2} = \frac{N \alpha}{3 \varepsilon_0} \equiv y,
   \qquad n^2 - 1 = \chi = \frac{3 y}{1 - y},

with ``alpha`` in SI (C m^2/V). The dilute limit ``chi = N alpha / eps0`` is the linear
density approximation; it differs from Lorentz-Lorenz by the factor
``1 / (1 - y) ~ 1 + chi/3``, a relative error of about 2e-4 for argon at 1 bar but about
2 % at 100 bar.

Beyond Lorentz-Lorenz, the refractivity virial expansion ``(n^2 - 1)/(n^2 + 2) = A_R rho
+ B_R rho^2 + ...`` (``rho`` the molar density, ``A_R`` the molar refractivity) adds the
pair term. Every function here takes an optional ratio ``virial = B_R / A_R``
(m^3/mol) and uses ``y = (N alpha / 3 eps0)(1 + virial rho)``; the ratio is of order
1e-7 to 1e-6 m^3/mol, a correction of 1e-6 to 1e-5 at 1 atm.

Atomic-unit polarizabilities convert with
:data:`~optmat.constants.AU_POLARIZABILITY`; the molar refractivity
``A_R = N_A alpha / (3 eps0)``, used by gas metrology, with :func:`molar_refractivity`.
"""

from typing import Final

import numpy as np
from jaxtyping import Complex, Float
from numpy.typing import ArrayLike

from optmat.constants import EPS0, N_A

__all__ = [
    "EVALUATION_VERSION",
    "dpolarizability_dsusceptibility",
    "molar_refractivity",
    "polarizability_from_refractivity",
    "polarizability_from_susceptibility",
    "refractive_index",
    "refractivity",
    "susceptibility",
    "susceptibility_dilute",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""


def _pair_factor(density: float, virial: float) -> float:
    """Return ``1 + b rho``, the refractivity virial correction at a density.

    Parameters
    ----------
    density : float
        Number density, m^-3.
    virial : float
        ``B_R / A_R``, m^3/mol.

    Returns
    -------
    float
        The factor multiplying ``N alpha / (3 eps0)``.
    """
    return 1.0 + virial * density / N_A


def susceptibility(
    polarizability: ArrayLike | float, density: float, virial: float = 0.0
) -> Float[np.ndarray, "*shape"]:
    """Return ``chi = n^2 - 1`` from the Lorentz-Lorenz relation.

    Parameters
    ----------
    polarizability : array_like or float
        Per-particle dipole polarizability ``alpha`` in C m^2/V (real).
    density : float
        Number density ``N`` in m^-3, >= 0.
    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        ``chi = 3 y / (1 - y)`` with ``y = (N alpha / 3 eps0)(1 + b rho)``, the shape of
        ``polarizability``.

    Examples
    --------
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> f"{float(susceptibility(11.08 * AU_POLARIZABILITY, 2.687e25)):.6e}"
    '5.544992e-04'
    """
    y = density * np.asarray(polarizability, dtype=float) / (3.0 * EPS0)
    y = y * _pair_factor(density, virial)
    return np.asarray(3.0 * y / (1.0 - y))


def susceptibility_dilute(
    polarizability: ArrayLike | float, density: float
) -> Float[np.ndarray, "*shape"]:
    """Return the dilute-limit susceptibility ``chi = N alpha / eps0``.

    This is linear density scaling, exact only as ``N -> 0``; it is provided to measure
    the error of that approximation against :func:`susceptibility`.

    Parameters
    ----------
    polarizability : array_like or float
        Per-particle dipole polarizability in C m^2/V.
    density : float
        Number density in m^-3.

    Returns
    -------
    numpy.ndarray
        ``N alpha / eps0``, the shape of ``polarizability``.

    Examples
    --------
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> f"{float(susceptibility_dilute(11.08 * AU_POLARIZABILITY, 2.687e25)):.6e}"
    '5.543968e-04'
    """
    return np.asarray(density * np.asarray(polarizability, dtype=float) / EPS0)


def polarizability_from_susceptibility(
    chi: ArrayLike | float, density: float, virial: float = 0.0
) -> Float[np.ndarray, "*shape"]:
    """Return the polarizability from ``chi = n^2 - 1`` (inverse Lorentz-Lorenz).

    Parameters
    ----------
    chi : array_like or float
        Susceptibility ``n^2 - 1``, > -1.
    density : float
        Number density in m^-3, > 0.

    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        ``alpha = 3 eps0 y / N`` with ``y = chi / (chi + 3)``, in C m^2/V.

    Examples
    --------
    >>> alpha = polarizability_from_susceptibility(5.544992e-4, 2.687e25)
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> f"{float(alpha / AU_POLARIZABILITY):.4f}"
    '11.0800'
    """
    chi = np.asarray(chi, dtype=float)
    alpha = 3.0 * EPS0 * (chi / (chi + 3.0)) / density
    return np.asarray(alpha / _pair_factor(density, virial))


def dpolarizability_dsusceptibility(
    chi: ArrayLike | float, density: float, virial: float = 0.0
) -> Float[np.ndarray, "*shape"]:
    """Return d alpha / d chi of the inverse Lorentz-Lorenz relation.

    With ``alpha = 3 eps0 chi / [N P (chi + 3)]``
    (:func:`polarizability_from_susceptibility`, ``P`` the pair factor of the
    refractivity virial at ``N``), ``d alpha / d chi = 9 eps0 / [N P (chi + 3)^2]``.
    For a formula in ``n - 1``, ``d chi / d(n - 1) = 2 n``.

    Parameters
    ----------
    chi : array_like or float
        Susceptibility ``n^2 - 1``, > -1; complex for an absorbing medium.
    density : float
        Number density in m^-3, > 0.
    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        ``d alpha / d chi``, C m^2/V.

    Examples
    --------
    >>> h = 1e-9
    >>> numeric = (polarizability_from_susceptibility(5e-4 + h, 2.687e25)
    ...            - polarizability_from_susceptibility(5e-4 - h, 2.687e25)) / (2 * h)
    >>> bool(abs(dpolarizability_dsusceptibility(5e-4, 2.687e25) / numeric - 1) < 1e-6)
    True
    """
    chi = np.asarray(chi)
    return np.asarray(
        9.0 * EPS0 / (density * (chi + 3.0) ** 2) / _pair_factor(density, virial)
    )


def refractivity(
    polarizability: ArrayLike | float, density: float, virial: float = 0.0
) -> Float[np.ndarray, "*shape"]:
    """Return ``n - 1`` from the Lorentz-Lorenz relation, without cancellation.

    Parameters
    ----------
    polarizability : array_like or float
        Per-particle dipole polarizability in C m^2/V.
    density : float
        Number density in m^-3.

    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        ``n - 1 = chi / (1 + sqrt(1 + chi))``, the shape of ``polarizability``.

    Examples
    --------
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> f"{float(refractivity(11.08 * AU_POLARIZABILITY, 2.687e25)):.6e}"
    '2.772112e-04'
    """
    chi = susceptibility(polarizability, density, virial)
    return np.asarray(chi / (1.0 + np.sqrt(1.0 + chi)))


def polarizability_from_refractivity(
    refractivity_value: ArrayLike | float, density: float, virial: float = 0.0
) -> Float[np.ndarray, "*shape"]:
    """Return the polarizability from ``n - 1`` at number density ``N``.

    Parameters
    ----------
    refractivity_value : array_like or float
        ``n - 1``.
    density : float
        Number density in m^-3, > 0.

    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        The polarizability in C m^2/V, the shape of ``refractivity_value``.

    Examples
    --------
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> alpha = polarizability_from_refractivity(2.772112e-4, 2.687e25)
    >>> f"{float(alpha / AU_POLARIZABILITY):.5f}"
    '11.08000'
    """
    dn = np.asarray(refractivity_value, dtype=float)
    return polarizability_from_susceptibility(dn * (dn + 2.0), density, virial)


def molar_refractivity(
    polarizability: ArrayLike | float,
) -> Float[np.ndarray, "*shape"]:
    """Return the molar refractivity ``A_R = N_A alpha / (3 eps0)``.

    Parameters
    ----------
    polarizability : array_like or float
        Per-particle dipole polarizability in C m^2/V.

    Returns
    -------
    numpy.ndarray
        Molar refractivity in m^3/mol (multiply by 1e6 for the cm^3/mol of gas
        metrology).

    Examples
    --------
    Helium's static polarizability, 1.383760 a.u., gives the static molar
    refractivity of gas metrology:

    >>> from optmat.constants import AU_POLARIZABILITY
    >>> f"{float(molar_refractivity(1.383760 * AU_POLARIZABILITY)) * 1e6:.5f}"
    '0.51725'
    """
    return np.asarray(N_A * np.asarray(polarizability, dtype=float) / (3.0 * EPS0))


def refractive_index(
    polarizability: ArrayLike | complex, density: float, virial: float = 0.0
) -> Complex[np.ndarray, "*shape"]:
    """Return the complex refractive index from a complex polarizability.

    Parameters
    ----------
    polarizability : array_like or complex
        Per-particle dipole polarizability in C m^2/V; ``Im alpha > 0`` is absorption.
    density : float
        Number density in m^-3.
    virial : float, optional
        ``B_R / A_R`` in m^3/mol (0 by default: no pair term).

    Returns
    -------
    numpy.ndarray
        ``n = sqrt(1 + chi)`` with ``chi = 3 y / (1 - y)`` and ``y = (N alpha / 3
        eps0)(1 + b rho)`` (complex; ``Im n > 0`` is absorption), the shape of
        ``polarizability``.

    Examples
    --------
    >>> from optmat.constants import AU_POLARIZABILITY
    >>> n = refractive_index((11.08 + 0.01j) * AU_POLARIZABILITY, 2.687e25)
    >>> f"{float(n.real - 1):.6e} {float(n.imag):.3e}"
    '2.772112e-04 2.502e-07'
    """
    y = density * np.asarray(polarizability, dtype=complex) / (3.0 * EPS0)
    y = y * _pair_factor(density, virial)
    return np.asarray(np.sqrt(1.0 + 3.0 * y / (1.0 - y)))

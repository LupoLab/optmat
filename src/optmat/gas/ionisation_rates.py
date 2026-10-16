# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Strong-field ionisation rates of the noble gases (TC-652 to TC-662).

The rate at which an atom loses an electron in a linearly polarised field of amplitude
:math:`F` and angular frequency :math:`\omega` is built here from three factors, each
validated against its own references (the manual page *Ionisation rates*):

.. math::

   W(F, \omega) = N_\text{eff} \sum_c p_c\;
       w_c^\text{PPT}(F, \omega)\; S_c(F).

- :math:`w^\text{PPT}` is the rate of **one electron** in the orbital
  :math:`(l, m = 0)`: Perelomov, Popov and Terent'ev's rate in a short-range potential
  (their eq. 54) times a Coulomb factor (:func:`coulomb_factor`). The default's is
  Perelomov and Popov's (1967) static factor :math:`(2F_0/F)^{2n^*}`, eq. (42), times
  :math:`(1 + 2\gamma/e)^{-2n^* c s}` with :math:`s = n^*\omega/2I_p`, fitted to
  one-electron calculations (``frequency``; TC-831, TC-848). For :math:`\omega \to 0`
  it is the cycle average of the static tunnelling rate (:func:`tunnelling_rate_au`).
- :math:`S` is the correction of that rate at finite field, fitted by Tong and Lin
  (2005) or by Zhang, Lan and Lu (2014) to exact static rates of one electron in a model
  potential (:func:`tong_lin_factor`, :func:`zhang_factor`).
- :math:`N_\text{eff}` is the rate of the atom over that of the one electron, a
  constant per gas set from all-electron calculations and measurements (``data/gas/
  ionisation_rates.toml``). Its weak-field limit is 2, the two electrons of the
  :math:`m = 0` orbital (Tolstikhin, Madsen and Morishita 2014, section II G); at
  working fields it is 0.8 to 1.5: about 1 for argon, 1.35 for krypton, 1.5 for
  xenon, 1.47 for helium (TC-847, TC-848, TC-850).
- The sum is over the fine-structure levels :math:`c` of the ion a :math:`p` shell
  leaves, :math:`^2P_{3/2}` and :math:`^2P_{1/2}`, each with its own ionisation
  energy, with the weights :math:`p_c = 2/3` and :math:`1/3` of the :math:`m = 0`
  orbital.

Atomic units throughout the functions whose names end in ``_au``: :math:`\kappa =
\sqrt{2 I_p}`, :math:`n^* = Z/\kappa`, :math:`F_0 = \kappa^3`, and the Keldysh
parameter :math:`\gamma = \omega\kappa/F`. :class:`IonisationRateModel` takes and
returns SI: fields in V/m, rates in s\ :sup:`-1`, wavelengths in m.

**Two rates** are given, of the same model:

- the **cycle-averaged** rate as a function of the field's amplitude, which is what
  the theory defines (:meth:`IonisationRateModel.cycle_averaged_rate`);
- the **instantaneous** rate as a function of the field itself, for a calculation
  that resolves the carrier (:meth:`IonisationRateModel.instantaneous_rate`). Its
  average over a cycle is the cycle-averaged rate, to 4 %, at any Keldysh parameter
  (TC-662, TC-828), so that a pulse ionises the same fraction
  whichever rate is used; in a static field it is the static rate. It does not give
  the timing of the ionisation within the cycle at :math:`\gamma \gtrsim 1`.

**What is not here**, and why (the manual page has the evidence):

- the factor :math:`(1 + \gamma^2)^{-n^*}` that Ilkov, Decker and Chin (1992, eq. 8)
  and Bergé et al. (2007, eq. B14) carry in the Coulomb term. No paper derives it, and
  it makes the rate too low by 2 at :math:`\gamma = 1` and by 8 at :math:`\gamma = 3`
  against time-dependent Schrödinger calculations. It is available
  (``coulomb="ilkov"``), as is Popruzhenko et al.'s factor (``"popruzhenko"``), to
  reproduce rates computed with them;
- a Stark shift of the ionisation energy and a factor for the polarisation of the ion.
  For the active electron they are inside :math:`S`; for the rest of the atom they are
  inside :math:`N_\text{eff}`.

References
----------
- A. M. Perelomov, V. S. Popov and M. V. Terent'ev, "Ionization of atoms in an
  alternating electric field", Sov. Phys. JETP 23, 924 (1966).
- A. M. Perelomov and V. S. Popov, "Ionization of atoms in an alternating
  electrical field. III", Sov. Phys. JETP 25, 336 (1967).
- A. M. Perelomov, V. S. Popov and V. P. Kuznetsov, "Allowance for the Coulomb
  interaction in multiphoton ionization", Sov. Phys. JETP 27, 451 (1968).
- S. F. J. Larochelle, A. Talebpour and S. L. Chin, "Coulomb effect in
  multiphoton ionization of rare-gas atoms", J. Phys. B 31, 1215 (1998).
- X. M. Tong, Z. X. Zhao and C. D. Lin, "Theory of molecular tunneling ionization",
  Phys. Rev. A 66, 033402 (2002).
- X. M. Tong and C. D. Lin, "Empirical formula for static field ionization rates
  of atoms and molecules by lasers in the barrier-suppression regime", J. Phys. B 38,
  2593 (2005).
- Q. Zhang, P. Lan and P. Lu, "Empirical formula for over-barrier strong-field
  ionization", Phys. Rev. A 90, 043410 (2014).
- F. A. Ilkov, J. E. Decker and S. L. Chin, "Ionization of atoms in the tunnelling
  regime with experimental evidence using Hg atoms", J. Phys. B 25, 4005 (1992).
- L. Bergé, S. Skupin, R. Nuter, J. Kasparian and J.-P. Wolf, "Ultrashort
  filaments of light in weakly ionized, optically transparent media", Rep. Prog. Phys.
  70, 1633 (2007).
- S. V. Popruzhenko, V. D. Mur, V. S. Popov and D. Bauer, "Strong field ionization
  rate for arbitrary laser frequencies", Phys. Rev. Lett. 101, 193003 (2008).
- O. I. Tolstikhin, L. B. Madsen and T. Morishita, "Weak-field asymptotic theory
  of tunneling ionization in many-electron atomic and molecular systems", Phys. Rev. A
  89, 013421 (2014).
- M. V. Ammosov, N. B. Delone and V. P. Krainov, "Tunnel ionization of complex
  atoms and of atomic ions in an alternating electromagnetic field", Sov. Phys. JETP
  64, 1191 (1986).
"""

import functools
import math
import tomllib
import warnings
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from typing import Any, ClassVar, Final, Literal, get_args

import numpy as np
from jaxtyping import Float
from numpy.typing import ArrayLike
from scipy import special

from optmat import _catalogue
from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_ENERGY,
    AU_FIELD,
    E_CHARGE,
    HBAR,
    C,
)
from optmat.datasets import Dataset
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.gas.ionisation import (
    IonisationParameters,
    hartree_coefficient_squared,
    ionisation_parameters,
)
from optmat.identity import Identified, parse_status, source_ids

__all__ = [
    "EVALUATION_VERSION",
    "IONISATION_RATES_FILE",
    "Averaging",
    "Channel",
    "ChannelSum",
    "CoulombFactor",
    "IonisationRateModel",
    "StaticCorrection",
    "Structure",
    "channel_factor",
    "channel_function",
    "circular_enhancement",
    "circular_saddle_point",
    "coulomb_factor",
    "cycle_average",
    "cycle_average_factor",
    "ionisation_rate_model",
    "ionisation_rate_models",
    "keldysh_parameter",
    "power_law_average",
    "ppt_exponent_function",
    "short_range_rate_au",
    "tong_lin_factor",
    "tunnelling_rate_au",
    "zhang_factor",
]

EVALUATION_VERSION: Final = 2
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

IONISATION_RATES_FILE: Final = GAS_DATA_DIRECTORY / "ionisation_rates.toml"
"""The named rate models and their constants per gas, each with its reason."""

CoulombFactor = Literal["static", "ilkov", "popruzhenko", "frequency"]
"""The Coulomb factor of the rate (:func:`coulomb_factor`)."""

StaticCorrection = Literal["tong-lin", "zhang", "none"]
"""The correction of the one-electron static rate at finite field."""

Averaging = Literal["consistent", "tunnelling"]
"""How the cycle-averaged and the instantaneous rate are related.

``consistent``: the average of the instantaneous rate over a cycle is the
cycle-averaged rate, at any Keldysh parameter (the averages over a cycle computed by
quadrature; TC-662, TC-828).
``tunnelling``: the two differ by the factor of the tunnelling limit,
:func:`cycle_average_factor`, as in the literature's forms.
"""

ChannelSum = Literal["integral", "sum"]
"""How the sum over the above-threshold channels is taken (:func:`channel_factor`)."""

Structure = Literal["tong_2002", "hartree"]
"""Where the asymptotic coefficient of the orbital comes from."""

_SERIES_LIMIT: Final = 1e-2
"""Keldysh parameter below which the functions of gamma come from their series.

The closed forms are differences of terms that agree to leading order, and lose
``eps / gamma^2`` of relative accuracy. The series are truncated at relative order
``gamma^6``, 1e-12 at the switch, so the two are equally accurate where they meet.
"""

_SLOPE_STEP: Final = 1e-3
"""Relative step of the field in the logarithmic slope that scales the average's nodes.

The slope only places the nodes of :func:`_log_cycle_averages`; its rounding error,
1e-16 x 30 / 1e-3 = 3e-12 for a logarithm of 30, moves the average by far less.
"""

_AVERAGE_NODES: Final = 32
"""Gauss-Legendre nodes of the average over a cycle (:func:`_log_cycle_averages`)."""

_AVERAGE_SPAN: Final = 5.0
"""The nodes span t sqrt(1 + s) in [0, 5], where the integrand has fallen to e^-25."""

_FREQUENCY_CONSTANT: Final = 5.25
"""The constant :math:`c` of the ``frequency`` Coulomb factor (:func:`coulomb_factor`).

Fitted to the one-electron ionised fractions of argon for one 10 fs pulse at 400, 800
and 1030 nm (``travers_2026_trecx``, yields below 0.75): the least root-mean-square
logarithm of model over calculation, a factor of 1.17 against 1.78 for the static
factor (TC-831).
"""

_SUM_BLOCK: Final = 256
"""Terms of the channel sum evaluated per block."""

_SUM_MAXIMUM: Final = 1 << 22
"""Most terms of the channel sum before it is refused."""

_ZHANG_RANGE: Final = (0.5, 4.5)
"""The range of F/F_b over which Zhang et al. fitted their factor."""

_TONG_LIN_MAXIMUM: Final = 2.0
"""F/F_b up to which Tong and Lin's factor holds (Zhang et al. 2014)."""

_MINIMUM_WAVELENGTH: Final = 350e-9
"""Shortest wavelength, m, at which the rates were compared with references (390 nm).

In the vacuum ultraviolet the rate fails: for argon at 140 nm it is an order of
magnitude above a time-dependent Schrödinger calculation (TC-653).
"""


def keldysh_parameter(
    field_au: Float[np.ndarray, " *shape"], omega_au: float, ionization_energy_au: float
) -> Float[np.ndarray, " *shape"]:
    r"""Return the Keldysh parameter :math:`\gamma = \omega\kappa / F`.

    Perelomov et al. (1966), after eq. (1): the ratio of the field's frequency to the
    frequency of tunnelling, :math:`\omega_t = F/\kappa`.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency of the field, atomic units, >= 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.

    Returns
    -------
    numpy.ndarray
        :math:`\gamma`, with the shape of ``field_au``.

    Examples
    --------
    Argon at 800 nm and 1e14 W/cm\ :sup:`2` (0.05338 au):

    >>> import numpy as np
    >>> gamma = keldysh_parameter(np.array(0.05338), 0.056954, 0.57916)
    >>> f"{float(gamma):.3f}"
    '1.148'
    """
    return np.asarray(omega_au * math.sqrt(2 * ionization_energy_au) / field_au)


def ppt_exponent_function(
    gamma: Float[np.ndarray, " *shape"],
) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`g(\gamma)` of the exponent of the rate.

    Perelomov et al. (1966), eq. (33):

    .. math::

       g(\gamma) = \frac{3}{2\gamma} \left[ \left(1 + \frac{1}{2\gamma^2}\right)
           \operatorname{arsinh}\gamma - \frac{\sqrt{1 + \gamma^2}}{2\gamma} \right],

    so that the rate is proportional to :math:`\exp[-(2F_0/3F)\, g(\gamma)]`. It is 1
    in the tunnelling limit, where it has the series :math:`1 - \gamma^2/10 +
    9\gamma^4/280` (their eq. 59 gives the first two terms).

    Parameters
    ----------
    gamma : numpy.ndarray
        Keldysh parameter, >= 0.

    Returns
    -------
    numpy.ndarray
        :math:`g(\gamma)`, in (0, 1].

    Examples
    --------
    >>> import numpy as np
    >>> f"{float(ppt_exponent_function(np.array(1.0))):.6f}"
    '0.922430'
    """
    gamma = np.asarray(gamma, dtype=float)
    small = gamma < _SERIES_LIMIT
    safe = np.where(small, 1.0, gamma)
    closed = (
        1.5
        / safe
        * ((1 + 0.5 / safe**2) * np.arcsinh(safe) - np.sqrt(1 + safe**2) / (2 * safe))
    )
    series = 1 - gamma**2 / 10 + 9 * gamma**4 / 280
    return np.asarray(np.where(small, series, closed))


def _alpha(gamma: Float[np.ndarray, " *shape"]) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`\alpha(\gamma)`, the decay constant of the channel sum.

    Perelomov et al. (1966), eq. (34): :math:`\alpha = 2[\operatorname{arsinh}\gamma
    - \gamma/\sqrt{1+\gamma^2}]`. Its
    series is :math:`(2/3)\gamma^3 (1 - 9\gamma^2/10 + 45\gamma^4/56)`.

    Parameters
    ----------
    gamma : numpy.ndarray
        Keldysh parameter, > 0.

    Returns
    -------
    numpy.ndarray
        :math:`\alpha(\gamma)`, > 0.
    """
    small = gamma < _SERIES_LIMIT
    closed = 2 * (np.arcsinh(gamma) - gamma / np.sqrt(1 + gamma**2))
    series = 2 / 3 * gamma**3 * (1 - 0.9 * gamma**2 + 45 / 56 * gamma**4)
    return np.asarray(np.where(small, series, closed))


def _beta(gamma: Float[np.ndarray, " *shape"]) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`\beta(\gamma) = 2\gamma/\sqrt{1 + \gamma^2}`.

    Perelomov et al. (1966), after eq. (55).

    Parameters
    ----------
    gamma : numpy.ndarray
        Keldysh parameter, > 0.

    Returns
    -------
    numpy.ndarray
        :math:`\beta(\gamma)`, in (0, 2).
    """
    return np.asarray(2 * gamma / np.sqrt(1 + gamma**2))


def channel_function(
    m: int, x: Float[np.ndarray, " *shape"]
) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`w_m(x)`, the weight of one above-threshold channel.

    Perelomov et al. (1966), eq. (56):

    .. math::

       w_m(x) = e^{-x^2} \int_0^x e^{y^2} (x^2 - y^2)^{|m|}\, dy.

    For :math:`m = 0` it is Dawson's integral. For any :math:`m` the integral is
    :math:`\tfrac{\sqrt\pi}{2} \tfrac{|m|!}{\Gamma(|m| + 3/2)} x^{2|m|+1}
    \,{}_1F_1(\tfrac12; |m| + \tfrac32; x^2)`, and Kummer's transformation,
    :math:`e^{-z}\,{}_1F_1(a; b; z) = {}_1F_1(b - a; b; -z)`, removes the exponential,
    so that nothing overflows at large :math:`x`.

    Parameters
    ----------
    m : int
        Magnetic quantum number; only :math:`|m|` enters.
    x : numpy.ndarray
        Argument, >= 0.

    Returns
    -------
    numpy.ndarray
        :math:`w_m(x)`, which rises as :math:`x^{2|m|+1}` and falls as
        :math:`|m|!/(2x)`.

    Examples
    --------
    >>> import numpy as np
    >>> f"{float(channel_function(0, np.array(1.0))):.6f}"
    '0.538080'
    """
    order = abs(m)
    x = np.asarray(x, dtype=float)
    if order == 0:
        return np.asarray(special.dawsn(x))
    prefactor = (
        math.sqrt(math.pi) * math.factorial(order) / (2 * math.gamma(order + 1.5))
    )
    return np.asarray(
        prefactor
        * x ** (2 * order + 1)
        * special.hyp1f1(order + 1.0, order + 1.5, -(x**2))
    )


def _channel_sum(
    m: int,
    alpha: Float[np.ndarray, " n"],
    beta: Float[np.ndarray, " n"],
    nu: Float[np.ndarray, " n"],
    tolerance: float,
) -> Float[np.ndarray, " n"]:
    r"""Return :math:`\sum_{n \ge \nu} e^{-\alpha(n - \nu)} w_m(\sqrt{\beta(n - \nu)})`.

    The sum of Perelomov et al.'s eq. (55), over the numbers of photons above the
    threshold :math:`\nu`. Past the maximum of :math:`w_m` the terms fall at least as
    fast as :math:`e^{-\alpha k}`, so the remainder after a term :math:`t` is below
    :math:`t/\alpha`; blocks of terms are added until that bound is below
    ``tolerance`` of the sum.

    Parameters
    ----------
    m : int
        Magnetic quantum number.
    alpha, beta : numpy.ndarray
        The functions of gamma, shape ``(n,)``.
    nu : numpy.ndarray
        Threshold in photons, shape ``(n,)``.
    tolerance : float
        Relative bound on the truncated remainder.

    Returns
    -------
    numpy.ndarray
        The sum, shape ``(n,)``.

    Raises
    ------
    DataError
        If the sum needs more than :data:`_SUM_MAXIMUM` terms (a Keldysh parameter
        below about 0.02; use the integral there).
    """
    offset = np.ceil(nu) - nu
    total = np.zeros_like(alpha)
    # w_m peaks below x = 2 (|m| + 1): the bound on the remainder holds past it.
    minimum = int(np.max(4 * (abs(m) + 1) ** 2 / beta)) + 2
    start = 0
    while start < _SUM_MAXIMUM:
        index = np.arange(start, start + _SUM_BLOCK, dtype=float)
        distance = offset[:, None] + index[None, :]
        terms = np.exp(-alpha[:, None] * distance) * channel_function(
            m, np.sqrt(beta[:, None] * distance)
        )
        total += terms.sum(axis=1)
        start += _SUM_BLOCK
        if start >= minimum and np.all(terms[:, -1] <= tolerance * alpha * total):
            return total
    msg = (
        f"the channel sum needs more than {_SUM_MAXIMUM} terms (smallest alpha "
        f"{float(np.min(alpha)):.2g}); use channel_sum='integral'"
    )
    raise DataError(msg)


def channel_factor(
    m: int,
    gamma: Float[np.ndarray, " *shape"],
    threshold: Float[np.ndarray, " *shape"] | None = None,
    *,
    tolerance: float = 1e-6,
) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`A_m(\omega, \gamma)`, the factor of the above-threshold channels.

    Perelomov et al. (1966), eq. (55):

    .. math::

       A_m(\omega, \gamma) = \frac{4}{\sqrt{3\pi}} \frac{1}{|m|!}
           \frac{\gamma^2}{1 + \gamma^2}
           \sum_{n \ge \nu} e^{-\alpha (n - \nu)}\, w_m\!\left(\sqrt{\beta (n - \nu)}
           \right),
       \qquad \nu = \frac{I_p}{\omega} \left(1 + \frac{1}{2\gamma^2}\right),

    which tends to 1 for :math:`\gamma \to 0`. With ``threshold`` (:math:`\nu`) the
    sum is taken term by term. It then has a step wherever :math:`\nu` passes an
    integer (a channel closes as the ponderomotive energy grows). Without it the sum is
    replaced by its integral over the number of photons, which is smooth:

    .. math::

       \int_0^\infty e^{-\alpha x}\, w_m(\sqrt{\beta x})\, dx
       = \frac{\sqrt\pi}{2}\, |m|!\, \frac{\beta^{|m|}}{(\alpha + \beta)^{|m| + 1}}
         \sqrt{\frac{\beta}{\alpha}}.

    The integral is within 10 % of the sum for :math:`\gamma < 1.5` and within 50 % up
    to :math:`\gamma = 3` (TC-654).

    Parameters
    ----------
    m : int
        Magnetic quantum number.
    gamma : numpy.ndarray
        Keldysh parameter, > 0.
    threshold : numpy.ndarray or None, optional
        :math:`\nu`, the threshold in photons, with the shape of ``gamma``: take the
        sum. None: take the integral.
    tolerance : float, optional
        Relative bound on the truncated remainder of the sum.

    Returns
    -------
    numpy.ndarray
        :math:`A_m`, with the shape of ``gamma``.

    Examples
    --------
    The tunnelling limit of the integral:

    >>> import numpy as np
    >>> f"{float(channel_factor(0, np.array(1e-3))):.4f}"
    '1.0000'
    """
    gamma = np.asarray(gamma, dtype=float)
    order = abs(m)
    alpha = _alpha(gamma)
    beta = _beta(gamma)
    if threshold is None:
        total = (
            0.5
            * math.sqrt(math.pi)
            * math.factorial(order)
            * beta**order
            / (alpha + beta) ** (order + 1)
            * np.sqrt(beta / alpha)
        )
    else:
        shape = gamma.shape
        total = _channel_sum(
            order,
            alpha.ravel(),
            beta.ravel(),
            np.broadcast_to(np.asarray(threshold, dtype=float), shape).ravel(),
            tolerance,
        ).reshape(shape)
    return np.asarray(
        4
        / math.sqrt(3 * math.pi)
        / math.factorial(order)
        * gamma**2
        / (1 + gamma**2)
        * total
    )


def _angular_factor(l: int, m: int) -> float:  # noqa: E741 - l is the orbital's symbol
    r"""Return :math:`f(l, m) = (2l+1)(l+|m|)! / (2^{|m|} |m|! (l-|m|)!)`.

    Perelomov et al. (1966), eq. (4).

    Parameters
    ----------
    l : int
        Orbital angular momentum, >= 0.
    m : int
        Magnetic quantum number, ``|m| <= l``.

    Returns
    -------
    float
        The angular factor.

    Raises
    ------
    ValueError
        If ``|m| > l``.
    """
    order = abs(m)
    if order > l:
        msg = f"|m| = {order} exceeds l = {l}"
        raise ValueError(msg)
    return (
        (2 * l + 1)
        * math.factorial(l + order)
        / (2**order * math.factorial(order) * math.factorial(l - order))
    )


def tunnelling_rate_au(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    coefficient_squared: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
    m: int = 0,
    charge: float = 1.0,
) -> Float[np.ndarray, " *shape"]:
    r"""Return the static tunnelling rate of one electron, to leading order.

    Perelomov et al. (1966), eq. (4):

    .. math::

       w_\text{stat}(F) = I_p\, |C_{\kappa l}|^2\, f(l, m)
           \left(\frac{2F_0}{F}\right)^{2n^* - |m| - 1}
           \exp\left(-\frac{2F_0}{3F}\right),

    with the orbital's asymptotic form :math:`\psi \to C_{\kappa l}\, \kappa^{3/2}
    (\kappa r)^{n^* - 1} e^{-\kappa r}\, Y_{lm}` (their eq. 6). Tong, Zhao and Lin's
    (2002) eq. (7) is the same rate with their coefficient :math:`C_l =
    C_{\kappa l}\, \kappa^{n^* + 1/2}`. It is exact for :math:`F \to 0` and too high at
    any field of practical interest: by 1.9 for argon at 0.64 of the barrier-suppression
    field (TC-655).

    Parameters
    ----------
    field_au : numpy.ndarray
        Static field, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.
    coefficient_squared : float
        :math:`|C_{\kappa l}|^2`, in Perelomov et al.'s normalisation.
    l : int
        Orbital angular momentum of the electron.
    m : int, optional
        Its magnetic quantum number along the field.
    charge : float, optional
        Charge :math:`Z` of the ion left behind.

    Returns
    -------
    numpy.ndarray
        Rate, atomic units, with the shape of ``field_au``.

    Examples
    --------
    Hydrogen (:math:`C = 2`) gives :math:`(4/F) e^{-2/3F}`:

    >>> import numpy as np
    >>> w = tunnelling_rate_au(np.array(0.05), 0.5, 4.0, 0)
    >>> f"{float(w) / (4 / 0.05 * np.exp(-2 / 0.15)):.12f}"
    '1.000000000000'
    """
    return np.asarray(
        np.exp(
            _log_tunnelling_rate(
                np.asarray(field_au, dtype=float),
                ionization_energy_au,
                coefficient_squared,
                l,
                m,
                charge,
            )
        )
    )


def _log_tunnelling_rate(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    coefficient_squared: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
    m: int,
    charge: float,
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of :func:`tunnelling_rate_au`.

    Parameters
    ----------
    field_au : numpy.ndarray
        Static field, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    coefficient_squared : float
        ``|C_kappa l|^2``.
    l : int
        Orbital angular momentum.
    m : int
        Magnetic quantum number.
    charge : float
        Charge of the ion.

    Returns
    -------
    numpy.ndarray
        The logarithm of the rate in atomic units.
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    n_star = charge / kappa
    ratio = 2 * kappa**3 / field_au
    return np.asarray(
        math.log(ionization_energy_au * coefficient_squared * _angular_factor(l, m))
        + (2 * n_star - abs(m) - 1) * np.log(ratio)
        - ratio / 3
    )


def cycle_average_factor(
    field_au: Float[np.ndarray, " *shape"], ionization_energy_au: float
) -> Float[np.ndarray, " *shape"]:
    r"""Return :math:`(3F/\pi F_0)^{1/2}`, the cycle average in the tunnelling limit.

    Perelomov et al. (1966), eq. (7): for :math:`\omega \ll \omega_t` the rate in a
    field of amplitude :math:`F` is the static rate at :math:`F` times this factor,
    the average of :math:`w_\text{stat}(F |\cos\omega t|)` over a cycle.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.

    Returns
    -------
    numpy.ndarray
        The factor, with the shape of ``field_au``.

    Examples
    --------
    >>> import numpy as np
    >>> f"{float(cycle_average_factor(np.array(0.05), 0.5)):.6f}"
    '0.218510'
    """
    f0 = (2 * ionization_energy_au) ** 1.5
    return np.asarray(np.sqrt(3 * np.asarray(field_au, dtype=float) / (math.pi * f0)))


def power_law_average(
    slope: Float[np.ndarray, " *shape"],
) -> Float[np.ndarray, " *shape"]:
    r"""Return the average over a cycle of a rate that goes as a power of the field.

    For :math:`W \propto F^s`, the average of :math:`W(F|\cos\omega t|)` over a
    cycle is :math:`c(s)\, W(F)` with

    .. math::

       c(s) = \frac{2}{\pi} \int_0^{\pi/2} \cos^s\varphi\, d\varphi
            = \frac{\Gamma\!\left(\frac{s + 1}{2}\right)}
                   {\sqrt\pi\, \Gamma\!\left(\frac{s}{2} + 1\right)}
            \approx \sqrt{\frac{2}{\pi s}} .

    A rate that is steep in the field is a power law over the part of the cycle
    that contributes, with the local slope :math:`s = d\ln W/d\ln F`; for the
    tunnelling rate :math:`s \to 2F_0/3F` and :math:`c \to (3F/\pi F_0)^{1/2}`,
    Perelomov et al.'s (1966) eq. (7). :func:`cycle_average` adds the curvature.

    Parameters
    ----------
    slope : numpy.ndarray
        Logarithmic slope of the rate, > -1.

    Returns
    -------
    numpy.ndarray
        :math:`c(s)`, in (0, 1].

    Examples
    --------
    >>> import numpy as np
    >>> f"{float(power_law_average(np.array(2.0))):.4f}"
    '0.5000'
    """
    slope = np.asarray(slope, dtype=float)
    return np.asarray(
        np.exp(special.gammaln((slope + 1) / 2) - special.gammaln(slope / 2 + 1))
        / math.sqrt(math.pi)
    )


def cycle_average(
    slope: Float[np.ndarray, " *shape"], curvature: Float[np.ndarray, " *shape"]
) -> Float[np.ndarray, " *shape"]:
    r"""Return the average over a cycle of a rate, over the rate at the amplitude.

    With :math:`L(\ln F) = \ln W`, its slope :math:`s = L'` and its curvature
    :math:`q = L''`, and :math:`u = \ln|\cos\varphi|`,

    .. math::

       \frac{\overline{W}}{W(F)}
       = \frac{2}{\pi} \int_0^{\pi/2} e^{s u + q u^2/2 + \dots}\, d\varphi
       \approx c(s) \left[1 + \frac{q}{2}\, \frac{c''(s)}{c(s)}\right],

    with :math:`c` of :func:`power_law_average`, whose derivatives with respect to
    :math:`s` are the moments of :math:`u`. For the tunnelling rate :math:`q = -s`
    and the bracket is :math:`1 - 3/8s`. The average of the default model's
    instantaneous rate over a cycle is its cycle-averaged rate to 1 % up to the
    barrier-suppression field (TC-662). The models no longer use it: where the rate
    bends sharply the expansion fails, 24 % low for argon's average at twice the
    barrier-suppression field, and they average by quadrature instead (TC-828).

    Parameters
    ----------
    slope : numpy.ndarray
        :math:`d\ln W/d\ln F`, > 0.
    curvature : numpy.ndarray
        :math:`d^2\ln W/d(\ln F)^2`.

    Returns
    -------
    numpy.ndarray
        The ratio, in (0, 1].

    Examples
    --------
    A power law has no curvature:

    >>> import numpy as np
    >>> f"{float(cycle_average(np.array(2.0), np.array(0.0))):.4f}"
    '0.5000'
    """
    slope = np.asarray(slope, dtype=float)
    first = 0.5 * (special.digamma((slope + 1) / 2) - special.digamma(slope / 2 + 1))
    second = 0.25 * (
        special.polygamma(1, (slope + 1) / 2) - special.polygamma(1, slope / 2 + 1)
    )
    # c'' / c = (ln c)'^2 + (ln c)'': the mean square of u under the weight cos^s.
    bracket = 1 + 0.5 * np.asarray(curvature, dtype=float) * (first**2 + second)
    # The expansion is of use where its correction is small; where it is not (a
    # slope below 1, which no rate here has at a field that ionises) it is held.
    return np.asarray(power_law_average(slope) * np.maximum(bracket, 0.5))


def _log_sum_exp(terms: Float[np.ndarray, " c *shape"]) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of the sum of the exponentials over the first axis.

    The arithmetic of ``scipy.special.logsumexp(terms, axis=0)`` (SciPy 1.18), step by
    step, so that the result is the same to the last bit, without its dispatch, which
    costs some 70 µs a call, more than the rates' arithmetic on a few fields: the
    largest term is taken out of the sum (its ties counted), and a result that is not
    finite is the sum taken directly.

    Parameters
    ----------
    terms : numpy.ndarray
        Logarithms, shape ``(c, ...)``.

    Returns
    -------
    numpy.ndarray
        The logarithm of the sum over the first axis, shape ``(...)``.
    """
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        peak = np.max(terms, axis=0)
        ties = terms == peak
        count = np.sum(ties, axis=0, dtype=terms.dtype)
        rest = np.sum(np.exp(np.where(ties, -np.inf, terms) - peak), axis=0)
        rest = np.where(rest == 0, rest, rest / count)
        out = np.log1p(rest) + np.log(count) + peak
        finite = np.isfinite(out)
        if not np.all(finite):
            out = np.where(finite, out, np.log(np.sum(np.exp(terms), axis=0)))
    return np.asarray(out)


@functools.cache
def _average_rule() -> tuple[np.ndarray, np.ndarray]:
    """Return the nodes and weights of the average over a cycle, on t in [0, 5].

    Returns
    -------
    tuple of numpy.ndarray
        :data:`_AVERAGE_NODES` Gauss-Legendre nodes and weights scaled to
        :data:`_AVERAGE_SPAN`, read-only.
    """
    nodes, weights = np.polynomial.legendre.leggauss(_AVERAGE_NODES)
    nodes = (nodes + 1) * (_AVERAGE_SPAN / 2)
    weights = weights * (_AVERAGE_SPAN / 2)
    nodes.flags.writeable = False
    weights.flags.writeable = False
    return nodes, weights


def _log_cycle_averages(
    log_terms: Callable[[Float[np.ndarray, " n"]], tuple[np.ndarray, ...]],
    field_au: Float[np.ndarray, " n"],
    *,
    centre: tuple[np.ndarray, ...] | None = None,
    log_term: Callable[[int, Float[np.ndarray, " m"]], np.ndarray] | None = None,
) -> tuple[Float[np.ndarray, " n"], ...]:
    r"""Return the logarithms of the averages over a cycle of rates, over their peaks.

    For each rate :math:`W` of ``log_terms``, :math:`\overline{W}/W(F) = (2/\pi)
    \int_0^{\pi/2} W(F\cos\varphi)/W(F)\, d\varphi`. With :math:`\cos\varphi =
    e^{-t^2}` the integrand,
    :math:`W(F e^{-t^2})/W(F) \cdot 2t\,e^{-t^2}/(1 - e^{-2t^2})^{1/2}`, is smooth at
    :math:`t = 0` and falls as :math:`e^{-(1+s)t^2}` with the rate's logarithmic slope
    :math:`s`. Gauss-Legendre quadrature over :math:`t \sqrt{1 + s} \in [0, 5]`, with
    the slope of each rate (or 0 where it falls), integrates it to 3e-5 or better up to
    three times the maximum field of the noble gases' models (TC-828). Unlike the
    expansion of :func:`cycle_average`, it holds where the rate bends sharply, near
    the barrier-suppression field and above, and it is free of the rounding of a
    second difference.

    The nodes of each rate are its own, so each rate is evaluated alone at them, by
    ``log_term`` where given: at 32 nodes for each of two rates that is most of the
    work, and the expensive rate need not be evaluated at the nodes of the cheap one.

    Parameters
    ----------
    log_terms : callable
        Fields (atomic units, shape ``(m,)``) to a tuple of logarithms of rates, each
        of shape ``(m,)``.
    field_au : numpy.ndarray
        Amplitudes, atomic units, shape ``(n,)``, > 0.
    centre : tuple of numpy.ndarray, optional
        ``log_terms(field_au)``, where the caller has it.
    log_term : callable, optional
        ``(k, fields)`` to ``log_terms(fields)[k]``, the same numbers, evaluating that
        rate alone.

    Returns
    -------
    tuple of numpy.ndarray
        The logarithm of each rate's average over its value at the amplitude, shape
        ``(n,)``.
    """
    if centre is None:
        centre = log_terms(field_au)
    evaluate = log_term or (lambda k, fields: log_terms(fields)[k])
    above = log_terms(field_au * (1 + _SLOPE_STEP))
    nodes, weights = _average_rule()
    averages = []
    for k, log_peak in enumerate(centre):
        slope = (above[k] - log_peak) / math.log1p(_SLOPE_STEP)
        scale = 1 / np.sqrt(1 + np.maximum(slope, 0.0))
        t = scale[:, None] * nodes[None, :]
        square = t**2
        # dphi/dt for cos(phi) = exp(-t^2), whose limit at t = 0 is sqrt(2).
        cosine = np.exp(-square)
        jacobian = 2 * t * cosine / np.sqrt(-np.expm1(-2 * square))
        values = evaluate(k, (field_au[:, None] * cosine).ravel())
        ratio = np.exp(values.reshape(t.shape) - log_peak[:, None])
        average = (2 / math.pi) * scale * np.sum(weights * jacobian * ratio, axis=1)
        averages.append(np.log(average))
    return tuple(averages)


def short_range_rate_au(  # noqa: PLR0913 - the arguments of the formula
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    coefficient_squared: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
    m: int = 0,
    *,
    channel_sum: ChannelSum = "integral",
) -> Float[np.ndarray, " *shape"]:
    r"""Return the cycle-averaged rate of one electron in a short-range potential.

    Perelomov et al. (1966), eq. (54):

    .. math::

       w_{lm}(F, \omega) = I_p\, |C_{\kappa l}|^2 \sqrt{\frac{6}{\pi}}\, f(l, m)
           \left(\frac{F \sqrt{1 + \gamma^2}}{2F_0}\right)^{|m| + 3/2}
           A_m(\omega, \gamma)\,
           \exp\left[-\frac{2F_0}{3F}\, g(\gamma)\right].

    Ilkov et al. (1992) print the prefactor of their eq. (8) as
    :math:`(3/2\pi)^{1/2}`, which is half of this and does not reduce to their own
    eq. (2) (TC-656); Larochelle et al.'s eq. (10) and Bergé et al.'s (B14) agree with
    Perelomov et al.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.
    coefficient_squared : float
        :math:`|C_{\kappa l}|^2`.
    l : int
        Orbital angular momentum of the electron.
    m : int, optional
        Its magnetic quantum number along the field.
    channel_sum : {"integral", "sum"}, optional
        How :func:`channel_factor` takes the sum over channels.

    Returns
    -------
    numpy.ndarray
        Rate, atomic units, with the shape of ``field_au``.

    Examples
    --------
    In the tunnelling limit it is the cycle average of the static rate without its
    Coulomb factor :math:`(2F_0/F)^{2n^*}`; for hydrogen at :math:`\gamma = 0.002`:

    >>> import numpy as np
    >>> field = np.array(0.05)
    >>> short = short_range_rate_au(field, 1e-4, 0.5, 4.0, 0)
    >>> static = tunnelling_rate_au(field, 0.5, 4.0, 0) / (2 / 0.05) ** 2
    >>> f"{float(short / (static * cycle_average_factor(field, 0.5))):.4f}"
    '1.0000'
    """
    return np.asarray(
        np.exp(
            _log_short_range_rate(
                np.asarray(field_au, dtype=float),
                omega_au,
                ionization_energy_au,
                coefficient_squared,
                l,
                m,
                channel_sum,
            )
        )
    )


def _log_short_range_rate(  # noqa: PLR0913, PLR0917 - the arguments of the formula
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    coefficient_squared: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
    m: int,
    channel_sum: str,
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of :func:`short_range_rate_au`.

    The rate spans hundreds of decades over the fields of a pulse, and its exponential
    underflows where the product with the prefactor is still representable, so the
    rates are assembled as logarithms and exponentiated once.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    coefficient_squared : float
        ``|C_kappa l|^2``.
    l : int
        Orbital angular momentum.
    m : int
        Magnetic quantum number.
    channel_sum : str
        ``integral`` or ``sum``.

    Returns
    -------
    numpy.ndarray
        The logarithm of the rate in atomic units.
    """
    return _log_short_range_strength(
        field_au,
        omega_au,
        ionization_energy_au,
        math.log(coefficient_squared * _angular_factor(l, m)),
        m,
        channel_sum,
    )


def _log_short_range_strength(
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    log_strength: float,
    m: int,
    channel_sum: str,
) -> Float[np.ndarray, " *shape"]:
    r"""Return the logarithm of the short-range rate with its prefactor given.

    :func:`short_range_rate_au` with :math:`|C_{\kappa l}|^2 f(l, m)` replaced by
    :math:`e^\text{log\_strength}`: the orbital enters the rate only through that
    product, and a molecule's orbital through its structure factor
    (:mod:`optmat.gas.molecular_ionisation_rates`).

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units, > 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    log_strength : float
        The logarithm of :math:`|C_{\kappa l}|^2 f(l, m)`.
    m : int
        Magnetic quantum number along the field.
    channel_sum : str
        ``integral`` or ``sum``.

    Returns
    -------
    numpy.ndarray
        The logarithm of the rate in atomic units.
    """
    order = abs(m)
    f0 = (2 * ionization_energy_au) ** 1.5
    gamma = keldysh_parameter(field_au, omega_au, ionization_energy_au)
    square = gamma**2
    if channel_sum == "integral" and not np.any(gamma < _SERIES_LIMIT):
        factor, exponent = _closed_forms(order, gamma, square)
    else:
        threshold = None
        if channel_sum == "sum":
            threshold = ionization_energy_au / omega_au * (1 + 0.5 / gamma**2)
        factor = channel_factor(m, gamma, threshold)
        exponent = ppt_exponent_function(gamma)
    return np.asarray(
        math.log(ionization_energy_au * math.sqrt(6 / math.pi))
        + log_strength
        + (order + 1.5) * (np.log(field_au / (2 * f0)) + 0.5 * np.log1p(square))
        + np.log(factor)
        - 2 * f0 / (3 * field_au) * exponent
    )


def _closed_forms(
    order: int,
    gamma: Float[np.ndarray, " *shape"],
    square: Float[np.ndarray, " *shape"],
) -> tuple[Float[np.ndarray, " *shape"], Float[np.ndarray, " *shape"]]:
    r"""Return :func:`channel_factor`'s integral and :func:`ppt_exponent_function`.

    Their closed forms, for Keldysh parameters all at or above
    :data:`_SERIES_LIMIT`, with :math:`\operatorname{arsinh}\gamma` and
    :math:`(1 + \gamma^2)^{1/2}` evaluated once for both: the same arithmetic,
    to the last bit, at a third less cost (the rates evaluate them some 34 times for
    each field, in their average over a cycle).

    Parameters
    ----------
    order : int
        :math:`|m|`.
    gamma : numpy.ndarray
        Keldysh parameter, >= :data:`_SERIES_LIMIT`.
    square : numpy.ndarray
        ``gamma**2``.

    Returns
    -------
    tuple of numpy.ndarray
        :math:`A_m(\omega, \gamma)` with the integral, and :math:`g(\gamma)`.
    """
    root = np.sqrt(1 + square)
    arsinh = np.arcsinh(gamma)
    alpha = 2 * (arsinh - gamma / root)
    beta = 2 * gamma / root
    total = (
        0.5
        * math.sqrt(math.pi)
        * math.factorial(order)
        * beta**order
        / (alpha + beta) ** (order + 1)
        * np.sqrt(beta / alpha)
    )
    factor = 4 / math.sqrt(3 * math.pi) / math.factorial(order) * square / (1 + square)
    exponent = 1.5 / gamma * ((1 + 0.5 / square) * arsinh - root / (2 * gamma))
    return np.asarray(factor * total), np.asarray(exponent)


def coulomb_factor(
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    form: CoulombFactor = "static",
    charge: float = 1.0,
) -> Float[np.ndarray, " *shape"]:
    r"""Return the factor by which the Coulomb field of the ion raises the rate.

    Three forms are in the literature (TC-653), and a fourth is optmat's fit:

    ``static``
        :math:`Q = (2F_0/F)^{2n^*}`: Perelomov and Popov (1967), eq. (42), for linear
        polarisation, of which they write that "the dependence ... on the frequency
        has dropped out"; Perelomov, Popov and Kuznetsov (1968), eqs. (1) and (11), that
        it "does not depend on :math:`\gamma` at all". It holds while the Coulomb field
        does not distort the path of the electron under the barrier,
        :math:`\gamma \ll \gamma_k = 2I_p/n^*\omega` (their eq. 20; 22 for argon and 43
        for helium at 800 nm), and more exactly while :math:`n^*\omega/2I_p` is below a
        critical value that falls with :math:`\gamma` (their Table I; TC-670). It was
        the default until TC-848. Against time-dependent Schrödinger calculations of one
        electron it is the closest of the three below a Keldysh parameter of about 1
        but too high above it, by 1.1 to 1.7 for argon at 800 nm, 1.4 to 3.9 at 400
        nm, 1.75 to 2.2 for krypton and xenon at 800 nm and, with Tong and Lin's static
        correction, 1.1 to 1.5 for hydrogen at 1200 to 600 nm (TC-812, TC-824, TC-825,
        TC-827);
        for the helium atom it agrees with Parker et al. (2007) at 780 nm.
    ``ilkov``
        :math:`Q\,(1 + \gamma^2)^{-n^*}`: the form of Ilkov et al. (1992), eq. (8), and
        Bergé et al. (2007), eq. (B14), in which the field of the Coulomb factor is
        multiplied by :math:`\sqrt{1 + \gamma^2}` together with that of the short-range
        rate. No paper derives it.
    ``popruzhenko``
        :math:`Q\,(1 + 2\gamma/e)^{-2n^*}`: Popruzhenko et al. (2008), eq. (5), derived
        for any :math:`\gamma`. Against the one-electron calculations it is too low
        above :math:`\gamma \approx 1`, by 2 to 4.
    ``frequency``
        :math:`Q\,(1 + 2\gamma/e)^{-2n^* c s}` with :math:`s = n^*\omega/2I_p`, the
        parameter of the static factor's range of validity (Perelomov, Popov and
        Kuznetsov 1968), and :math:`c = 5.25` fitted to one-electron calculations of
        argon at 400, 800 and 1030 nm: an empirical form, not a derived one. Every
        first-order Coulomb correction of the theory depends on :math:`\gamma` alone;
        the calculations' excess over the static factor also grows with the frequency
        at fixed :math:`\gamma`, as a term of order :math:`s` would. It tends to the
        static factor as :math:`\omega \to 0` (TC-831). It is the default (TC-848).

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units, >= 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.
    form : {"static", "ilkov", "popruzhenko", "frequency"}, optional
        The form.
    charge : float, optional
        Charge :math:`Z` of the ion left behind.

    Returns
    -------
    numpy.ndarray
        The factor, with the shape of ``field_au``.

    Raises
    ------
    ValueError
        If the form is not one of the four.

    Examples
    --------
    At :math:`\gamma = 1` the other two forms are below the static one:

    >>> import numpy as np
    >>> field, omega, energy = np.array(0.05), 0.05, 0.5
    >>> static = coulomb_factor(field, omega, energy)
    >>> f"{float(coulomb_factor(field, omega, energy, 'ilkov') / static):.4f}"
    '0.5000'
    >>> f"{float(coulomb_factor(field, omega, energy, 'popruzhenko') / static):.4f}"
    '0.3319'
    >>> f"{float(coulomb_factor(field, omega, energy, 'frequency') / static):.4f}"
    '0.7486'
    """
    return np.asarray(
        np.exp(
            _log_coulomb_factor(
                np.asarray(field_au, dtype=float),
                omega_au,
                ionization_energy_au,
                form,
                charge,
            )
        )
    )


def _log_coulomb_factor(
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    form: str,
    charge: float,
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of :func:`coulomb_factor`.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field amplitude, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units.
    ionization_energy_au : float
        Ionisation energy, hartree.
    form : str
        ``static``, ``ilkov``, ``popruzhenko`` or ``frequency``.
    charge : float
        Charge of the ion.

    Returns
    -------
    numpy.ndarray
        The logarithm of the factor.

    Raises
    ------
    ValueError
        If the form is not one of the four.
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    n_star = charge / kappa
    gamma = keldysh_parameter(field_au, omega_au, ionization_energy_au)
    static = 2 * n_star * np.log(2 * kappa**3 / field_au)
    static = np.asarray(static)
    if form == "static":
        return static
    if form == "ilkov":
        return np.asarray(static - n_star * np.log1p(gamma**2))
    if form == "popruzhenko":
        return np.asarray(static - 2 * n_star * np.log1p(2 * gamma / math.e))
    if form == "frequency":
        s = n_star * omega_au / (2 * ionization_energy_au)
        exponent = 2 * n_star * _FREQUENCY_CONSTANT * s
        return np.asarray(static - exponent * np.log1p(2 * gamma / math.e))
    msg = f"unknown Coulomb factor {form!r}; use one of {get_args(CoulombFactor)}"
    raise ValueError(msg)


def circular_saddle_point(gamma: ArrayLike) -> Float[np.ndarray, " *shape"]:
    r"""Return the saddle point :math:`\zeta_0(\gamma)` of a circularly polarised field.

    The root in :math:`(0, 1)` of Barth & Smirnova's (2013) eq. (94),

    .. math:: \operatorname{artanh} u = \frac{u}{1 - \zeta_0}, \qquad
       u = \sqrt{\frac{\zeta_0^2 + \gamma^2}{1 + \gamma^2}},

    found by bisection; :math:`\zeta_0 \approx \gamma^2/3` for :math:`\gamma \ll 1`
    (their Appendix F). It fixes the photon number of the most probable n-photon
    process, :math:`2 n_0/(1 + \zeta_0)`.

    Parameters
    ----------
    gamma : array_like
        Keldysh parameter, > 0.

    Returns
    -------
    numpy.ndarray
        :math:`\zeta_0`, with the shape of ``gamma``.

    Examples
    --------
    >>> f"{float(circular_saddle_point(0.05)) / (0.05**2 / 3):.4f}"
    '0.9984'
    """
    g = np.asarray(gamma, dtype=float)
    low = np.zeros_like(g)
    high = np.ones_like(g)
    for _ in range(80):
        middle = 0.5 * (low + high)
        u = np.sqrt((middle**2 + g**2) / (1 + g**2))
        # artanh u = ln(1 + u) - ln(1 - u^2)/2, with 1 - u^2 = (1 - zeta^2)/(1 + g^2):
        # finite where u rounds to 1 (g of 1e7 in a pulse's far wings).
        artanh = np.log1p(u) - 0.5 * (np.log1p(-(middle**2)) - np.log1p(g**2))
        positive = artanh > u / (1 - middle)
        low = np.where(positive, middle, low)
        high = np.where(positive, high, middle)
    return np.asarray(0.5 * (low + high))


def _log_circular_terms(
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
) -> Float[np.ndarray, " *shape"]:
    r"""Return the logarithm of the shell's rate in a circular field, short range.

    Barth & Smirnova's (2013) simple formulas, eqs. (97) and (100) to (105), without
    the factor :math:`|C_{\kappa l}|^2 I_p F/2F_0` common to every orbital: for an s
    orbital :math:`h^s e^{-(2F_0/3F) g}`, for a p shell the sum over its three
    orbitals :math:`(h^{p_+} + h^{p_-} + h^{p_0}) e^{-(2F_0/3F) g}`, the magnetic
    quantum numbers along the direction of propagation.

    Parameters
    ----------
    field_au : numpy.ndarray
        Magnitude of the turning field, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units, >= 0.
    ionization_energy_au : float
        Ionisation energy, hartree.
    l : int
        0 or 1.

    Returns
    -------
    numpy.ndarray
        The logarithm.

    Raises
    ------
    ValueError
        If ``l`` is neither 0 nor 1.
    """
    if l not in {0, 1}:
        msg = (
            f"the circular-polarisation rate is given for s and p orbitals, not l = {l}"
        )
        raise ValueError(msg)
    field = np.asarray(field_au, dtype=float)
    ratio = (2 * ionization_energy_au) ** 1.5 / field  # F_0 / F
    gamma = np.asarray(keldysh_parameter(field, omega_au, ionization_energy_au))
    # q = zeta_0 / gamma^2, finite at gamma = 0: below 1e-2 from Appendix F's series
    # (bisection cannot resolve zeta_0 ~ gamma^2/3 there in double precision).
    small = gamma < 1e-2
    safe = np.where(small, 1.0, gamma)
    zeta = np.where(small, gamma**2 / 3, circular_saddle_point(safe))
    q = np.where(small, (1 - 28 * gamma**2 / 45) / 3, zeta / safe**2)
    over = q * gamma  # zeta_0 / gamma
    g = 3 * q / (1 - zeta**2) * np.sqrt((1 + gamma**2) * (over**2 + 1))
    h_s = (1 - zeta) * np.sqrt(
        (1 + gamma**2) * (1 - zeta**2) / ((over**2 + 1) * (2 * over**2 + zeta**2 + 1))
    )
    if l == 0:
        prefactor = h_s
    else:
        root = np.sqrt((over**2 + 1) / (1 + gamma**2))
        h_pm = (
            h_s
            * 3
            * (1 + gamma**2)
            / (2 * (1 - zeta**2))
            * ((root - over) ** 2 + (root + over) ** 2)
        )
        h_0 = h_s * 1.5 / ratio * (1 - zeta) * np.sqrt((1 + gamma**2) / (over**2 + 1))
        prefactor = h_pm + h_0
    return np.asarray(np.log(prefactor) - 2 * ratio / 3 * g)


def circular_enhancement(
    field_au: ArrayLike,
    omega_au: float,
    ionization_energy_au: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
) -> Float[np.ndarray, " *shape"]:
    r"""Return the rate of a shell in a circular field over its adiabatic limit.

    The short-range rate of Barth & Smirnova (2013, eqs. 100 to 105), summed over the
    shell's orbitals, at the frequency over the same at zero frequency and the same
    magnitude of the field: the factor by which the field's turning raises the rate
    above the static rate at its magnitude. It tends to 1 for :math:`\gamma \to 0`
    (TC-820).

    Parameters
    ----------
    field_au : array_like
        Magnitude of the turning field, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units.
    ionization_energy_au : float
        Ionisation energy, hartree.
    l : int
        0 (s) or 1 (p).

    Returns
    -------
    numpy.ndarray
        The factor, with the shape of ``field_au``.

    Examples
    --------
    Argon at 800 nm and 2e14 W/cm^2 (a turning field of 0.0534 au):

    >>> f"{float(circular_enhancement(0.0534, 0.05695, 0.5792, 1)):.3f}"
    '3.709'
    """
    return np.asarray(
        np.exp(
            _log_circular_enhancement(
                np.asarray(field_au, dtype=float), omega_au, ionization_energy_au, l
            )
        )
    )


def _log_circular_enhancement(
    field_au: Float[np.ndarray, " *shape"],
    omega_au: float,
    ionization_energy_au: float,
    l: int,  # noqa: E741 - l is the orbital's symbol
) -> Float[np.ndarray, " *shape"]:
    """Return the logarithm of :func:`circular_enhancement`.

    It stays finite where the factor itself overflows (weak fields, large Keldysh
    parameters), so that it can be added to the logarithm of a small static rate.

    Parameters
    ----------
    field_au : numpy.ndarray
        Magnitude of the turning field, atomic units, > 0.
    omega_au : float
        Angular frequency, atomic units.
    ionization_energy_au : float
        Ionisation energy, hartree.
    l : int
        0 (s) or 1 (p).

    Returns
    -------
    numpy.ndarray
        The logarithm of the factor.
    """
    return np.asarray(
        _log_circular_terms(field_au, omega_au, ionization_energy_au, l)
        - _log_circular_terms(field_au, 0.0, ionization_energy_au, l)
    )


def tong_lin_factor(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    alpha: float,
    charge: float = 1.0,
) -> Float[np.ndarray, " *shape"]:
    r"""Return Tong and Lin's correction of the tunnelling rate at finite field.

    Tong and Lin (2005), eq. (2):

    .. math::

       S(F) = \exp\left(-\alpha\, \frac{Z^2}{I_p}\, \frac{F}{\kappa^3}\right),

    with :math:`\alpha` fitted to exact static rates of one electron in a model
    potential (their Table 2: 7 for He, 9 for Ne and Ar), relative to the tunnelling
    rate with Tong et al.'s (2002) coefficient :math:`C_l` and :math:`m = 0`. With
    :math:`\alpha = 9` it reproduces Zhang et al.'s (2014) exact rates of argon in Tong
    and Lin's potential within 7 % up to 1.7 times the barrier-suppression field
    (TC-655). The name "barrier
    suppression" misleads: the factor is 0.5 for argon at 0.64 of that field, and is
    needed at every field.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field, atomic units, >= 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.
    alpha : float
        Tong and Lin's constant.
    charge : float, optional
        Charge :math:`Z` of the ion left behind.

    Returns
    -------
    numpy.ndarray
        The factor, in (0, 1].

    Examples
    --------
    Argon at 1e14 W/cm\ :sup:`2`:

    >>> import numpy as np
    >>> f"{float(tong_lin_factor(np.array(0.05338), 0.57916, 9.0)):.3f}"
    '0.514'
    """
    kappa = math.sqrt(2 * ionization_energy_au)
    return np.asarray(
        np.exp(
            -alpha
            * charge**2
            / ionization_energy_au
            * np.asarray(field_au, dtype=float)
            / kappa**3
        )
    )


def zhang_factor(
    field_au: Float[np.ndarray, " *shape"],
    ionization_energy_au: float,
    coefficients: tuple[float, float, float],
    charge: float = 1.0,
) -> Float[np.ndarray, " *shape"]:
    r"""Return Zhang, Lan and Lu's correction of the tunnelling rate.

    Zhang et al. (2014), eq. (8), with the sign their own Table IX requires (the
    equation is printed with a minus sign before the polynomial; TC-509):

    .. math::

       S(F) = \exp\left(a_1 x^2 + a_2 x + a_3\right), \qquad x = F / F_b,
       \quad F_b = \frac{I_p^2}{4 Z},

    fitted for :math:`0.5 \le x \le 4.5`. The polynomial is positive at small
    :math:`x` for every gas but neon, where a correction of a rate that is already
    exact must be 1; the factor is held at 1 there, which leaves a corner in the rate.

    Parameters
    ----------
    field_au : numpy.ndarray
        Field, atomic units, >= 0.
    ionization_energy_au : float
        Ionisation energy, hartree, > 0.
    coefficients : tuple of float
        :math:`a_1, a_2, a_3` (their Table XI).
    charge : float, optional
        Charge :math:`Z` of the ion left behind.

    Returns
    -------
    numpy.ndarray
        The factor, in (0, 1].

    Examples
    --------
    Argon at its barrier-suppression field:

    >>> import numpy as np
    >>> barrier = 0.57916**2 / 4
    >>> coefficients = (0.16178, -1.50441, 0.32127)
    >>> f"{float(zhang_factor(np.array(barrier), 0.57916, coefficients)):.3f}"
    '0.360'
    """
    a1, a2, a3 = coefficients
    x = np.asarray(field_au, dtype=float) * 4 * charge / ionization_energy_au**2
    return np.asarray(np.exp(np.minimum(a1 * x**2 + a2 * x + a3, 0.0)))


@dataclass(frozen=True)
class Channel:
    """One level of the ion that the ionisation can leave, with its weight.

    Attributes
    ----------
    label : str
        The level, as the source names it.
    ionization_energy : float
        Energy of the level above the atom's ground state, J.
    weight : float
        Share of the ionising orbital that leaves the ion in this level; the weights
        of a model's channels sum to 1.
    source : str
        The curated datasets the energy comes from.
    """

    label: str
    ionization_energy: float
    weight: float
    source: str

    @property
    def ionization_energy_au(self) -> float:
        """The channel's ionisation energy in hartree.

        Returns
        -------
        float
            ``ionization_energy / E_h``.
        """
        return self.ionization_energy / AU_ENERGY


@dataclass(frozen=True)
class IonisationRateModel(Identified):
    r"""The ionisation rate of a gas at one wavelength, in one named model.

    Built by :func:`ionisation_rate_model` from the curated parameters, or by
    :meth:`from_parameters` from a caller's; the attributes say what the model is made
    of. :meth:`~optmat.identity.Identified.with_parameters` rebuilds it with any of
    :attr:`INPUTS` changed (OEP 1, section 6.3). Every field is an amplitude or an
    instantaneous value of a **linearly polarised** field.

    Attributes
    ----------
    gas : Gas or str
        The gas, or a custom species' name.
    wavelength : float
        Vacuum wavelength of the field, m; ``math.inf`` for a static field.
    model : str
        The model's name in ``ionisation_rates.toml``, or ``"user"``.
    channels : tuple of Channel
        The levels of the ion and their weights.
    orbital_l : int
        Orbital angular momentum of the ionising electron; :math:`m = 0`.
    structure : str
        ``tong_2002``: the coefficient :math:`C_l` of Tong et al. (2002), Table II.
        ``hartree``: Ammosov et al.'s (1986) eq. (19) with :math:`l^* = n^* - 1`.
    structure_coefficient : float or None
        :math:`C_l` in Tong et al.'s normalisation; None for ``hartree``.
    electrons : float
        :math:`N_\text{eff}`, the rate of the atom over that of one electron.
    coulomb : str
        The Coulomb factor (:func:`coulomb_factor`).
    static_correction : str
        ``tong-lin``, ``zhang`` or ``none``.
    static_constants : tuple of float
        Tong and Lin's :math:`\alpha`, or Zhang et al.'s :math:`a_1, a_2, a_3`.
    channel_sum : str
        ``integral`` (smooth) or ``sum`` (:func:`channel_factor`).
    averaging : str
        ``consistent`` or ``tunnelling``: how the cycle-averaged and the
        instantaneous rate are related (:data:`Averaging`).
    uncertainty_factor : float
        The factor within which the rate is expected to be right between
        :math:`10^{-6}` of saturation and saturation, at the wavelengths of
        :attr:`wavelength_range`: a judgement from the comparisons of the manual page,
        not a statistical interval.
    sources : dict of str to str
        The curated dataset or decision behind each constant.
    reasons : dict of str to str
        Why each choice was made (from ``ionisation_rates.toml``).

    Examples
    --------
    >>> model = ionisation_rate_model("Ar", 800e-9)
    >>> f"{float(model.cycle_averaged_rate(2.744e10)):.3e}"  # 1e14 W/cm^2
    '6.112e+11'
    >>> f"{float(model.instantaneous_rate(2.744e10)):.3e}"
    '2.576e+12'
    """

    gas: Gas | str
    wavelength: float
    model: str
    channels: tuple[Channel, ...]
    orbital_l: int
    structure: str
    structure_coefficient: float | None
    electrons: float
    coulomb: str
    static_correction: str
    static_constants: tuple[float, ...]
    channel_sum: str
    averaging: str
    uncertainty_factor: float
    sources: Mapping[str, str] = field(default_factory=dict)
    reasons: Mapping[str, str] = field(default_factory=dict)

    IDENTITY_KIND: ClassVar[str] = "ionisation_rate"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.ionisation_rates",)
    INPUTS: ClassVar[tuple[str, ...]] = (
        "averaging",
        "channel_sum",
        "channels",
        "coulomb",
        "electrons",
        "ionization_energy",
        "orbital_l",
        "static_constants",
        "static_correction",
        "structure",
        "structure_coefficient",
        "uncertainty_factor",
    )

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and source records.

        Returns
        -------
        tuple
            ``("rate/<model>/<gas>@<wavelength>", gas, sources)``, the wavelength
            in nm or ``static``.
        """
        at = (
            "static"
            if math.isinf(self.wavelength)
            else f"{self.wavelength * 1e9:.6g}nm"
        )
        name = f"rate/{self.model}/{self.gas}@{at}"
        return (
            name,
            str(self.gas),
            source_ids(*self.sources.values(), *(c.source for c in self.channels)),
        )

    def __post_init__(self) -> None:
        """Check the settings.

        Raises
        ------
        ValueError
            If a setting is not one the model knows; the wavelength, the number of
            electrons, a channel's energy or a Tong-Lin constant is not positive; the
            orbital's l is negative; or the weights do not sum to 1.
        """
        super().__post_init__()
        for name, value, allowed in (
            ("coulomb", self.coulomb, get_args(CoulombFactor)),
            ("static_correction", self.static_correction, get_args(StaticCorrection)),
            ("channel_sum", self.channel_sum, get_args(ChannelSum)),
            ("averaging", self.averaging, get_args(Averaging)),
            ("structure", self.structure, get_args(Structure)),
        ):
            if value not in allowed:
                msg = f"{name} = {value!r} is not one of {allowed}"
                raise ValueError(msg)
        if not self.wavelength > 0:
            msg = f"the wavelength must be positive, not {self.wavelength!r}"
            raise ValueError(msg)
        if not self.electrons > 0:
            msg = f"electrons must be positive, not {self.electrons!r}"
            raise ValueError(msg)
        if self.orbital_l < 0:
            msg = f"orbital_l must not be negative, not {self.orbital_l}"
            raise ValueError(msg)
        if not self.channels or not all(c.ionization_energy > 0 for c in self.channels):
            msg = "every channel needs a positive ionisation energy"
            raise ValueError(msg)
        if not math.isclose(sum(c.weight for c in self.channels), 1.0, abs_tol=1e-12):
            msg = "the weights of the channels must sum to 1"
            raise ValueError(msg)
        expected = {"tong-lin": 1, "zhang": 3, "none": 0}[self.static_correction]
        if len(self.static_constants) != expected:
            msg = f"{self.static_correction} takes {expected} constants"
            raise ValueError(msg)
        if self.static_correction == "tong-lin" and not self.static_constants[0] > 0:
            msg = (
                f"Tong and Lin's alpha must be positive, not {self.static_constants[0]}"
            )
            raise ValueError(msg)
        if (self.structure == "tong_2002") != (self.structure_coefficient is not None):
            msg = "structure_coefficient goes with structure = 'tong_2002' only"
            raise ValueError(msg)
        if self.wavelength < _MINIMUM_WAVELENGTH:
            warnings.warn(
                f"{self.wavelength * 1e9:.0f} nm is below the shortest wavelength, "
                f"{_MINIMUM_WAVELENGTH * 1e9:.0f} nm, at which the ionisation rates "
                "were compared with references",
                ModelRangeWarning,
                stacklevel=3,
            )

    @property
    def omega_au(self) -> float:
        """The angular frequency of the field in atomic units.

        Returns
        -------
        float
            ``2 pi c / wavelength`` over the atomic unit of angular frequency; 0 for a
            static field.
        """
        return 2 * math.pi * C / self.wavelength / AU_ANGULAR_FREQUENCY

    @property
    def wavelength_range(self) -> tuple[float, float]:
        """The wavelengths over which the model was compared with references.

        Returns
        -------
        tuple of (float, float)
            The shortest and the longest wavelength, m (the longest is a static
            field).
        """
        return (_MINIMUM_WAVELENGTH, math.inf)

    @property
    def barrier_suppression_field(self) -> float:
        r"""The field at which the barrier falls to the lowest channel's level.

        :math:`F_b = I_p^2 / 4Z` in atomic units (Tong and Lin 2005, eq. 1, as
        :math:`\kappa^4/16Z`).

        Returns
        -------
        float
            The field, V/m.
        """
        energy = min(c.ionization_energy_au for c in self.channels)
        return energy**2 / 4 * AU_FIELD

    @property
    def maximum_field(self) -> float:
        """The largest field at which the static correction was fitted.

        Returns
        -------
        float
            V/m: twice the barrier-suppression field for Tong and Lin's factor, 4.5
            times for Zhang et al.'s, and the barrier-suppression field itself
            without a correction.
        """
        factor = {
            "tong-lin": _TONG_LIN_MAXIMUM,
            "zhang": _ZHANG_RANGE[1],
            "none": 1.0,
        }[self.static_correction]
        return factor * self.barrier_suppression_field

    @property
    def citation(self) -> str:
        """The references of the model as it is set.

        Returns
        -------
        str
            One line naming the papers of each factor.
        """
        coulomb = {
            "static": "Perelomov & Popov, Sov. Phys. JETP 25, 336 (1967), eq. (42), "
            "as in Larochelle et al., J. Phys. B 31, 1215 (1998), eq. (10)",
            "ilkov": "Ilkov et al., J. Phys. B 25, 4005 (1992), eq. (8); Bergé et "
            "al., Rep. Prog. Phys. 70, 1633 (2007), eq. (B14)",
            "frequency": "Perelomov & Popov, Sov. Phys. JETP 25, 336 (1967), eq. (42), "
            "as in Larochelle et al., J. Phys. B 31, 1215 (1998), eq. (10), times "
            "(1 + 2 gamma/e)^(-2 n* c s), s = n* omega/2 I_p, c = 5.25, optmat's fit "
            "to one-electron calculations (TC-831)",
            "popruzhenko": "Popruzhenko et al., Phys. Rev. Lett. 101, 193003 (2008), "
            "eq. (5)",
        }[self.coulomb]
        static = {
            "tong-lin": "; static correction of Tong & Lin, J. Phys. B 38, 2593 (2005)",
            "zhang": "; static correction of Zhang et al., Phys. Rev. A 90, 043410 "
            "(2014)",
            "none": "",
        }[self.static_correction]
        structure = {
            "tong_2002": "C_l of Tong et al., Phys. Rev. A 66, 033402 (2002)",
            "hartree": "C of Ammosov et al., Sov. Phys. JETP 64, 1191 (1986), eq. (19)",
        }[self.structure]
        given = f"{self.sources['citation']}; " if "citation" in self.sources else ""
        return (
            f"{given}Perelomov, Popov & Terent'ev, Sov. Phys. JETP 23, 924 (1966), eq. "
            f"(54); Coulomb factor of {coulomb}{static}; {structure}; "
            f"{self.electrons:g} effective electrons (optmat, model {self.model!r})"
        )

    def keldysh(self, field_strength: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the Keldysh parameter of the lowest channel.

        Parameters
        ----------
        field_strength : array_like
            Field amplitude, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            :math:`\gamma`, with the shape of ``field_strength``; ``inf`` at zero
            field.

        Examples
        --------
        >>> model = ionisation_rate_model("Ar", 800e-9)
        >>> f"{float(model.keldysh(2.744e10)):.3f}"
        '1.149'
        """
        amplitude = np.asarray(
            np.abs(np.asarray(field_strength, dtype=float)) / AU_FIELD
        )
        energy = min(c.ionization_energy_au for c in self.channels)
        with np.errstate(divide="ignore"):
            return keldysh_parameter(np.asarray(amplitude), self.omega_au, energy)

    def _coefficient_squared(self, energy_au: float) -> float:
        r"""Return :math:`|C_{\kappa l}|^2` of a channel.

        Tong et al.'s :math:`C_l` multiplies :math:`r^{n^*-1} e^{-\kappa r}` and
        Perelomov et al.'s :math:`C_{\kappa l}` multiplies :math:`\kappa^{3/2}
        (\kappa r)^{n^*-1} e^{-\kappa r}`, so :math:`|C_{\kappa l}|^2 = C_l^2 /
        \kappa^{2n^* + 1}`. The same :math:`C_l` is used for both levels of the ion.

        Parameters
        ----------
        energy_au : float
            The channel's ionisation energy, hartree.

        Returns
        -------
        float
            The coefficient squared.
        """
        kappa = math.sqrt(2 * energy_au)
        n_star = 1 / kappa
        if self.structure_coefficient is None:
            return hartree_coefficient_squared(n_star, n_star - 1)
        return self.structure_coefficient**2 / kappa ** (2 * n_star + 1)

    def _log_static_correction(
        self, field_au: Float[np.ndarray, " n"], energy_au: float
    ) -> Float[np.ndarray, " n"]:
        """Return the logarithm of the static correction of a channel.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``.
        energy_au : float
            The channel's ionisation energy, hartree.

        Returns
        -------
        numpy.ndarray
            The logarithm of the factor, <= 0.
        """
        if self.static_correction == "tong-lin":
            (alpha,) = self.static_constants
            kappa = math.sqrt(2 * energy_au)
            return np.asarray(-alpha / energy_au * field_au / kappa**3)
        if self.static_correction == "zhang":
            a1, a2, a3 = self.static_constants
            return np.log(zhang_factor(field_au, energy_au, (a1, a2, a3)))
        return np.zeros_like(field_au)

    def _log_terms(
        self, field_au: Float[np.ndarray, " n"]
    ) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
        """Return the logarithms of the channel sums, static and with the frequency.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        tuple of numpy.ndarray
            The logarithms of the sums over the channels of :meth:`_log_static_terms`
            and :meth:`_log_enhanced_terms`, shape ``(n,)``.
        """
        return self._log_term(0, field_au), self._log_term(1, field_au)

    def _log_term(
        self, which: int, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " n"]:
        """Return one of the two logarithms of :meth:`_log_terms`, evaluating it alone.

        Parameters
        ----------
        which : int
            0, the static sum, or 1, the sum with the frequency.
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        numpy.ndarray
            The logarithm, shape ``(n,)``.
        """
        if which == 0:
            return _log_sum_exp(self._log_static_terms(field_au))
        return _log_sum_exp(self._log_enhanced_terms(field_au))

    def _log_static_term(
        self, channel: Channel, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " n"]:
        """Return the logarithm of one channel's static term.

        Parameters
        ----------
        channel : Channel
            The channel.
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        numpy.ndarray
            The logarithm, shape ``(n,)``.
        """
        energy = channel.ionization_energy_au
        return (
            math.log(channel.weight)
            + _log_tunnelling_rate(
                field_au,
                energy,
                self._coefficient_squared(energy),
                self.orbital_l,
                0,
                1.0,
            )
            + self._log_static_correction(field_au, energy)
        )

    def _log_static_terms(
        self, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " c n"]:
        """Return each channel's static term, as a logarithm.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        numpy.ndarray
            Shape ``(channels, n)``.
        """
        return np.stack([self._log_static_term(c, field_au) for c in self.channels])

    def _log_enhanced_terms(
        self, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " c n"]:
        r"""Return each channel's term with the frequency, as a logarithm.

        The static rate of the atom's one electron is :math:`\sum_c p_c\,
        w_{\text{stat},c} S_c`. The second is :math:`\sum_c p_c\, w_{\text{stat},c}
        S_c E_c` with the enhancement of the rate by the field's frequency,

        .. math::

           E(F, \omega) = \frac{Q\, w^\text{sr}(F, \omega)}
               {w_\text{stat}(F)\, (3F/\pi F_0)^{1/2}},

        the cycle-averaged rate of the theory over its own tunnelling limit
        (Perelomov et al. 1966, eq. 7), which tends to 1 for :math:`\gamma \to 0`.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        numpy.ndarray
            Shape ``(channels, n)``; the static terms in a static field.
        """
        omega = self.omega_au
        if omega == 0:
            return self._log_static_terms(field_au)
        enhanced = []
        for channel in self.channels:
            energy = channel.ionization_energy_au
            squared = self._coefficient_squared(energy)
            log_cycle_averaged = _log_short_range_rate(
                field_au, omega, energy, squared, self.orbital_l, 0, self.channel_sum
            ) + _log_coulomb_factor(field_au, omega, energy, self.coulomb, 1.0)
            enhanced.append(
                math.log(channel.weight)
                + log_cycle_averaged
                - np.log(cycle_average_factor(field_au, energy))
                + self._log_static_correction(field_au, energy)
            )
        return np.stack(enhanced)

    def _log_channel_rates_au(
        self, field_au: Float[np.ndarray, " n"]
    ) -> Float[np.ndarray, " c n"]:
        """Return the logarithm of each channel's instantaneous rate, atomic units.

        Each channel's term times the one factor the model applies to their sum
        (N_eff and the averaging), so that the channel rates sum to the rate.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.

        Returns
        -------
        numpy.ndarray
            Shape ``(channels, n)``.
        """
        enhanced = self._log_enhanced_terms(field_au)
        total = self._log_rate_au(field_au, cycle_averaged=False)
        return enhanced + (total - _log_sum_exp(enhanced))[None, :]

    def _log_rate_au(
        self, field_au: Float[np.ndarray, " n"], *, cycle_averaged: bool
    ) -> Float[np.ndarray, " n"]:
        r"""Return the logarithm of the rate, atomic units, at positive fields.

        With :math:`B = \sum_c p_c\, w_{\text{stat},c} S_c E_c` of
        :meth:`_log_terms`, its static limit :math:`B_0` (:math:`E = 1`) and their
        logarithmic slopes :math:`s = d\ln B/d\ln F` and :math:`s_0`:

        ``averaging = "consistent"``
            the cycle-averaged rate is :math:`N_\text{eff}\, B\, a_0` and the
            instantaneous rate :math:`N_\text{eff}\, B\, a_0/a`, with :math:`a`
            and :math:`a_0` the averages over a cycle of :math:`B` and :math:`B_0`
            over their values at the amplitude, by quadrature (TC-828);
        ``averaging = "tunnelling"``
            the cycle-averaged rate is :math:`N_\text{eff}\, B\, (3F/\pi
            F_0)^{1/2}`, the theory's own, and the instantaneous rate
            :math:`N_\text{eff}\, B`.

        Parameters
        ----------
        field_au : numpy.ndarray
            Field, atomic units, shape ``(n,)``, > 0.
        cycle_averaged : bool
            The rate against the amplitude, or the instantaneous rate.

        Returns
        -------
        numpy.ndarray
            The logarithm of the rate, shape ``(n,)``.
        """
        centre = self._log_terms(field_au)
        log_static, log_enhanced = centre
        log_rate = log_enhanced + math.log(self.electrons)
        if self.averaging == "tunnelling":
            if cycle_averaged:
                energy = min(c.ionization_energy_au for c in self.channels)
                log_rate = log_rate + np.log(cycle_average_factor(field_au, energy))
            return np.asarray(log_rate)
        if not cycle_averaged and self.omega_au == 0:
            # A static field: the rate itself.
            return np.asarray(log_static + math.log(self.electrons))
        log_static_average, log_average = _log_cycle_averages(
            self._log_terms, field_au, centre=centre, log_term=self._log_term
        )
        log_rate = log_rate + log_static_average
        if not cycle_averaged:
            log_rate = log_rate - log_average
        return np.asarray(log_rate)

    def _held(self, strength: np.ndarray, *, stacklevel: int) -> np.ndarray:
        """Return the fields at which to evaluate, warning of any above the range.

        Above :attr:`maximum_field` a statically corrected rate is held at its value
        there: the correction is extrapolated beyond its fit, and the rate would turn
        over (TC-816, TC-828). A model without a correction keeps its formula, which
        rises with the field.

        Parameters
        ----------
        strength : numpy.ndarray
            Magnitudes of the field, V/m.
        stacklevel : int
            The warning's stack level, at the caller of the public method.

        Returns
        -------
        numpy.ndarray
            The fields, those above the range replaced by its end where held.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.
        """
        if not np.any(strength > self.maximum_field):
            return strength
        held = self.static_correction != "none"
        warnings.warn(
            f"fields above {self.maximum_field:.3g} V/m are beyond the range over "
            f"which the static correction ({self.static_correction}) of "
            f"{self.gas} was fitted"
            + ("; the rate is held at its value there" if held else ""),
            ModelRangeWarning,
            stacklevel=stacklevel,
        )
        return (
            np.asarray(np.minimum(strength, self.maximum_field)) if held else strength
        )

    def _rate(
        self, field_strength: ArrayLike, *, cycle_averaged: bool
    ) -> Float[np.ndarray, " *shape"]:
        """Return a rate in SI, zero at zero field, held above the fitted range.

        Parameters
        ----------
        field_strength : array_like
            Field, V/m.
        cycle_averaged : bool
            The rate against the amplitude, or the instantaneous rate.

        Returns
        -------
        numpy.ndarray
            Rate per atom, 1/s.
        """
        strength = np.asarray(np.abs(np.asarray(field_strength, dtype=float)))
        strength = self._held(strength, stacklevel=4)
        flat = strength.ravel() / AU_FIELD
        rate = np.zeros_like(flat)
        positive = flat > 0
        if np.any(positive):
            with np.errstate(under="ignore"):
                rate[positive] = np.exp(
                    self._log_rate_au(flat[positive], cycle_averaged=cycle_averaged)
                )
        return (rate * AU_ANGULAR_FREQUENCY).reshape(strength.shape)

    def cycle_averaged_rate(self, amplitude: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the rate averaged over a cycle of a field of constant amplitude.

        Parameters
        ----------
        amplitude : array_like
            Amplitude of the field, V/m, :math:`(2I/\varepsilon_0 c)^{1/2}` for an
            intensity :math:`I`; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rate per atom, 1/s, with the shape of ``amplitude``; 0 at zero field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        Helium at 800 nm and 1e15 W/cm\ :sup:`2`:

        >>> model = ionisation_rate_model("He", 800e-9)
        >>> f"{float(model.cycle_averaged_rate(8.68e10)):.3e}"
        '1.568e+13'
        """
        return self._rate(amplitude, cycle_averaged=True)

    def instantaneous_rate(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the rate at the instantaneous value of a carrier-resolved field.

        The cycle-averaged rate at an amplitude equal to the instantaneous field,
        divided by the average over a cycle that turned the rate into it: with the
        default ``averaging = "consistent"``, the model's own average
        (the average over a cycle, by quadrature, of its rate as a function of the
        field; TC-828); with ``averaging = "tunnelling"``, :func:`cycle_average_factor`
        (TC-733). The Keldysh parameter is that of the instantaneous field.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rate per atom, 1/s, with the shape of ``field_strength``; 0 at zero field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        >>> model = ionisation_rate_model("He", 800e-9)
        >>> f"{float(model.instantaneous_rate(8.68e10)):.3e}"
        '5.968e+13'
        """
        return self._rate(field_strength, cycle_averaged=False)

    def _positive_fields(
        self, field_strength: ArrayLike
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the fields' magnitudes, their positive ones in atomic units, a mask.

        Parameters
        ----------
        field_strength : array_like
            Field, V/m; its sign is ignored.

        Returns
        -------
        tuple of numpy.ndarray
            ``|F|`` (V/m), the positive fields flattened (atomic units), and the mask
            of the positive ones in the flattened array.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.
        """
        strength = np.asarray(np.abs(np.asarray(field_strength, dtype=float)))
        strength = self._held(strength, stacklevel=5)
        flat = strength.ravel() / AU_FIELD
        positive = flat > 0
        return strength, flat[positive], positive

    def circular_rate(self, field_strength: ArrayLike) -> Float[np.ndarray, " *shape"]:
        r"""Return the rate in a circularly polarised field of a given magnitude.

        The field of circular polarisation turns at constant magnitude :math:`F =
        (I/\varepsilon_0 c)^{1/2}`, so the rate is the same at every instant. Each
        channel's static rate (its tunnelling rate with its coefficient, the static
        correction and the effective number of electrons, as for linear
        polarisation) is multiplied by the factor by which the field's turning
        raises the shell's short-range rate, :func:`circular_enhancement` (Barth &
        Smirnova 2013), and by the model's Coulomb factor over its static limit. In
        a static field it is :meth:`instantaneous_rate` at the magnitude (TC-820).

        Parameters
        ----------
        field_strength : array_like
            Magnitude of the turning field, V/m.

        Returns
        -------
        numpy.ndarray
            Rate per atom, 1/s, with the shape of ``field_strength``; 0 at zero
            field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        Argon at 800 nm, 2e14 W/cm^2 in circular polarisation, over the static rate
        at the same magnitude of the field:

        >>> model = ionisation_rate_model("Ar", 800e-9)
        >>> static = ionisation_rate_model("Ar", float("inf"))
        >>> field = (2e18 / (8.8541878128e-12 * 299792458.0)) ** 0.5
        >>> ratio = model.circular_rate(field) / static.instantaneous_rate(field)
        >>> f"{float(ratio):.2f}"
        '2.85'
        """
        strength = np.asarray(np.abs(np.asarray(field_strength, dtype=float)))
        strength = self._held(strength, stacklevel=3)
        flat = strength.ravel() / AU_FIELD
        rate = np.zeros_like(flat)
        positive = flat > 0
        if np.any(positive):
            field = flat[positive]
            terms = []
            for channel in self.channels:
                energy = channel.ionization_energy_au
                log_term = (
                    math.log(channel.weight)
                    + _log_tunnelling_rate(
                        field,
                        energy,
                        self._coefficient_squared(energy),
                        self.orbital_l,
                        0,
                        1.0,
                    )
                    + self._log_static_correction(field, energy)
                )
                if self.omega_au > 0:
                    log_term = (
                        log_term
                        + _log_circular_enhancement(
                            field, self.omega_au, energy, self.orbital_l
                        )
                        + _log_coulomb_factor(
                            field, self.omega_au, energy, self.coulomb, 1.0
                        )
                        - _log_coulomb_factor(field, 0.0, energy, "static", 1.0)
                    )
                terms.append(log_term)
            with np.errstate(under="ignore"):
                rate[positive] = self.electrons * np.exp(_log_sum_exp(np.stack(terms)))
        return (rate * AU_ANGULAR_FREQUENCY).reshape(strength.shape)

    def channel_rates(self, field_strength: ArrayLike) -> Float[np.ndarray, " c ..."]:
        """Return each channel's instantaneous rate; they sum to the rate.

        For a linearly polarised field: along an elliptical trajectory the
        instantaneous rate is a quasi-static approximation optmat does not validate.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            Rates per atom, 1/s, shape ``(len(channels),) + shape(field)``; 0 at zero
            field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        >>> argon = ionisation_rate_model("Ar", 800e-9)
        >>> argon.channel_rates([2e10, 5e10]).shape
        (2, 2)
        """
        strength, positive_au, positive = self._positive_fields(field_strength)
        out = np.zeros((len(self.channels), strength.size))
        if positive_au.size:
            with np.errstate(under="ignore"):
                out[:, positive] = np.exp(self._log_channel_rates_au(positive_au))
        return (out * AU_ANGULAR_FREQUENCY).reshape(
            (len(self.channels), *strength.shape)
        )

    def energy_loss_rate(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        """Return the energy the ionisation takes from the field, per atom and second.

        ``L = sum_c I_c Gamma_c``, the quantity a consumer averages over a cycle
        beside the rate (OEP 1, section 4.6). For a linearly polarised field.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            The energy-loss rate, W per atom, the shape of the field; 0 at zero field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.
        """
        energies = np.array([c.ionization_energy for c in self.channels])
        rates = self.channel_rates(field_strength)
        return np.tensordot(energies, rates, axes=1)

    def energy_per_event(
        self, field_strength: ArrayLike
    ) -> Float[np.ndarray, " *shape"]:
        """Return the mean ionisation energy of an event, ``L / Gamma``.

        The rate-weighted mean of the channels' energies, from the logarithms of
        their terms, so that it is defined where the rates underflow: the lowest
        channel's energy at zero field. For a linearly polarised field.

        Parameters
        ----------
        field_strength : array_like
            The field, V/m; its sign is ignored.

        Returns
        -------
        numpy.ndarray
            The energy, J, the shape of the field.

        Warns
        -----
        ModelRangeWarning
            If a field exceeds :attr:`maximum_field`.

        Examples
        --------
        >>> from optmat.constants import E_CHARGE
        >>> argon = ionisation_rate_model("Ar", 800e-9)
        >>> [f"{e / E_CHARGE:.3f}" for e in argon.energy_per_event([5e9, 5e10])]
        ['15.798', '15.814']
        """
        strength, positive_au, positive = self._positive_fields(field_strength)
        energies = np.array([c.ionization_energy for c in self.channels])
        out = np.full(strength.size, float(np.min(energies)))
        if positive_au.size:
            terms = self._log_enhanced_terms(positive_au)
            weights = np.exp(terms - _log_sum_exp(terms)[None, :])
            out[positive] = np.tensordot(energies, weights, axes=1)
        return out.reshape(strength.shape)

    @property
    def energy_kind(self) -> str:
        """The kind of the ionisation energies the exponent uses.

        Returns
        -------
        str
            ``adiabatic`` or ``vertical`` as the curated dataset states it (the two
            coincide for an atom); ``not stated`` for a custom species or an
            overridden energy.
        """
        overridden = {name for name, _ in self.modified} & {
            "channels",
            "ionization_energy",
        }
        if not isinstance(self.gas, Gas) or overridden:
            return "not stated"
        return ionisation_parameters(self.gas).ionization_energy_kind

    @classmethod
    def from_parameters(  # noqa: PLR0913 - the model's parameters
        cls,
        species: Gas | str,
        wavelength: float,
        *,
        ionization_energy: float,
        orbital_l: int,
        electrons: float,
        structure: Structure = "hartree",
        structure_coefficient: float | None = None,
        coulomb: CoulombFactor = "static",
        static_correction: StaticCorrection = "tong-lin",
        static_constants: tuple[float, ...] = (),
        channels: tuple[Channel, ...] | None = None,
        channel_sum: ChannelSum = "integral",
        averaging: Averaging = "consistent",
        uncertainty_factor: float = math.nan,
        citation: str = "",
        model: str = "user",
        sources: Mapping[str, str] | None = None,
        reasons: Mapping[str, str] | None = None,
    ) -> "IonisationRateModel":
        """Return a rate model from explicit parameters, checked as a curated one is.

        The curated models are built through it (:func:`ionisation_rate_model`).

        Parameters
        ----------
        species : Gas or str
            A curated gas, or a custom species' name, the identity's material.
        wavelength : float
            Vacuum wavelength of the field, m; ``math.inf`` for a static field.
        ionization_energy : float
            The lowest channel's ionisation energy, J. Given channels are moved to it,
            keeping their splittings.
        orbital_l : int
            The ionising orbital's l.
        electrons : float
            The effective number of electrons.
        structure : str, optional
            ``hartree`` (the default) or ``tong_2002``, with ``structure_coefficient``.
        structure_coefficient : float or None, optional
            Tong et al.'s C_l.
        coulomb : str, optional
            The Coulomb factor, ``static`` by default.
        static_correction : str, optional
            ``tong-lin`` (the default, with its alpha), ``zhang`` or ``none``.
        static_constants : tuple of float, optional
            The correction's constants.
        channels : tuple of Channel or None, optional
            The ion's levels and weights; one ground channel by default.
        channel_sum : str, optional
            ``integral`` or ``sum``.
        averaging : str, optional
            ``consistent`` or ``tunnelling``.
        uncertainty_factor : float, optional
            A stated uncertainty factor; NaN (none) by default.
        citation : str, optional
            Where the parameters come from.
        model : str, optional
            The model's name, ``"user"`` by default.
        sources : Mapping or None, optional
            The source of each constant.
        reasons : Mapping or None, optional
            Why each choice was made.

        Returns
        -------
        IonisationRateModel
            The model, status ``user``.

        Raises
        ------
        ValueError
            As the model's checks do.

        Examples
        --------
        >>> from optmat.constants import E_CHARGE
        >>> model = IonisationRateModel.from_parameters(
        ...     "Ar*", 800e-9, ionization_energy=15.76 * E_CHARGE, orbital_l=1,
        ...     electrons=2.0, static_constants=(9.0,), citation="a test")
        >>> model.identity.material, model.identity.status
        ('Ar*', 'user')
        """
        ground = Channel("ground level of the ion", ionization_energy, 1.0, citation)
        given = (ground,) if channels is None else channels
        return cls(
            gas=_species(species),
            wavelength=wavelength,
            model=model,
            channels=_moved(given, ionization_energy),
            orbital_l=orbital_l,
            structure=structure,
            structure_coefficient=structure_coefficient,
            electrons=electrons,
            coulomb=coulomb,
            static_correction=static_correction,
            static_constants=tuple(static_constants),
            channel_sum=channel_sum,
            averaging=averaging,
            uncertainty_factor=uncertainty_factor,
            sources=({"citation": citation} if citation else {}) | dict(sources or {}),
            reasons=dict(reasons or {}),
        )

    def _override(self, changes: Mapping[str, Any]) -> "IonisationRateModel":
        """Return the model with some inputs changed, through :meth:`from_parameters`.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values; ``ionization_energy`` moves every channel.

        Returns
        -------
        IonisationRateModel
            The rebuilt model.
        """
        inputs = {
            "orbital_l": self.orbital_l,
            "electrons": self.electrons,
            "structure": self.structure,
            "structure_coefficient": self.structure_coefficient,
            "coulomb": self.coulomb,
            "static_correction": self.static_correction,
            "static_constants": self.static_constants,
            "channels": self.channels,
            "channel_sum": self.channel_sum,
            "averaging": self.averaging,
            "uncertainty_factor": self.uncertainty_factor,
            "ionization_energy": min(c.ionization_energy for c in self.channels),
        } | dict(changes)
        return type(self).from_parameters(
            self.gas,
            self.wavelength,
            model=self.model,
            sources=self.sources,
            reasons=self.reasons,
            **inputs,
        )

    def _stated_uncertainty(self) -> Mapping[str, Any]:
        """Return the stated uncertainty an override drops.

        Returns
        -------
        Mapping
            ``uncertainty_factor`` to NaN.
        """
        return {"uncertainty_factor": math.nan}

    def at_wavelength(self, wavelength: float) -> "IonisationRateModel":
        """Return the same model at another wavelength.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength, m; ``math.inf`` for a static field.

        Returns
        -------
        IonisationRateModel
            The model at that wavelength.
        """
        return replace(self, wavelength=wavelength)


def _species(name: Gas | str) -> Gas | str:
    """Return a curated gas by its name, or a custom species' name as given.

    Parameters
    ----------
    name : Gas or str
        The species.

    Returns
    -------
    Gas or str
        The :class:`Gas` member it names, else the name.
    """
    try:
        return Gas(name)
    except ValueError:
        return str(name)


def _moved(channels: tuple[Channel, ...], lowest: float) -> tuple[Channel, ...]:
    """Return channels moved so the lowest has an energy, their splittings kept.

    Parameters
    ----------
    channels : tuple of Channel
        The channels.
    lowest : float
        The lowest channel's new ionisation energy, J.

    Returns
    -------
    tuple of Channel
        The channels, unchanged when the lowest already has it.
    """
    shift = lowest - min(c.ionization_energy for c in channels)
    if shift == 0:
        return channels
    return tuple(
        replace(c, ionization_energy=c.ionization_energy + shift) for c in channels
    )


def _value(by_key: Mapping[str, Dataset], reference: str) -> tuple[Dataset, int]:
    """Return a dataset and a row from ``source/dataset`` or ``source/dataset#row``.

    Parameters
    ----------
    by_key : mapping of str to Dataset
        The curated datasets.
    reference : str
        The reference.

    Returns
    -------
    tuple of (Dataset, int)
        The dataset and the index of the row.

    Raises
    ------
    DataError
        If the dataset or the row does not exist.
    """
    key, _, row = reference.partition("#")
    if key not in by_key:
        msg = f"ionisation_rates.toml names an unknown dataset {key!r}"
        raise DataError(msg)
    dataset = by_key[key]
    if not row:
        if len(dataset) != 1:
            msg = f"{key} has {len(dataset)} rows: name one as {key}#<row>"
            raise DataError(msg)
        return dataset, 0
    if row not in dataset.rows:
        msg = f"{key} has no row {row!r}"
        raise DataError(msg)
    return dataset, dataset.rows.index(row)


def _level_energy(dataset: Dataset, index: int) -> float:
    """Return the energy of an ion's level above its ground level, J.

    Parameters
    ----------
    dataset : Dataset
        A ``level_energy`` dataset.
    index : int
        The row.

    Returns
    -------
    float
        The energy, J.

    Raises
    ------
    DataError
        If the dataset is not a level energy in a unit the loader knows.
    """
    units = {"cm-1": 2 * math.pi * HBAR * C * 100, "eV": E_CHARGE}
    if dataset.quantity != "level_energy" or dataset.unit not in units:
        msg = f"{dataset.key} is not a level energy in cm-1 or eV"
        raise DataError(msg)
    return float(dataset.scaled_values()[index]) * units[dataset.unit]


def _channels(
    parameters: IonisationParameters,
    entry: Mapping[str, object],
    by_key: Mapping[str, Dataset],
) -> tuple[Channel, ...]:
    """Return the levels of the ion and their weights.

    Parameters
    ----------
    parameters : IonisationParameters
        The gas's curated parameters.
    entry : mapping
        The gas's table in a model of ``ionisation_rates.toml``.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    tuple of Channel
        One channel for an s shell or a model without fine structure; the two levels
        of the ion's ground term for a p shell.

    Raises
    ------
    DataError
        If the level is not of the gas's ion.
    """
    ground = Channel(
        label="ground level of the ion",
        ionization_energy=parameters.ionization_energy,
        weight=1.0,
        source=parameters.ionization_energy_source,
    )
    if "fine_structure_level" not in entry:
        return (ground,)
    reference = str(entry["fine_structure_level"])
    dataset, row = _value(by_key, reference)
    if dataset.material != f"{parameters.gas.value}+":
        msg = f"{dataset.key} is not a level of {parameters.gas.value}+"
        raise DataError(msg)
    # The m = 0 orbital of a p shell with either spin is sqrt(2/3) |J = 3/2> and
    # sqrt(1/3) |J = 1/2> of the hole it leaves (the Clebsch-Gordan coefficients of
    # l = 1 and s = 1/2 at m_j = 1/2), so the ion is left in 2P3/2 with weight 2/3.
    return (
        replace(ground, label="2P3/2", weight=2 / 3),
        Channel(
            label="2P1/2",
            ionization_energy=parameters.ionization_energy
            + _level_energy(dataset, row),
            weight=1 / 3,
            source=f"{parameters.ionization_energy_source} + {reference}",
        ),
    )


def _text(table: Mapping[str, object], name: str) -> dict[str, str]:
    """Return a sub-table of strings.

    Parameters
    ----------
    table : mapping
        A table of ``ionisation_rates.toml``.
    name : str
        The sub-table.

    Returns
    -------
    dict of str to str
        Its keys and values; empty when absent.
    """
    value = table.get(name, {})
    if not isinstance(value, Mapping):
        return {}
    return {str(k): str(v) for k, v in value.items()}  # pyright: ignore[reportUnknownVariableType]


def _build(
    gas: Gas,
    name: str,
    model: Mapping[str, object],
    by_key: Mapping[str, Dataset],
) -> IonisationRateModel:
    """Build a gas's model from its entry in ``ionisation_rates.toml``.

    Parameters
    ----------
    gas : Gas
        The gas.
    name : str
        The model's name.
    model : mapping
        The model's table.
    by_key : mapping of str to Dataset
        The curated datasets.

    Returns
    -------
    IonisationRateModel
        The model, at 800 nm.

    Raises
    ------
    DataError
        If the model does not cover the gas, or a constant it needs is not curated.
    """
    gases = model.get("gases", {})
    if not isinstance(gases, Mapping) or gas.value not in gases:
        msg = f"the ionisation-rate model {name!r} does not cover {gas.value}"
        raise DataError(msg)
    entry: Mapping[str, object] = gases[gas.value]  # pyright: ignore[reportUnknownVariableType]
    parameters = ionisation_parameters(gas)
    if parameters.orbital_l is None:
        msg = f"{gas.value} is not an atom: no rate model"
        raise DataError(msg)
    sources = {"ionization_energy": parameters.ionization_energy_source}
    structure = str(model["structure"])
    coefficient = None
    if structure != "hartree":
        tabulation = parameters.structure[structure]
        coefficient = float(tabulation.coefficients[parameters.orbital_l])
        sources["structure_coefficient"] = tabulation.source
    correction = str(model["static_correction"])
    constants: tuple[float, ...] = ()
    if correction == "tong-lin":
        if "tong_lin_alpha" in entry:
            constants = (float(str(entry["tong_lin_alpha"])),)
            sources["static_correction"] = f"ionisation_rates.toml ({name})"
        elif parameters.barrier_suppression_alpha is not None:
            constants = (parameters.barrier_suppression_alpha,)
            sources["static_correction"] = parameters.barrier_suppression_source
        else:
            msg = f"no Tong-Lin constant for {gas.value} in model {name!r}"
            raise DataError(msg)
    elif correction == "zhang":
        if parameters.barrier_suppression_zhang is None:
            msg = f"no Zhang constants for {gas.value}"
            raise DataError(msg)
        constants = parameters.barrier_suppression_zhang
        sources["static_correction"] = parameters.barrier_suppression_zhang_source
    channels = _channels(parameters, entry, by_key)
    if len(channels) > 1:
        sources["fine_structure_level"] = str(entry["fine_structure_level"])
    sources["electrons"] = f"ionisation_rates.toml ({name})"
    return IonisationRateModel.from_parameters(
        gas,
        800e-9,
        model=name,
        channels=channels,
        ionization_energy=min(c.ionization_energy for c in channels),
        orbital_l=parameters.orbital_l,
        structure=structure,  # pyright: ignore[reportArgumentType] - checked on construction
        structure_coefficient=coefficient,
        electrons=float(str(entry["electrons"])),
        coulomb=str(model["coulomb"]),  # pyright: ignore[reportArgumentType]
        static_correction=correction,  # pyright: ignore[reportArgumentType]
        static_constants=constants,
        channel_sum=str(model.get("channel_sum", "integral")),  # pyright: ignore[reportArgumentType]
        averaging=str(model.get("averaging", "consistent")),  # pyright: ignore[reportArgumentType]
        uncertainty_factor=float(str(entry.get("uncertainty_factor", math.nan))),
        sources=sources,
        reasons={**_text(model, "reasons"), **_text(entry, "reasons")},
    )


@functools.cache
def ionisation_rate_models() -> dict[str, dict[Gas, IonisationRateModel]]:
    """Return every named model of every gas it covers, at 800 nm.

    Returns
    -------
    dict of str to dict of Gas to IonisationRateModel
        The models by name and gas; :meth:`IonisationRateModel.at_wavelength` moves
        one to another wavelength.
    """
    with IONISATION_RATES_FILE.open("rb") as handle:
        table = tomllib.load(handle)
    by_key = _catalogue.by_key()
    built = {}
    for name, model in table["models"].items():
        # The default is validated (Kr and Xe with their suggested constant, which
        # Wallace et al.'s yields test); the literature forms are kept for comparison.
        status, reason = parse_status(
            model, str(name), default="validated" if name == "default" else "comparison"
        )
        built[str(name)] = {
            as_gas(gas): replace(
                _build(as_gas(gas), str(name), model, by_key),
                status=status,
                status_reason=reason,
            )
            for gas in model["gases"]
        }
    return built


def ionisation_rate_model(
    gas: Gas | str, wavelength: float, model: str = "default"
) -> IonisationRateModel:
    """Return the ionisation rate of a gas at a wavelength.

    Parameters
    ----------
    gas : Gas or str
        The gas: He, Ne, Ar, Kr or Xe.
    wavelength : float
        Vacuum wavelength of the field, m; ``math.inf`` for a static field.
    model : str, optional
        A model of ``data/gas/ionisation_rates.toml``. ``default`` is the validated
        one; the others reproduce forms of the literature, for comparison.

    Returns
    -------
    IonisationRateModel
        The model.

    Raises
    ------
    DataError
        If the model does not exist or does not cover the gas.

    Warns
    -----
    ModelRangeWarning
        If the wavelength is below the range over which the rates were compared.

    Examples
    --------
    The static rate of helium at 0.1 au, in atomic units (Parker et al.'s two-electron
    calculation gives 2.94e-6):

    >>> import math
    >>> from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD
    >>> model = ionisation_rate_model("He", math.inf)
    >>> rate = model.instantaneous_rate(0.1 * AU_FIELD) / AU_ANGULAR_FREQUENCY
    >>> f"{float(rate):.2e}"
    '2.79e-06'
    """
    models = ionisation_rate_models()
    if model not in models:
        msg = f"no ionisation-rate model {model!r}; the models are {sorted(models)}"
        raise DataError(msg)
    species = as_gas(gas)
    if species not in models[model]:
        msg = f"the ionisation-rate model {model!r} does not cover {species.value}"
        raise DataError(msg)
    return models[model][species].at_wavelength(wavelength)

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Franta's universal dispersion model (UDM) of a solid's dielectric function.

The UDM builds the complex dielectric function from contributions whose imaginary
parts are simple models of the joint density of states, each normalised to unit
transition strength (its Thomas-Reiche-Kuhn integral), and whose real parts are their
Kramers-Kronig transforms in closed form (Franta, Nečas, Ohlídal & Giglia, Proc. SPIE
9628, 96281U (2015), Eqs. (1)-(47); the curated record ``franta_2015b``):

.. math::

    \varepsilon(E) = 1 + N_{vc}\,\varepsilon^0_{vc} + N_{vx}\,\varepsilon^0_{vx}
        + N_{ut}\,\varepsilon^0_{ut} + \sum_p N_{loc,p}\,\varepsilon^0_{loc,p}
        + \sum_p N_{ph,p}\,\varepsilon^0_{ph,p},

with the interband transitions (a symmetric term and excitonic terms, Eqs. 4-27), the
higher-energy excitations of the valence electrons (Eqs. 28-31), the Urbach tail
(Eqs. 32-40), and Gaussian-broadened peaks for localised states and phonons
(Eqs. 41-47). Energies are photon energies in eV; the transition strengths N are in
eV^2. The model is causal and passive by construction, and every term is a
closed-form, piecewise-analytic function of E: it can be differentiated to any order
away from the few band edges where its joint density of states has a kink.

The parameters come from Franta et al.'s newAD2 parameter files, whose values the
source records hold as the terms of their formulas (TC-320). Two details of those
files follow the authors' program, which reproduces their published tables, rather
than the paper (TC-321): the Urbach tail's upper energy, and the wavenumber-to-energy
conversion of the phonon terms.

First derivatives are exact (complex-step differentiation of the closed forms).
"""

import functools
import itertools
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar, Final

import numpy as np
from jaxtyping import Complex, Float
from scipy import special

from optmat.exceptions import DataError

__all__ = [
    "UDM_WAVENUMBER_TO_EV",
    "Term",
    "UniversalDispersionModel",
    "from_file",
    "from_parameters",
    "parse_attributes",
    "parse_parameters",
    "read_parameter_file",
]

# The conversion of the phonon terms' wavenumbers and widths to energy (eV per cm^-1),
# as newAD2's program computes it: 1e2 h c / e with its constants h = 6.6260693e-34 J s
# and e = 1.60217653e-19 C (CODATA 2002; Franta's udm.cxx, the program published with
# the parameter files). The parameters were fitted in these units (TC-321).
UDM_WAVENUMBER_TO_EV: Final = 1e2 * 6.6260693e-34 * 299792458.0 / 1.60217653e-19
"""newAD2's conversion of a wavenumber in cm^-1 to a photon energy in eV."""

_FWHM_TO_RMS: Final = 1 / (2 * math.sqrt(2 * math.log(2)))

# Below a fraction of a term's lowest characteristic energy E_0 the closed-form real
# part loses digits to cancellation: its leading terms cancel, losing about
# 2 log10(E_0/E) digits for the interband, exciton and Urbach terms and 4 log10(E_0/E)
# for the higher-energy term. There the real part is the even series
# sum_k M_k E^2k of the term's moments (_moments), whose neglected E^24 term is below
# 1e-14 of it at the switch.
_SMALL: Final = {"interband": 0.05, "exciton": 0.05, "urbach": 0.05, "higher": 0.25}
_ORDERS: Final = 12

# Beyond this argument exp(-x) Ei(x) is evaluated by its asymptotic series (exp and Ei
# overflow or underflow separately); at |x| = 700 six terms reach 1e-16.
_ASYMPTOTIC: Final = 700.0

type Array = np.ndarray


def _real(x: Array) -> Array:
    """Return the real part (the value, under complex-step differentiation).

    Parameters
    ----------
    x : numpy.ndarray
        Real, or complex with an infinitesimal imaginary part.

    Returns
    -------
    numpy.ndarray
        Its real part.
    """
    return np.real(x)


def _abs(x: Array) -> Array:
    """Return |x| analytically in the complex step: x times the sign of Re x.

    Parameters
    ----------
    x : numpy.ndarray
        Real, or complex with an infinitesimal imaginary part.

    Returns
    -------
    numpy.ndarray
        The modulus |Re x|, carrying the derivative of |x| in its imaginary part.
    """
    return x * np.sign(_real(x))


def _log_abs(x: Array) -> Array:
    """Return ln|x|, analytically in the complex step.

    Parameters
    ----------
    x : numpy.ndarray
        Real, or complex with an infinitesimal imaginary part.

    Returns
    -------
    numpy.ndarray
        The logarithm ln|x|.
    """
    return np.log(_abs(x))


def _ei_scaled_real(x: Float[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
    """Return the scaled exponential integral exp(-x) Ei(x) of a real argument.

    Parameters
    ----------
    x : numpy.ndarray
        Real arguments, nonzero.

    Returns
    -------
    numpy.ndarray
        The scaled integral exp(-x) Ei(x), written Ei~ in Franta et al. 2015, Eq. (40).
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    near = np.abs(x) < _ASYMPTOTIC
    out[near] = np.exp(-x[near]) * special.expi(x[near])
    far = x[~near]
    # exp(-x) Ei(x) ~ (1/x) sum_k k! / x^k, for large |x| of either sign.
    term = 1 / far
    total = term.copy()
    for k in range(1, 8):
        term = term * k / far
        total += term
    out[~near] = total
    return out


def _ei_scaled(z: Array) -> Array:
    """Return exp(-z) Ei(z), carrying the complex-step derivative.

    Parameters
    ----------
    z : numpy.ndarray
        Real, or complex with an infinitesimal imaginary part.

    Returns
    -------
    numpy.ndarray
        The scaled integral exp(-z) Ei(z); for complex input, its value plus i Im(z)
        times its derivative 1/x - exp(-x) Ei(x).
    """
    x = np.asarray(_real(z), dtype=float)
    value = _ei_scaled_real(x.ravel()).reshape(x.shape)
    if not np.iscomplexobj(z):
        return value
    return value + 1j * np.imag(z) * (1 / x - value)


def _step(e: Array, low: float, high: float) -> Array:
    """Return the window Pi of Franta et al. 2015, Eq. (6): 1 for low < |E| < high.

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies.
    low, high : float
        The window.

    Returns
    -------
    numpy.ndarray
        1.0 inside, 0.0 outside.
    """
    a = np.abs(_real(e))
    return ((a > low) & (a < high)).astype(float)


def _off_edges(e: Array, edges: tuple[float, ...]) -> Array:
    """Move energies that sit on a band edge 1e-9 (relative) above it.

    At a term's edges (E_g, E_h, E_xi, ...) the closed forms' logarithmic
    singularities cancel analytically but not numerically; the function is continuous
    there, so the value 1e-9 away stands for it (the authors' program takes the limits
    instead).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies (complex under complex-step differentiation).
    edges : tuple of float
        The term's edges, eV.

    Returns
    -------
    numpy.ndarray
        The energies, moved where they sat within 1e-12 of an edge.
    """
    x = np.abs(_real(e))
    shift = np.zeros_like(x)
    for edge in edges:
        shift = np.where(np.abs(x - edge) <= 1e-12 * edge, 1e-9 * edge, shift)
    return e + np.sign(_real(e)) * shift if np.any(shift) else e


# --- The normalised contributions. Each returns (real part, imaginary part) for unit
# --- transition strength, as functions of E > 0 (odd imaginary, even real).


def _interband(e: Array, eg: float, eh: float) -> tuple[Array, Array]:
    """Return the symmetric interband term, Franta et al. 2015 Eqs. (4)-(12).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies, eV.
    eg, eh : float
        The band gap and the maximum transition energy, eV.

    Returns
    -------
    tuple of numpy.ndarray
        The real and imaginary parts of the normalised dielectric function.
    """
    e = _off_edges(e, (eg, eh))
    cn = (eh**2 - eg**2) * (eh**2 + eg**2 - 8 * eg * eh) / 12 + (
        eg * eh
    ) ** 2 * math.log(eh / eg)
    e2 = e * e
    ae = _abs(e)
    imag = (np.sign(_real(e)) * (ae - eg) ** 2 * (ae - eh) ** 2 / (cn * e2)) * _step(
        e, eg, eh
    )
    x = 2 * e * (eh * (eg**2 + e2) + eg * (eh**2 + e2))
    y = e2 * (eh**2 + eg**2 + 4 * eg * eh + e2) + (eg * eh) ** 2
    c = -2 * (eg * eh) ** 2 * math.log(eh / eg) - 3 * (eh**2 - eg**2) * e2
    real = (
        (y - x) * _log_abs((eh - e) / (eg - e))
        + (y + x) * _log_abs((eh + e) / (eg + e))
        + c
    ) / (math.pi * cn * e2)
    return real, imag


def _exciton(
    e: Array, eg: float, eh: float, ec: float, bc: float
) -> tuple[Array, Array]:
    """Return an excitonic interband term, Franta et al. 2015 Eqs. (13)-(25).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies, eV.
    eg, eh : float
        The band gap and the maximum transition energy, eV.
    ec, bc : float
        The Lorentzian's centre and half width, eV.

    Returns
    -------
    tuple of numpy.ndarray
        The real and imaginary parts of the normalised dielectric function.
    """
    e = _off_edges(e, (eg, eh))
    alpha = ec**2 + bc**2
    beta = eg + eh
    gamma = eg**2 + eh**2
    gh = eg * eh
    delta = alpha**2 * (beta - ec) - alpha * beta * gh + gh**2 * ec
    eta = alpha**2 - alpha * gamma + 4 * alpha * (beta * ec - ec**2 - gh) + gh**2
    arc = math.atan((ec - eg) / bc) + math.atan((eh - ec) / bc)
    log_l = math.log(((ec - eh) ** 2 + bc**2) / ((ec - eg) ** 2 + bc**2))
    log_hg = math.log(eh / eg)
    cn = (
        (eh - eg) * (4 * ec - 3 * beta) / 2
        + gh**2 / alpha * log_hg
        + (4 * delta * arc - eta * (bc * log_l + 2 * ec * arc)) / (2 * alpha * bc)
    )
    e2 = e * e
    ae = _abs(e)
    imag = (
        np.sign(_real(e))
        * (ae - eg) ** 2
        * (ae - eh) ** 2
        / (cn * e2 * ((ae - ec) ** 2 + bc**2))
    ) * _step(e, eg, eh)
    x = 2 * beta * e * (e2 + gh)
    y = e2 * e2 + (gamma + 4 * gh) * e2 + gh**2
    a = (y - x) / (e2 - 2 * ec * e + alpha)
    b = (y + x) / (e2 + 2 * ec * e + alpha)
    d = (
        delta * e2
        + alpha**2 * (alpha * beta - beta * gh - gamma * ec)
        + 4 * alpha * gh * ec * (beta * ec - alpha)
        + gh**2 * (2 * alpha * ec - 4 * ec**3)
    )
    ee = (
        eta * e2
        + alpha**3
        - alpha**2 * gamma
        - 4 * alpha**2 * gh
        + alpha * gh**2
        + 4 * alpha * beta * ec * gh
        - 4 * ec**2 * gh**2
    )
    c = (4 * arc * d - (bc * log_l + 2 * ec * arc) * ee) / (
        alpha * bc * (4 * ec**2 * e2 - (e2 + alpha) ** 2)
    ) * e2 - 2 * gh**2 / alpha * log_hg
    real = (
        a * _log_abs((eh - e) / (eg - e)) + b * _log_abs((eh + e) / (eg + e)) + c
    ) / (math.pi * cn * e2)
    return real, imag


def _higher(e: Array, exi: float) -> tuple[Array, Array]:
    """Return the higher-energy excitations, Franta et al. 2015 Eqs. (28)-(31).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies, eV.
    exi : float
        The onset E_xi = E_x + E_g / 2, eV.

    Returns
    -------
    tuple of numpy.ndarray
        The real and imaginary parts of the normalised dielectric function.
    """
    e = _off_edges(e, (exi,))
    e2 = e * e
    e3 = e2 * e
    ae = _abs(e)
    imag = 3 * exi * np.sign(_real(e)) * (ae - exi) ** 2 / (e3 * e2)
    imag = imag * _step(e, exi, math.inf)
    a = -((exi - e) ** 2) / e3
    b = (exi + e) ** 2 / e3
    real = (
        3
        * exi
        / (math.pi * e2)
        * (
            a * _log_abs(1 - e / exi)
            + b * _log_abs(1 + e / exi)
            - 2 / (3 * exi)
            - 2 * exi / e2
        )
    )
    return real, imag


def _urbach(e: Array, eu: float, eg: float, em: float) -> tuple[Array, Array]:
    """Return the Urbach (exponential) tail, Franta et al. 2015 Eqs. (32)-(40).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies, eV.
    eu : float
        The Urbach energy, eV.
    eg : float
        The band gap, eV.
    em : float
        The energy where the tail's middle part ends (E_m of Eq. 32), eV.

    Returns
    -------
    tuple of numpy.ndarray
        The real and imaginary parts of the normalised dielectric function.
    """
    e = _off_edges(e, (eg / 2, eg, em, em + eg / 2))
    half = eg / 2
    top = em + half
    floor = math.exp(-eg / (2 * eu))

    def ei_scaled(x: float) -> float:
        """Return exp(-x) Ei(x) of a scalar.

        Parameters
        ----------
        x : float
            The argument, nonzero.

        Returns
        -------
        float
            The scaled exponential integral.
        """
        return float(_ei_scaled_real(np.array([x]))[0])

    # Eq. (36), with exp(a) Ei(b) written through exp(-x) Ei(x) so that nothing
    # overflows when E_m / E_u is large.
    cn = (
        ei_scaled(eg / eu)
        - floor * ei_scaled(eg / (2 * eu))
        + floor * ei_scaled(-top / eu)
        - ei_scaled(-em / eu)
        + floor * math.log(half / top)
        + math.log(em / eg) * (1 - eg * em / (eu * (em - eg)))
        + (em + eg) / (2 * eu)
    )
    e2 = e * e
    ae = _abs(e)
    sign = np.sign(_real(e))
    # Each exponential is clipped outside its window, where it is not used.
    rise = (ae - eg) / eu
    rise = np.where(_real(rise) > 0, 0.0, rise)
    low = (np.exp(rise) - floor) * _step(e, half, eg)
    # The middle part includes E = E_g, where the pieces join (the program's rule).
    a_re = np.abs(_real(e))
    middle_window = ((a_re >= eg) & (a_re <= em)).astype(float)
    middle = (
        1
        + (em - eg) / (4 * eu)
        - (2 * ae - em - eg) ** 2 / (4 * eu * (em - eg))
        - floor
    ) * middle_window
    fall = (em - ae) / eu
    fall = np.where(_real(fall) > 0, 0.0, fall)
    high = (np.exp(fall) - floor) * _step(e, em, top)
    imag = sign * (low + middle + high) / (cn * e2)
    x = -e * (em + eg) / (eu * (em - eg))
    y = 1 - (em * eg + e2) / (eu * (em - eg))
    c = 2 * (
        (em * eg / (eu * (em - eg)) - 1) * math.log(em / eg)
        - ei_scaled(eg / eu)
        + ei_scaled(-em / eu)
        + floor
        * (ei_scaled(half / eu) - ei_scaled(-top / eu) + math.log(em / half + 1))
    )
    real = (
        (y - x) * _log_abs((em - e) / (eg - e))
        + _ei_scaled((eg - e) / eu)
        - _ei_scaled(-(em - e) / eu)
        + (y + x) * _log_abs((em + e) / (eg + e))
        + _ei_scaled((eg + e) / eu)
        - _ei_scaled(-(em + e) / eu)
        + c
        - floor
        * (
            _log_abs(em / (half - e) + 1)
            + _ei_scaled((half - e) / eu)
            - _ei_scaled((e - top) / eu)
        )
        - floor
        * (
            _log_abs(em / (half + e) + 1)
            + _ei_scaled((half + e) / eu)
            - _ei_scaled(-(e + top) / eu)
        )
    ) / (math.pi * cn * e2)
    return real, imag


def _gaussian(e: Array, ep: float, bp: float) -> tuple[Array, Array]:
    """Return a Gaussian-broadened peak, Franta et al. 2015 Eqs. (41)-(43).

    Parameters
    ----------
    e : numpy.ndarray
        Photon energies, eV.
    ep : float
        The peak energy, eV.
    bp : float
        The rms broadening, eV.

    Returns
    -------
    tuple of numpy.ndarray
        The real and imaginary parts of the normalised dielectric function.
    """
    s = math.sqrt(2) * bp
    um = (e - ep) / s
    up = (e + ep) / s
    imag = (np.exp(-um * um) - np.exp(-up * up)) / (math.sqrt(math.pi) * s * ep)
    real = 2 * (special.dawsn(up) - special.dawsn(um)) / (math.pi * s * ep)
    return real, imag


# --- The model.


@dataclass(frozen=True)
class Term:
    """One contribution of Eq. (1).

    Attributes
    ----------
    kind : str
        "interband", "exciton", "higher", "urbach" or "gaussian".
    strength : float
        Its transition strength N, eV^2.
    shape : tuple of float
        Its shape parameters, eV: (E_g, E_h); (E_g, E_h, E_c, B_c); (E_xi,);
        (E_u, E_g, E_h); or (centre, rms width).
    lowest : float
        Its lowest characteristic energy, below which the real part is the moment
        series (0 for a Gaussian peak).
    """

    kind: str
    strength: float
    shape: tuple[float, ...]
    lowest: float = field(default=0.0)


_SHAPES: Final = {
    "interband": _interband,
    "exciton": _exciton,
    "higher": _higher,
    "urbach": _urbach,
    "gaussian": _gaussian,
}


def parse_attributes(line: str, where: str) -> tuple[str, ...]:
    """Return the attributes of a newAD2 parameter file's first line.

    Parameters
    ----------
    line : str
        The line, such as ``he:ex=10:ut:ph=52``.
    where : str
        What the line belongs to, for the error message.

    Returns
    -------
    tuple of str
        The attributes.

    Raises
    ------
    DataError
        If an attribute is one this implementation does not support.

    Examples
    --------
    >>> parse_attributes("he:ex=10:ut:ph=52", "Lithosil Q2")
    ('he', 'ex=10', 'ut', 'ph=52')
    """
    attributes = tuple(a for a in line.strip().split(":") if a)
    for attribute in attributes:
        name = attribute.split("=")[0]
        if name not in {"he", "ut", "ex", "loc", "ph"}:
            msg = f"{where}: UDM attribute {attribute!r} is not supported"
            raise DataError(msg)
    return attributes


def parse_parameters(text: str, where: str) -> tuple[tuple[str, ...], dict[str, float]]:
    r"""Read the text of a newAD2 UDM parameter file.

    Parameters
    ----------
    text : str
        The file's text: a first line of attributes (``he:ex=10:ut:ph=52``), then
        ``name = value`` lines, then the site's free-text description.
    where : str
        What the text belongs to, for the error message.

    Returns
    -------
    tuple
        The attributes, and the parameters by name.

    Examples
    --------
    >>> attributes, values = parse_parameters("he:ut\nEg = 8.36\nEh = 39.6\n", "x")
    >>> attributes, values["Eg"]
    (('he', 'ut'), 8.36)
    """
    lines = text.splitlines()
    attributes = parse_attributes(lines[0], where)
    pattern = re.compile(r"^\s*(\w+)\s*=\s*([-+0-9.eE]+)\s*$")
    values = {
        match.group(1): float(match.group(2))
        for line in lines[1:]
        if (match := pattern.match(line))
    }
    return attributes, values


def read_parameter_file(path: Path) -> tuple[tuple[str, ...], dict[str, float]]:
    """Read a newAD2 UDM parameter file.

    Parameters
    ----------
    path : pathlib.Path
        The ``.par`` file (:func:`parse_parameters`).

    Returns
    -------
    tuple
        The attributes, and the parameters by name.
    """
    return parse_parameters(path.read_text(encoding="utf-8"), path.name)


def _count(attributes: tuple[str, ...], name: str) -> int:
    """Return the n of an ``name=n`` attribute (0 when absent).

    Parameters
    ----------
    attributes : tuple of str
        The attributes.
    name : str
        The counted attribute.

    Returns
    -------
    int
        The count.
    """
    for attribute in attributes:
        key, _, count = attribute.partition("=")
        if key == name:
            return int(count)
    return 0


@dataclass(frozen=True)
class UniversalDispersionModel:
    """A dielectric function in Franta et al.'s universal dispersion model.

    Build one from a newAD2 parameter file with :func:`from_file`.

    Attributes
    ----------
    terms : tuple
        The contributions, each with its kind, transition strength (eV^2) and shape
        parameters (eV).
    """

    terms: tuple[Term, ...]

    EVALUATION_VERSION: ClassVar[int] = 1

    def _parts(self, e: Array) -> tuple[Array, Array]:
        """Sum the terms' real and imaginary parts (the real part without the 1).

        Parameters
        ----------
        e : numpy.ndarray
            Photon energies, eV (complex under complex-step differentiation).

        Returns
        -------
        tuple of numpy.ndarray
            The real part less one, and the imaginary part.
        """
        real = np.zeros_like(e)
        imag = np.zeros_like(e)
        for term in self.terms:
            shape = _SHAPES[term.kind]
            part_re, part_im = shape(e, *term.shape)
            small = np.abs(_real(e)) < _SMALL.get(term.kind, 0.0) * term.lowest
            if np.any(small):
                series = sum(m * (e * e) ** k for k, m in enumerate(_moments(term)))
                part_re = np.where(small, series, part_re)
            real = real + term.strength * part_re
            imag = imag + term.strength * part_im
        return real, imag

    def permittivity(
        self, energy: Float[np.ndarray, " n"]
    ) -> Complex[np.ndarray, " n"]:
        """Return the complex dielectric function at photon energies.

        Parameters
        ----------
        energy : numpy.ndarray
            Photon energies, eV, > 0.

        Returns
        -------
        numpy.ndarray
            The permittivity epsilon_r + i epsilon_i, with epsilon_i >= 0 for
            absorption.

        Examples
        --------
        >>> from optmat.solid.models.formulas import published_models
        >>> model = published_models()["franta_2016/udm_lithosil_q2"].udm
        >>> f"{model.transition_strength():.1f}"
        '1116.0'
        >>> n = np.sqrt(model.permittivity(np.array([1.5498])))
        >>> f"{n.real[0]:.5f}"
        '1.45355'
        """
        e = np.asarray(energy, dtype=float)
        real, imag = self._parts(e)
        return (1 + real) + 1j * imag

    def dpermittivity_denergy(
        self, energy: Float[np.ndarray, " n"]
    ) -> Complex[np.ndarray, " n"]:
        """Return d epsilon / dE, exact to rounding (complex-step differentiation).

        Parameters
        ----------
        energy : numpy.ndarray
            Photon energies, eV, > 0, away from the terms' band edges.

        Returns
        -------
        numpy.ndarray
            The derivative, per eV.
        """
        e = np.asarray(energy, dtype=float)
        h = 1e-20 * np.maximum(np.abs(e), 1e-3)
        real, imag = self._parts(e + 1j * h)
        return (np.imag(real) + 1j * np.imag(imag)) / h

    def transition_strength(self) -> float:
        """Return the total transition strength, the sum of the terms' N (eV^2).

        Returns
        -------
        float
            The integral of epsilon_i(E) E over E (Franta et al. 2015, Eq. 3).
        """
        return sum(term.strength for term in self.terms)


def from_parameters(
    attributes: tuple[str, ...], values: Mapping[str, float]
) -> UniversalDispersionModel:
    """Build a model from a parameter file's attributes and values.

    Parameters
    ----------
    attributes : tuple of str
        The model's attributes (``he``, ``ut``, ``ex=n``, ``loc=n``, ``ph=n``).
    values : mapping of str to float
        The parameters, in newAD2's names and units (eV, eV^2, cm^-1).

    Returns
    -------
    UniversalDispersionModel
        The model.

    Raises
    ------
    DataError
        If a parameter the attributes call for is missing.
    """
    try:
        return UniversalDispersionModel(_terms(attributes, values))
    except KeyError as exc:
        msg = f"UDM parameter {exc.args[0]} is missing"
        raise DataError(msg) from exc


def from_file(path: Path) -> UniversalDispersionModel:
    """Build a model from a newAD2 parameter file.

    Parameters
    ----------
    path : pathlib.Path
        The ``.par`` file.

    Returns
    -------
    UniversalDispersionModel
        The model.
    """
    return from_parameters(*read_parameter_file(path))


def _terms(
    attributes: tuple[str, ...], values: Mapping[str, float]
) -> tuple[Term, ...]:
    """Build the terms of Franta et al. 2015 Eq. (1) from newAD2 parameters.

    Parameters
    ----------
    attributes : tuple of str
        The model's attributes.
    values : mapping of str to float
        The parameters.

    Returns
    -------
    tuple of Term
        The contributions.
    """
    eg, eh = values["Eg"], values["Eh"]
    terms: list[Term] = []
    for p in range(1, _count(attributes, "ph") + 1):
        energy = UDM_WAVENUMBER_TO_EV * values[f"nuph{p}"]
        width = UDM_WAVENUMBER_TO_EV * _FWHM_TO_RMS * values[f"betaph{p}"]
        terms.append(Term("gaussian", values[f"Nph{p}"], (energy, width)))
    if "he" in attributes:
        exi = values["Ex"] + eg / 2
        terms.append(Term("higher", values["Nvx"], (exi,), exi))
    if "ut" in attributes:
        # newAD2's program ends the tail's middle part at E_h, where Eq. (33) of the
        # paper defines E_m = (E_g + E_h)/2; the program reproduces the published
        # tables, so its definition is the one the parameters were fitted with (TC-321).
        terms.append(Term("urbach", values["Nut"], (values["Eu"], eg, eh), eg / 2))
    for p in range(1, _count(attributes, "loc") + 1):
        energy = values[f"Eloc{p}"]
        terms.append(Term("gaussian", values[f"Nloc{p}"], (energy, values[f"Bloc{p}"])))
    # Eqs. (26)-(27): the interband strength N_vc shared in proportion to the A's.
    excitons = _count(attributes, "ex")
    a0 = values.get("A0", 1.0)
    amplitudes = [values[f"Aex{p}"] for p in range(1, excitons + 1)]
    total = a0 + sum(amplitudes)
    terms.append(Term("interband", values["Nvc"] * a0 / total, (eg, eh), eg))
    for p, amplitude in enumerate(amplitudes, start=1):
        shape = (eg, eh, values[f"Eex{p}"], values[f"Bex{p}"])
        terms.append(Term("exciton", values["Nvc"] * amplitude / total, shape, eg))
    return tuple(terms)


def _edges(term: Term) -> list[float]:
    """Return the energies where a term's imaginary part has a kink or a peak.

    Parameters
    ----------
    term : Term
        An interband, exciton, higher-energy or Urbach term.

    Returns
    -------
    list of float
        The energies, eV, for the quadrature's breakpoints.
    """
    match term.kind:
        case "urbach":
            _, eg, em = term.shape
            return [eg / 2, eg, em, em + eg / 2]
        case "exciton":
            eg, eh, ec, bc = term.shape
            return sorted({eg, eh} | {x for x in (ec - bc, ec, ec + bc) if eg < x < eh})
        case "interband":
            return list(term.shape)
        case _:
            return [term.shape[0]]


# The moments' quadrature: each stretch between a term's edges is cut into pieces of
# equal logarithmic width, each integrated by a 24-point Gauss-Legendre rule; checked
# against adaptive quadrature to 1e-12 (tests/solid/test_udm.py).
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(24)
_PIECES_PER_STRETCH: Final = 24


def _quadrature_grid(
    edges: list[float], start: float, top: float
) -> tuple[Array, Array]:
    """Return composite Gauss-Legendre nodes and weights from ``start`` to ``top``.

    Parameters
    ----------
    edges : list of float
        Breakpoints inside the range, eV.
    start, top : float
        The range, eV.

    Returns
    -------
    tuple of numpy.ndarray
        The nodes and their weights.
    """
    cuts = sorted({start, top, *(x for x in edges if start < x < top)})
    bounds = np.concatenate(
        [
            np.geomspace(a, b, _PIECES_PER_STRETCH + 1)[:-1]
            for a, b in itertools.pairwise(cuts)
        ]
        + [np.array([top])]
    )
    low, high = bounds[:-1, None], bounds[1:, None]
    nodes = (low + high) / 2 + (high - low) / 2 * _NODES
    weights = (high - low) / 2 * _WEIGHTS
    return nodes.ravel(), (weights * np.ones_like(nodes)).ravel()


@functools.cache
def _moments(term: Term) -> tuple[float, ...]:
    """Return the low-energy series of a term's real part, eps_r = sum_k M_k E^2k.

    Below the absorption the Kramers-Kronig integral
    eps_r(E) = (2/pi) int eps_i(x) x / (x^2 - E^2) dx expands as
    M_k = (2/pi) int eps_i(x) / x^(2k+1) dx.

    Parameters
    ----------
    term : Term
        The term (normalised: its strength is applied by the caller).

    Returns
    -------
    tuple of float
        M_0 ... M_11, in eV^(-2k).
    """
    start = term.lowest / 2
    edges = _edges(term)
    top = 50 * max(edges)
    if term.kind == "exciton":
        # Refine around the Lorentzian, whose width can be far below its energy.
        _, _, ec, bc = term.shape
        edges = edges + [ec + s * bc for s in (-8, -4, -2, -1, 1, 2, 4, 8)]
    nodes, weights = _quadrature_grid(edges, start, top)
    imag = _SHAPES[term.kind](nodes, *term.shape)[1]
    moments = []
    for k in range(_ORDERS):
        power = 2 * k + 1
        value = float(np.sum(weights * imag / nodes**power))
        if term.kind == "higher":
            # Beyond `top` eps_i = 3 E_xi (x - E_xi)^2 / x^5 (Eq. 29): its integral
            # against x^-power from `top` to infinity, exactly.
            exi = term.shape[0]
            q = power + 2
            value += (
                3
                * exi
                * (
                    1 / (q * top**q)
                    - 2 * exi / ((q + 1) * top ** (q + 1))
                    + exi**2 / ((q + 2) * top ** (q + 2))
                )
            )
        moments.append(2 / math.pi * value)
    return tuple(moments)

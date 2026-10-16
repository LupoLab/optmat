# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The Raman response function of a molecular gas or a gas mixture (TC-685 to TC-690).

A molecule's delayed third-order response is a sum of damped oscillators, one per Raman
line (:mod:`optmat.gas.raman`). This module assembles them into the response function
itself, with its polarisation (tensor) structure and its dependence on the density,
so that a propagation code only samples it.

**The convention** (TC-442, TC-685). For a field linearly polarised along :math:`x`,
the Raman polarisation of a gas of number density :math:`N` is

.. math::

   P_x(t) = \varepsilon_0 N E_x(t) \int_0^\infty h(\tau) E_x^2(t - \tau)\, d\tau,
   \qquad h(\tau) = \sum_k A_k \sin(\Omega_k \tau)\, e^{-\Gamma_k \tau},

with :math:`A_k` in m^5 V^-2 s^-1 per molecule of the gas (of the mixture: a
constituent's lines carry its mole fraction), :math:`\Omega_k` in rad/s and
:math:`\Gamma_k = \pi \Delta\nu_k`, :math:`\Delta\nu_k` the line's full width at half
maximum in Hz. Its transform, with the sign convention of
:func:`optmat.gas.raman.stokes_gain`, is

.. math::

   H(\omega) = \int_0^\infty h(\tau) e^{i\omega\tau} d\tau
   = \sum_k \frac{A_k \Omega_k}
                  {\Omega_k^2 + \Gamma_k^2 - \omega^2 - 2 i \Gamma_k \omega},

so that :math:`H(\Omega_k) \approx i A_k / (2\Gamma_k)` at a resonance and
:math:`H(0) = \int h\, d\tau`.

**The tensor form** (Hellwarth). For any polarisation,

.. math::

   P_i(t) = \varepsilon_0 N \int_0^\infty R_{ijkl}(\tau)\,
       E_j(t) E_k(t - \tau) E_l(t - \tau)\, d\tau, \qquad
   R_{ijkl} = a(\tau)\, \delta_{ij}\delta_{kl}
       + \tfrac{1}{2} b(\tau) (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}),

and :math:`h = a + b` for a single linear polarisation. Each line has an isotropic
(rank-0) share, from the mean polarizability :math:`\bar\alpha`, and an anisotropic
(rank-2) share, from its anisotropy. The isotropic share enters :math:`a` only; the
anisotropic share :math:`h_2` has :math:`a = -h_2/2`, :math:`b = 3h_2/2` (the induced
change of the anisotropic polarizability is the traceless part of :math:`E E`). A weak
probe polarised across the pump sees :math:`a`: the rotational lines give it
:math:`-1/2` of the parallel response, the isotropic vibrational lines all of it.

- The rotational S lines are anisotropic (:math:`\Delta\alpha`).
- A vibrational Q(J) line has both shares. Its anisotropic share is
  :math:`(4/45)\, b_{J,J}\, \gamma'^2 P_J / (\hbar\varepsilon_0)`, with the
  Placzek-Teller coefficient :math:`b_{J,J} = J(J+1)/[(2J-1)(2J+3)]` and
  :math:`\gamma'` the transition anisotropy: from its components where they are
  curated (H2, D2), else from the Q branch's depolarization ratio
  :math:`\rho = 3 b\gamma'^2 / (45\bar\alpha'^2 + 4 b\gamma'^2)` (TC-686).
- The vibrational O and S branches (:math:`\Delta J = \mp 2` with :math:`v = 0 \to 1`),
  anisotropic only, are left out: they carry :math:`3\times` the anisotropic Q share,
  8 % of the Q branch for N2, 19 % for O2 and 3 % for H2 and D2, spread over some
  100 cm^-1 about the band
  (TC-686).

**The widths and the density.** Every curated width law has the form
:math:`\Delta\nu = c + d/n + \beta n` at a fixed temperature and composition (n in
amagat; :mod:`optmat.gas.linewidths`), so the response keeps the three coefficients
of each line (TC-687). A line's frequency shifts with the density where a shift is
curated (H2's Q branch; TC-688). The amplitudes depend on the temperature and the pump
wavelength only. The response is therefore fixed by (temperature, composition,
wavelength), and :meth:`GasRamanResponse.decay_rates` and
:meth:`GasRamanResponse.frequencies` give its density dependence cheaply, for any array
of densities: the part a propagation code evaluates along a density gradient. Lines
whose widths follow one law have the same ``damping_group``: they can be summed into
one damped kernel without approximation.

**Ranges** (``ModelRangeWarning``, TC-688). The lines are independent Lorentzians
with collisional (impact) widths, which holds between two densities, found for the
strongest line of each band:

- below, the motion of the molecules sets the width. For a line with a Dicke term the
  diffusion law holds above Herring et al.'s cutoff :math:`\rho_c = 3.33\, d /
  \Delta\nu_D` (their Eq. (3)), with :math:`\Delta\nu_D` the Doppler width for
  forward scattering (their Eq. (4)); below it :math:`d/n` is capped at
  :math:`\Delta\nu_D`, a Lorentzian stand-in for the Doppler limit. For a line without
  one, the limit is the density at which its collisional width falls to
  :math:`\Delta\nu_D`;
- above, the line's width reaches the spacing to its nearest neighbour: the overlap and
  mixing of the lines (the collisional narrowing of a Q branch) is not modelled
  (TC-527);
- the temperatures of each band's width measurements (warned when the response is
  built).

**What is not here.** The sampling of the response on a time grid, its taper and
padding, and any grouping or low-rank tabulation for speed belong to the propagation
code (the manual page *The Raman response*). The response is a coherent (Kubo)
response: the redistribution of the M sublevels within a level, which a slow
isothermal measurement also sees, is not part of it (TC-493).

References
----------
- R. W. Hellwarth, "Third-order optical susceptibilities of liquids and solids", Prog.
  Quantum Electron. 5, 1 (1977).
- G. Placzek and E. Teller, "Die Rotationsstruktur der Ramanbanden mehratomiger
  Moleküle", Z. Phys. 81, 209 (1933).
- D. A. Long, *The Raman Effect* (Wiley, 2002), chapters 5 and 6.
"""

import functools
import math
import tomllib
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any, ClassVar, Final, Literal, get_args

import numpy as np
from jaxtyping import Complex, Float, Int
from numpy.typing import ArrayLike

from optmat._catalogue import records as sources_of
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.composition import composition
from optmat.gas.gases import Gas, as_gas, info
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import (
    AnisotropyTable,
    BandLines,
    RamanLines,
    RotationalBand,
    VibrationalBand,
    band_lines,
    rotational_nonlinear_index,
    stokes_gain,
)
from optmat.identity import Identified, Status, parse_status, source_ids

__all__ = [
    "EVALUATION_VERSION",
    "RAMAN_RESPONSE_FILE",
    "Component",
    "GasRamanResponse",
    "RamanGroup",
    "placzek_teller_q",
    "raman_response",
    "raman_response_status",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

RAMAN_RESPONSE_FILE: Final = GAS_DATA_DIRECTORY / "raman_response.toml"
"""The range, shift and mixture choices of the responses, each with its reason."""

Component = Literal["parallel", "isotropic", "anisotropic", "a", "b"]
"""A part of the response.

``parallel``: :math:`h = a + b`, the response to a single linear polarisation.
``isotropic`` and ``anisotropic``: the rank-0 and rank-2 shares of ``parallel``.
``a`` and ``b``: Hellwarth's two functions; ``a`` is the response a weak probe
polarised across the pump sees.
"""

_PROBE_DENSITIES: Final = (1.0, 2.0, 4.0)
"""Densities, amagat, at which a width law is evaluated to find its coefficients.

``c + d/n + beta n`` at 1, 2 and 4 amagat is a well-conditioned 3 x 3 system (its
condition number is 20).
"""

_CHECK_DENSITY: Final = 8.0
"""A fourth density, amagat, at which the coefficients must reproduce the law."""

_LAW_TOLERANCE: Final = 1e-9
"""Relative residual at the fourth density above which a law is refused.

The solve loses about the condition number (20) times the double-precision epsilon;
a law of another form (a power of the density, say) misses by far more.
"""

_TERM_BLOCK: Final = 32
"""Lines summed per block in :meth:`GasRamanResponse.response`, bounding the memory."""


def placzek_teller_q(j: Int[np.ndarray, " n"]) -> Float[np.ndarray, " n"]:
    r"""Return the Placzek-Teller coefficient of the Q lines, :math:`b_{J,J}`.

    The share of a level's anisotropic scattering that keeps :math:`J`,
    :math:`J(J+1)/[(2J-1)(2J+3)]`; with :math:`b_{J+2,J}` and :math:`b_{J-2,J}` the
    three sum to one. It is 0 for :math:`J = 0` and tends to 1/4.

    Parameters
    ----------
    j : numpy.ndarray
        Rotational levels, >= 0.

    Returns
    -------
    numpy.ndarray
        :math:`b_{J,J}`.

    Examples
    --------
    >>> [f"{b:.4f}" for b in placzek_teller_q(np.array([0, 1, 2, 50]))]
    ['0.0000', '0.4000', '0.2857', '0.2501']
    """
    x = np.asarray(j, dtype=float)
    # J = 0 has no anisotropic Q line: the product is 0 / -3, kept at +0.
    return np.abs(x * (x + 1) / ((2 * x - 1) * (2 * x + 3)))


_COMPONENT_WEIGHTS: Final = {
    "parallel": (1.0, 1.0),
    "isotropic": (1.0, 0.0),
    "anisotropic": (0.0, 1.0),
    "a": (1.0, -0.5),
    "b": (0.0, 1.5),
}
"""The isotropic and anisotropic weights of each component (the module docstring)."""


@dataclass(frozen=True)
class RamanGroup:
    r"""One damping group of a Raman response: lines that share a width law and shift.

    The group's decay rate at a total number density :math:`N` (m^-3) is
    :math:`\Gamma_g(N) = \gamma_c + \min(\gamma_d/N,\ \gamma_D) + \gamma_\beta N`, 1/s,
    and every line's angular frequency :math:`\omega_k + \sigma_g N`, so that the
    group's response is :math:`e^{-\Gamma_g t}\sum_k A_k \sin((\omega_k + \sigma_g N)
    t)` (OEP 1, section 4.5). Where :math:`\sigma_g = 0` the row of amplitudes and
    frequencies is fixed and only the decay depends on the density.

    Attributes
    ----------
    index : int
        The group's number (``damping_group``).
    lines : tuple of int
        The indices of its lines in the response.
    line_ids : tuple of str
        Their identifiers (:attr:`GasRamanResponse.line_id`).
    gamma_constant : float
        :math:`\gamma_c = \pi c`, 1/s.
    gamma_dicke : float
        :math:`\gamma_d = \pi d N_\text{amg}`, 1/s m^-3; 0 without a Dicke term.
    gamma_doppler : float
        :math:`\gamma_D = \pi \Delta\nu_D`, the Dicke term's cap, 1/s; 0 without one.
    gamma_broadening : float
        :math:`\gamma_\beta = \pi \beta / N_\text{amg}`, 1/s m^3.
    shift : float
        :math:`\sigma_g = 2\pi s_g / N_\text{amg}`, rad s^-1 m^3.
    frequency : numpy.ndarray
        The lines' angular frequencies at zero density, rad/s.
    isotropic, anisotropic : numpy.ndarray
        Their amplitudes' two shares, m^5 V^-2 s^-1.
    """

    index: int
    lines: tuple[int, ...]
    line_ids: tuple[str, ...]
    gamma_constant: float
    gamma_dicke: float
    gamma_doppler: float
    gamma_broadening: float
    shift: float
    frequency: Float[np.ndarray, " m"]
    isotropic: Float[np.ndarray, " m"]
    anisotropic: Float[np.ndarray, " m"]

    @property
    def shifted(self) -> bool:
        """Whether the group's frequencies move with the density.

        Returns
        -------
        bool
            True where :attr:`shift` is not 0.
        """
        return self.shift != 0

    def decay_rate(self, density: ArrayLike) -> Float[np.ndarray, " *shape"]:
        """Return the group's decay rate, down to zero density.

        At ``N = 0`` it is the law's limit: the capped Dicke term tends to its cap.
        It does not warn outside the response's validated range; see
        :meth:`GasRamanResponse.check_densities`.

        Parameters
        ----------
        density : array_like
            Total number densities, m^-3, >= 0.

        Returns
        -------
        numpy.ndarray
            The decay rate, 1/s, the shape of ``density``.

        Raises
        ------
        DataError
            If a density is negative or not finite.
        """
        n = np.asarray(density, dtype=float)
        if not np.all(np.isfinite(n) & (n >= 0)):
            msg = "a decay rate needs finite densities >= 0"
            raise DataError(msg)
        with np.errstate(divide="ignore"):
            dicke = (
                np.minimum(self.gamma_dicke / n, self.gamma_doppler)
                if self.gamma_dicke > 0
                else np.zeros(n.shape)
            )
        return np.asarray(self.gamma_constant + dicke + self.gamma_broadening * n)

    def response(
        self, time: ArrayLike, density: float, component: Component = "parallel"
    ) -> Float[np.ndarray, " *shape"]:
        """Return the group's share of the response at a density.

        Parameters
        ----------
        time : array_like
            Delay after the impulse, s (zero before it).
        density : float
            Total number density, m^-3, >= 0.
        component : {"parallel", "isotropic", "anisotropic", "a", "b"}, optional
            The part of the response.

        Returns
        -------
        numpy.ndarray
            The group's response, m^5 V^-2 s^-1, the shape of ``time``.
        """
        tau = np.asarray(time, dtype=float)
        positive = np.maximum(tau, 0.0)[..., None]
        iso, aniso = _COMPONENT_WEIGHTS[component]
        amplitude = iso * self.isotropic + aniso * self.anisotropic
        omega = self.frequency + self.shift * density
        rate = float(self.decay_rate(density))
        total = np.exp(-rate * positive[..., 0]) * np.sum(
            amplitude * np.sin(omega * positive), axis=-1
        )
        return np.where(tau >= 0, total, 0.0)


@dataclass(frozen=True)
class GasRamanResponse(Identified):
    r"""The Raman response of a gas or a mixture at one temperature and composition.

    Built by :func:`raman_response`, or by :meth:`from_lines` from a caller's lines.
    Every array has one entry per line (a term of the response); the lines are
    grouped by species and band, the rotational S lines first. The damping groups and
    the density range are derived from the lines by one rule (:meth:`from_lines`).

    Its overrides (:meth:`~optmat.identity.Identified.with_parameters`) act on the
    assembled lines: ``amplitude_scale``, a factor on every amplitude;
    ``width_scale``, a factor on ``c``, ``d``, ``beta`` and the Doppler cap alike,
    which keeps the groups; and ``width_law``, ``{band: (c, d, beta, cap)}`` for
    every line of a band (``"rotation"``, a mode, or ``"<species> <band>"`` in a
    mixture). A different temperature, composition, pump wavelength or population
    floor is a new :func:`raman_response`, not an override.

    Attributes
    ----------
    composition : tuple of tuple of (str, float)
        The species and their mole fractions (a pure gas is one species at 1).
    temperature : float
        Temperature, K.
    wavelength : float or None
        The pump wavelength at which the vibrational strengths are dispersed, m;
        None for the strengths as curated.
    species : tuple of str
        The species of each line.
    band : tuple of str
        The band of each line (``rotation`` or the vibrational mode).
    branch : tuple of str
        ``S`` or ``Q``.
    j : numpy.ndarray
        The lower rotational level (-1 for an unresolved band).
    frequency : numpy.ndarray
        :math:`\Omega` at zero density, rad/s.
    isotropic, anisotropic : numpy.ndarray
        The rank-0 and rank-2 shares of each line's amplitude for a single linear
        polarisation, per molecule of the gas, m^5 V^-2 s^-1.
    width_constant : numpy.ndarray
        The density-independent width :math:`c`, Hz.
    width_dicke : numpy.ndarray
        The Dicke coefficient :math:`d`, Hz amagat.
    width_broadening : numpy.ndarray
        The broadening coefficient :math:`\beta`, Hz/amagat.
    doppler_width : numpy.ndarray
        The forward-scattering Doppler width that caps a Dicke term, Hz (0 where the
        line has no Dicke term).
    shift : numpy.ndarray
        The line's frequency shift per unit density, Hz/amagat.
    damping_group : numpy.ndarray
        Lines with the same number share their width and shift at every density.
    density_range : tuple of float
        The number densities between which the lines are independent collisional
        Lorentzians (module docstring), m^-3.
    temperature_range : tuple of float
        The temperatures covered by the width measurements of every band, K.
    sources : dict of str to str
        The datasets and formulas behind the response, by what they give.
    reasons : dict of str to str
        Why each choice was made.
    material : str
        What the response is of: a species, a named mixture (``"dry air"``) or, for
        mole fractions given directly, the species joined by ``+``.
    population_floor : float
        The relative population below which rotational levels were left out.
    status : str
        ``validated``, or ``provisional`` where a constituent's widths are a stopgap
        (CO2 and N2O, TC-542); ``user`` for a response a caller built.
    """

    IDENTITY_KIND: ClassVar[str] = "raman_response"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.raman_response",)
    INPUTS: ClassVar[tuple[str, ...]] = ("amplitude_scale", "width_law", "width_scale")
    DERIVED: ClassVar[tuple[str, ...]] = ("damping_group", "density_range")

    composition: tuple[tuple[str, float], ...]
    temperature: float
    wavelength: float | None
    species: tuple[str, ...]
    band: tuple[str, ...]
    branch: tuple[str, ...]
    j: Int[np.ndarray, " n"]
    frequency: Float[np.ndarray, " n"]
    isotropic: Float[np.ndarray, " n"]
    anisotropic: Float[np.ndarray, " n"]
    width_constant: Float[np.ndarray, " n"]
    width_dicke: Float[np.ndarray, " n"]
    width_broadening: Float[np.ndarray, " n"]
    doppler_width: Float[np.ndarray, " n"]
    shift: Float[np.ndarray, " n"]
    damping_group: Int[np.ndarray, " n"]
    density_range: tuple[float, float]
    temperature_range: tuple[float, float]
    sources: Mapping[str, str] = field(default_factory=dict)
    reasons: Mapping[str, str] = field(default_factory=dict)
    material: str = field(default="", kw_only=True)
    population_floor: float = field(default=math.nan, kw_only=True)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the response's name, material and source records.

        Returns
        -------
        tuple
            ``("raman_response/<material>@<T>K[,<pump>nm]", material, sources)``.
        """
        pump = "" if self.wavelength is None else f",{self.wavelength * 1e9:.6g}nm"
        name = f"raman_response/{self.material}@{self.temperature:g}K{pump}"
        return name, self.material, source_ids(*self.sources.values())

    @classmethod
    def from_lines(  # noqa: PLR0913, PLR0917 - the lines' arrays, as OEP 1 names them
        cls,
        species: str | Sequence[str],
        temperature: float,
        frequency: ArrayLike,
        isotropic: ArrayLike,
        anisotropic: ArrayLike,
        width_constant: ArrayLike,
        width_dicke: ArrayLike,
        width_broadening: ArrayLike,
        doppler_width: ArrayLike,
        shift: ArrayLike,
        *,
        band: Sequence[str] | None = None,
        branch: Sequence[str] | None = None,
        j: ArrayLike | None = None,
        mass: float | Mapping[str, float] | None = None,
        composition: Mapping[str, float] | None = None,
        wavelength: float | None = None,
        temperature_range: tuple[float, float] | None = None,
        citation: str = "",
        material: str = "",
        population_floor: float = math.nan,
        sources: Mapping[str, str] | None = None,
        reasons: Mapping[str, str] | None = None,
    ) -> "GasRamanResponse":
        """Return a Raman response from its lines, as :func:`raman_response` does.

        The damping groups number the lines by first appearance of their species,
        band and width law ``(c, d, beta, cap, shift)``. The density range is each
        band's (TC-688) for its strongest line, with that line's forward-scattering
        Doppler width (Herring et al.'s Eq. 4, from the species' mass): the
        narrowest over the bands.

        Parameters
        ----------
        species : str or sequence of str
            The species of every line, or of each.
        temperature : float
            K.
        frequency : array_like
            Each line's frequency at zero density, rad/s.
        isotropic, anisotropic : array_like
            The two shares of each line's amplitude, m^5 V^-2 s^-1.
        width_constant, width_dicke, width_broadening : array_like
            ``c`` (Hz), ``d`` (Hz amagat) and ``beta`` (Hz/amagat).
        doppler_width : array_like
            The Doppler cap of each line's Dicke term, Hz; 0 without one.
        shift : array_like
            The frequency shift per unit density, Hz/amagat.
        band, branch : sequence of str or None, optional
            Each line's band (``"rotation"`` by default) and branch (``"S"``).
        j : array_like or None, optional
            Each line's lower level; -1 (unresolved) by default.
        mass : float or Mapping or None, optional
            The molecular mass in u, one or by species; the curated one of a curated
            species by default.
        composition : Mapping or None, optional
            Species to mole fraction; each species at its share of the lines' by
            default (one species: 1).
        wavelength : float or None, optional
            The pump wavelength the strengths are dispersed at, m.
        temperature_range : tuple of float or None, optional
            The temperatures the widths hold over, K; the temperature by default.
        citation : str, optional
            Where the lines come from.
        material : str, optional
            What the response is of; the species joined by ``+`` by default.
        population_floor : float, optional
            The population below which levels were left out, if one was.
        sources, reasons : Mapping or None, optional
            The response's sources and reasons.

        Returns
        -------
        GasRamanResponse
            The response, status ``user``.

        Raises
        ------
        ValueError
            If the arrays differ in length, a Dicke term has no cap, or a species'
            mass is unknown.

        Examples
        --------
        >>> response = GasRamanResponse.from_lines(
        ...     "X2", 295.0, [1e13], [0.0], [1e-40], [1e9], [0.0], [1e8], [0.0], [0.0],
        ...     mass=30.0)
        >>> response.identity.material, int(response.damping_group[0])
        ('X2', 0)
        """
        freq = np.asarray(frequency, dtype=float)
        n = freq.size
        per_line = [species] * n if isinstance(species, str) else list(species)
        arrays = {
            name: np.asarray(value, dtype=float)
            for name, value in (
                ("isotropic", isotropic),
                ("anisotropic", anisotropic),
                ("width_constant", width_constant),
                ("width_dicke", width_dicke),
                ("width_broadening", width_broadening),
                ("doppler_width", doppler_width),
                ("shift", shift),
            )
        }
        bands = list(band) if band is not None else ["rotation"] * n
        branches = list(branch) if branch is not None else ["S"] * n
        levels = np.full(n, -1) if j is None else np.asarray(j, dtype=int)
        if any(
            len(x) != n for x in (per_line, bands, branches, levels, *arrays.values())
        ):
            msg = "every array of a Raman response needs one entry per line"
            raise ValueError(msg)
        uncapped = (arrays["width_dicke"] > 0) & ~(arrays["doppler_width"] > 0)
        if np.any(uncapped):
            msg = "a Dicke term needs its Doppler cap (doppler_width > 0)"
            raise ValueError(msg)
        masses = _masses(per_line, mass)
        present = list(dict.fromkeys(per_line))
        fractions = (
            dict(composition)
            if composition is not None
            else {s: per_line.count(s) / n for s in present}
        )
        return cls(
            composition=tuple((s, float(x)) for s, x in fractions.items()),
            temperature=temperature,
            wavelength=wavelength,
            species=tuple(per_line),
            band=tuple(bands),
            branch=tuple(branches),
            j=levels,
            frequency=freq,
            damping_group=_damping_groups(per_line, bands, arrays),
            density_range=_density_range(
                per_line, bands, freq, arrays, temperature, masses
            ),
            temperature_range=temperature_range or (temperature, temperature),
            sources=({"citation": citation} if citation else {}) | dict(sources or {}),
            reasons=dict(reasons or {}),
            material=material or "+".join(present),
            population_floor=population_floor,
            isotropic=arrays["isotropic"],
            anisotropic=arrays["anisotropic"],
            width_constant=arrays["width_constant"],
            width_dicke=arrays["width_dicke"],
            width_broadening=arrays["width_broadening"],
            doppler_width=arrays["doppler_width"],
            shift=arrays["shift"],
        )

    def _lines(self) -> dict[str, Any]:
        """Return the arguments of :meth:`from_lines` that rebuild this response.

        Returns
        -------
        dict
            The lines, the state and the descriptions.
        """
        return {
            "species": self.species,
            "temperature": self.temperature,
            "frequency": self.frequency,
            "isotropic": self.isotropic,
            "anisotropic": self.anisotropic,
            "width_constant": self.width_constant,
            "width_dicke": self.width_dicke,
            "width_broadening": self.width_broadening,
            "doppler_width": self.doppler_width,
            "shift": self.shift,
            "band": self.band,
            "branch": self.branch,
            "j": self.j,
            "composition": dict(self.composition),
            "wavelength": self.wavelength,
            "temperature_range": self.temperature_range,
            "material": self.material,
            "population_floor": self.population_floor,
            "sources": self.sources,
            "reasons": self.reasons,
        }

    def _override(self, changes: Mapping[str, Any]) -> "GasRamanResponse":
        """Return the response with its lines' amplitudes or widths changed.

        Parameters
        ----------
        changes : Mapping
            ``amplitude_scale``, ``width_scale`` and ``width_law``.

        Returns
        -------
        GasRamanResponse
            The response rebuilt by :meth:`from_lines`, its groups and range
            recomputed.

        Raises
        ------
        ValueError
            If ``width_law`` names a band the response does not have, or a law is
            not four numbers.
        """
        lines = self._lines()
        scale = float(changes.get("amplitude_scale", 1.0))
        lines["isotropic"] = self.isotropic * scale
        lines["anisotropic"] = self.anisotropic * scale
        width = float(changes.get("width_scale", 1.0))
        names = ("width_constant", "width_dicke", "width_broadening", "doppler_width")
        for name in names:
            lines[name] = np.array(getattr(self, name)) * width
        labels = np.array(
            [f"{s} {b}" for s, b in zip(self.species, self.band, strict=True)]
        )
        for key, law in dict(changes.get("width_law", {})).items():
            chosen = (np.array(self.band) == key) | (labels == key)
            if not np.any(chosen) or len(law) != len(names):
                msg = f"width_law: no band {key!r}, or its law is not (c, d, beta, cap)"
                raise ValueError(msg)
            for name, value in zip(names, law, strict=True):
                lines[name][chosen] = float(value)
        return type(self).from_lines(**lines)

    def __len__(self) -> int:
        """Return the number of lines.

        Returns
        -------
        int
            The number of terms of the response.
        """
        return len(self.branch)

    @property
    def amplitude(self) -> Float[np.ndarray, " n"]:
        """The amplitude of each line for a single linear polarisation, m^5 V^-2 s^-1.

        Returns
        -------
        numpy.ndarray
            ``isotropic + anisotropic``.
        """
        return self.isotropic + self.anisotropic

    @property
    def references(self) -> tuple[str, ...]:
        """The curated sources the response rests on.

        Returns
        -------
        tuple of str
            Source ids, sorted.
        """
        ids = {
            reference.strip().split("/")[0]
            for value in self.sources.values()
            for reference in value.split(",")
            if "/" in reference
        }
        return tuple(sorted(ids))

    def component_amplitudes(
        self, component: Component = "parallel"
    ) -> Float[np.ndarray, " n"]:
        """Return each line's amplitude in one part of the response.

        Parameters
        ----------
        component : {"parallel", "isotropic", "anisotropic", "a", "b"}, optional
            The part (:data:`Component`).

        Returns
        -------
        numpy.ndarray
            The amplitudes, m^5 V^-2 s^-1.

        Raises
        ------
        DataError
            If the component is not one of these.
        """
        if component not in get_args(Component):
            msg = f"unknown component {component!r}; one of {get_args(Component)}"
            raise DataError(msg)
        iso, aniso = _COMPONENT_WEIGHTS[component]
        return iso * self.isotropic + aniso * self.anisotropic

    def _densities(self, density: ArrayLike) -> Float[np.ndarray, " *shape"]:
        """Return densities in amagat, warning outside the model's range.

        Parameters
        ----------
        density : array_like
            Number densities, m^-3, > 0.

        Returns
        -------
        numpy.ndarray
            The densities in amagat.

        Raises
        ------
        DataError
            If a density is not positive.

        Warns
        -----
        ModelRangeWarning
            Outside :attr:`density_range`.
        """
        n = np.asarray(density, dtype=float)
        if not np.all(n > 0):
            msg = "the Raman response needs positive densities"
            raise DataError(msg)
        self._warn_outside(n)
        return np.asarray(n / AMAGAT)

    def _warn_outside(self, n: np.ndarray) -> None:
        """Warn for densities outside the validated range.

        Parameters
        ----------
        n : numpy.ndarray
            Number densities, m^-3.

        Warns
        -----
        ModelRangeWarning
            Outside :attr:`density_range`.
        """
        low, high = self.density_range
        if np.any(n < low):
            msg = (
                f"density below {low / AMAGAT:.3g} amagat, where the motion of the "
                "molecules sets the widths: the Dicke terms are capped at the Doppler "
                "widths (TC-688)"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)
        if np.any(n > high):
            msg = (
                f"density above {high / AMAGAT:.3g} amagat, where the lines overlap: "
                "their mixing is not modelled (TC-527, TC-688)"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=4)

    def check_densities(self, density: ArrayLike) -> None:
        """Warn once for the densities a consumer will evaluate the groups at.

        The group view (:meth:`groups`) does not warn; a propagation code declares
        its densities here instead.

        Parameters
        ----------
        density : array_like
            Number densities, m^-3, >= 0.

        Raises
        ------
        DataError
            If a density is negative or not finite.

        Warns
        -----
        ModelRangeWarning
            If any is outside :attr:`density_range`.
        """
        n = np.asarray(density, dtype=float)
        if not np.all(np.isfinite(n) & (n >= 0)):
            msg = "densities must be finite and >= 0"
            raise DataError(msg)
        self._warn_outside(n)

    @property
    def line_id(self) -> tuple[str, ...]:
        """A stable identifier of each line: ``species/band/branch/J``.

        ``J`` is the lower level, or ``-`` for an unresolved band
        (``"N2/rotation/S/6"``, ``"H2/v1/Q/1"``).

        Returns
        -------
        tuple of str
            One per line.
        """
        return tuple(
            f"{s}/{b}/{r}/{'-' if j < 0 else int(j)}"
            for s, b, r, j in zip(
                self.species, self.band, self.branch, self.j, strict=True
            )
        )

    def groups(self) -> tuple[RamanGroup, ...]:
        """Return the damping groups, with their laws in SI (OEP 1, section 6.7).

        Returns
        -------
        tuple of RamanGroup
            In the order of :attr:`damping_group`; their responses sum to
            :meth:`response`.

        Examples
        --------
        >>> h2 = raman_response("H2", 295.0)
        >>> [g.line_ids[0] for g in h2.groups() if g.shifted][:1]
        ['H2/v1/Q/0']
        """
        ids = self.line_id
        out = []
        for index in range(int(np.max(self.damping_group)) + 1):
            chosen = np.flatnonzero(self.damping_group == index)
            k = int(chosen[0])
            out.append(
                RamanGroup(
                    index=index,
                    lines=tuple(int(i) for i in chosen),
                    line_ids=tuple(ids[i] for i in chosen),
                    gamma_constant=math.pi * float(self.width_constant[k]),
                    gamma_dicke=math.pi * float(self.width_dicke[k]) * AMAGAT,
                    gamma_doppler=math.pi * float(self.doppler_width[k]),
                    gamma_broadening=math.pi * float(self.width_broadening[k]) / AMAGAT,
                    shift=2 * math.pi * float(self.shift[k]) / AMAGAT,
                    frequency=np.array(self.frequency[chosen]),
                    isotropic=np.array(self.isotropic[chosen]),
                    anisotropic=np.array(self.anisotropic[chosen]),
                )
            )
        return tuple(out)

    def _fwhm(
        self, amagat: Float[np.ndarray, " *shape"]
    ) -> Float[np.ndarray, " n *shape"]:
        """Return the widths at densities in amagat, without range checks.

        Parameters
        ----------
        amagat : numpy.ndarray
            Densities, amagat.

        Returns
        -------
        numpy.ndarray
            FWHM, Hz, shape ``(len(self),) + amagat.shape``.
        """
        n = amagat[None, ...]
        shape = (-1,) + (1,) * amagat.ndim
        dicke = self.width_dicke.reshape(shape) / n
        cap = self.doppler_width.reshape(shape)
        dicke = np.where(cap > 0, np.minimum(dicke, cap), dicke)
        return (
            self.width_constant.reshape(shape)
            + dicke
            + self.width_broadening.reshape(shape) * n
        )

    def fwhm(self, density: ArrayLike) -> Float[np.ndarray, " n *shape"]:
        """Return every line's full width at half maximum at number densities.

        Parameters
        ----------
        density : array_like
            Total number densities, m^-3, > 0.

        Returns
        -------
        numpy.ndarray
            FWHM, Hz, shape ``(len(self),) + shape(density)``.

        Warns
        -----
        ModelRangeWarning
            Outside :attr:`density_range`.
        """
        return self._fwhm(self._densities(density))

    def decay_rates(self, density: ArrayLike) -> Float[np.ndarray, " n *shape"]:
        r"""Return every line's dephasing rate :math:`\Gamma = \pi \Delta\nu`.

        Parameters
        ----------
        density : array_like
            Total number densities, m^-3, > 0.

        Returns
        -------
        numpy.ndarray
            :math:`\Gamma`, 1/s, shape ``(len(self),) + shape(density)``.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.

        Examples
        --------
        >>> from optmat.gas.linewidths import AMAGAT
        >>> response = raman_response("H2", 295.0, wavelength=800e-9)
        >>> rates = response.decay_rates(np.array([2.0, 10.0]) * AMAGAT)
        >>> rates.shape == (len(response), 2)
        True
        """
        return math.pi * self._fwhm(self._densities(density))

    def frequencies(self, density: ArrayLike) -> Float[np.ndarray, " n *shape"]:
        r"""Return every line's angular frequency :math:`\Omega` at number densities.

        Parameters
        ----------
        density : array_like
            Total number densities, m^-3, > 0.

        Returns
        -------
        numpy.ndarray
            :math:`\Omega`, rad/s, shape ``(len(self),) + shape(density)``.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.
        """
        amagat = self._densities(density)
        shape = (-1,) + (1,) * amagat.ndim
        return (
            self.frequency.reshape(shape)
            + 2 * math.pi * self.shift.reshape(shape) * amagat[None, ...]
        )

    def _state(self, density: float) -> tuple[np.ndarray, np.ndarray]:
        """Return the lines' frequencies and decay rates at one density.

        Parameters
        ----------
        density : float
            Number density, m^-3.

        Returns
        -------
        tuple of numpy.ndarray
            Omega (rad/s) and Gamma (1/s).
        """
        amagat = self._densities(np.asarray(float(density)))
        omega = self.frequency + 2 * math.pi * self.shift * amagat
        return omega, math.pi * self._fwhm(amagat)

    def response(
        self,
        time: ArrayLike,
        density: float,
        component: Component = "parallel",
    ) -> Float[np.ndarray, " *shape"]:
        r"""Return the response function :math:`h(\tau)` at a number density.

        Parameters
        ----------
        time : array_like
            Delay after the impulse, s (the response is zero before it).
        density : float
            Total number density, m^-3, > 0.
        component : {"parallel", "isotropic", "anisotropic", "a", "b"}, optional
            The part of the response (:data:`Component`).

        Returns
        -------
        numpy.ndarray
            :math:`h`, m^5 V^-2 s^-1, with the shape of ``time``.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.

        Examples
        --------
        Nitrogen's response first peaks 75 fs after the impulse, as the molecules
        begin to align:

        >>> from optmat.gas.linewidths import AMAGAT
        >>> response = raman_response("N2", 295.0)
        >>> t = np.linspace(0.0, 0.2e-12, 2001)
        >>> f"{t[np.argmax(response.response(t, AMAGAT))] * 1e15:.0f}"
        '75'
        """
        tau = np.asarray(time, dtype=float)
        omega, rate = self._state(density)
        amplitude = self.component_amplitudes(component)
        positive = np.maximum(tau, 0.0)
        total = np.zeros(tau.shape)
        for start in range(0, len(self), _TERM_BLOCK):
            block = slice(start, start + _TERM_BLOCK)
            t = positive[..., None]
            total += np.sum(
                amplitude[block] * np.sin(omega[block] * t) * np.exp(-rate[block] * t),
                axis=-1,
            )
        return np.where(tau >= 0, total, 0.0)

    def spectrum(
        self,
        omega: ArrayLike,
        density: float,
        component: Component = "parallel",
    ) -> Complex[np.ndarray, " *shape"]:
        r"""Return :math:`H(\omega) = \int_0^\infty h(\tau) e^{i\omega\tau} d\tau`.

        Parameters
        ----------
        omega : array_like
            Angular frequency, rad/s.
        density : float
            Total number density, m^-3, > 0.
        component : {"parallel", "isotropic", "anisotropic", "a", "b"}, optional
            The part of the response (:data:`Component`).

        Returns
        -------
        numpy.ndarray
            :math:`H`, m^5 V^-2, complex, with the shape of ``omega``. Its imaginary
            part is positive at a Stokes resonance (gain).

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.
        """
        w = np.asarray(omega, dtype=float)
        frequency, rate = self._state(density)
        amplitude = self.component_amplitudes(component)
        total = np.zeros(w.shape, dtype=complex)
        for start in range(0, len(self), _TERM_BLOCK):
            block = slice(start, start + _TERM_BLOCK)
            x = w[..., None]
            total += np.sum(
                amplitude[block]
                * frequency[block]
                / (
                    frequency[block] ** 2
                    + rate[block] ** 2
                    - x**2
                    - 2j * rate[block] * x
                ),
                axis=-1,
            )
        return total

    def lines(self, density: float) -> RamanLines:
        """Return the lines at a number density, for :mod:`optmat.gas.raman`'s formulas.

        Parameters
        ----------
        density : float
            Total number density, m^-3, > 0.

        Returns
        -------
        RamanLines
            Frequencies, amplitudes (parallel) and widths at that density.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.
        """
        omega, rate = self._state(density)
        return RamanLines(
            branch=self.branch,
            j=self.j,
            frequency=omega,
            amplitude=self.amplitude,
            fwhm=rate / math.pi,
        )

    def stokes_gain(
        self, pump_wavelength: float, density: float
    ) -> Float[np.ndarray, " n"]:
        """Return each line's steady-state Stokes gain, pump and Stokes polarised alike.

        :func:`optmat.gas.raman.stokes_gain` of :meth:`lines`.

        Parameters
        ----------
        pump_wavelength : float
            Vacuum wavelength of the pump, m.
        density : float
            Total number density, m^-3, > 0.

        Returns
        -------
        numpy.ndarray
            The gain of each line, m/W.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.
        """
        return stokes_gain(self.lines(density), pump_wavelength, density)

    def rotational_nonlinear_index(self, density: float) -> float:
        """Return the rotational nonlinear index for a long pulse.

        :func:`optmat.gas.raman.rotational_nonlinear_index` of :meth:`lines`.

        Parameters
        ----------
        density : float
            Total number density, m^-3, > 0.

        Returns
        -------
        float
            ``n2``, m^2/W.

        Warns
        -----
        ModelRangeWarning
            As :meth:`fwhm`.
        """
        return rotational_nonlinear_index(self.lines(density), density)


# ---------------------------------------------------------------------------------
# Building a response from the curated lines.


@dataclass(frozen=True)
class _GasChoices:
    """The response's choices for one gas, from ``raman_response.toml``.

    Attributes
    ----------
    temperature_ranges : dict of str to tuple of float
        By band kind (``rotation``, ``vibration``): the temperatures of the width
        measurements, K.
    mass : float
        The molecular mass in the Doppler width, u.
    shifts : dict of str to dict of str to str
        By band kind: the shift law's entry.
    sources : dict of str to str
        By what they give.
    reasons : dict of str to str
        Why.
    status : str
        The status the gas's lines give a response (``validated`` unless the file
        says otherwise).
    status_reason : str
        Why.
    """

    temperature_ranges: Mapping[str, tuple[float, float]]
    mass: float
    shifts: Mapping[str, Mapping[str, str]]
    sources: Mapping[str, str]
    reasons: Mapping[str, str]
    status: Status = "validated"
    status_reason: str = ""


def _pair(value: object, name: str) -> tuple[float, float]:
    """Return a TOML array of two increasing numbers.

    Parameters
    ----------
    value : object
        The value read.
    name : str
        Its key, for the message.

    Returns
    -------
    tuple of float
        The two numbers.

    Raises
    ------
    DataError
        If it is not such a pair.
    """
    if not isinstance(value, list) or len(value) != 2:  # pyright: ignore[reportUnknownArgumentType]
        msg = f"raman_response.toml: {name} must be [low, high]"
        raise DataError(msg)
    low, high = (float(str(v)) for v in value)  # pyright: ignore[reportUnknownVariableType]
    if not 0 < low < high:
        msg = f"raman_response.toml: {name} must increase from a positive value"
        raise DataError(msg)
    return low, high


def _strings(value: object) -> dict[str, str]:
    """Return a TOML table as strings.

    Parameters
    ----------
    value : object
        The table (or nothing).

    Returns
    -------
    dict of str to str
        Its keys and values.
    """
    if not isinstance(value, Mapping):
        return {}
    return {str(k): str(v) for k, v in value.items()}  # pyright: ignore[reportUnknownVariableType]


@functools.cache
def _config() -> tuple[dict[Gas, _GasChoices], dict[str, dict[str, object]]]:
    """Read ``raman_response.toml``.

    Returns
    -------
    tuple of (dict of Gas to _GasChoices, dict of str to dict)
        The gases' choices, and the mixtures' entries.
    """
    with RAMAN_RESPONSE_FILE.open("rb") as handle:
        table = tomllib.load(handle)
    gases: dict[Gas, _GasChoices] = {}
    for name, entry in dict(table.get("gases", {})).items():
        shifts = {
            str(kind): _strings(sub)
            for kind, sub in dict(entry.get("shifts", {})).items()
        }
        ranges = {
            str(kind): _pair(value, f"{name} {kind} temperatures")
            for kind, value in dict(entry["temperature_range_k"]).items()
        }
        status, status_reason = parse_status(
            entry, f"raman_response.toml [gases.{name}]", default="validated"
        )
        gases[as_gas(name)] = _GasChoices(
            temperature_ranges=ranges,
            mass=float(str(entry["mass_u"])),
            shifts=shifts,
            sources=_strings(entry.get("sources")),
            reasons=_strings(entry.get("reasons")),
            status=status,
            status_reason=status_reason,
        )
    mixtures = {str(k): dict(v) for k, v in dict(table.get("mixtures", {})).items()}
    return gases, mixtures


def _formula_coefficients(reference: str) -> dict[str, float]:
    """Return a curated formula's numeric coefficients.

    Parameters
    ----------
    reference : str
        ``source/formula``.

    Returns
    -------
    dict of str to float
        The coefficients that are numbers.

    Raises
    ------
    DataError
        If the formula does not exist.
    """
    source_id, _, formula_id = reference.partition("/")
    matches = [f for f in sources_of()[source_id].formulas if f["id"] == formula_id]
    if len(matches) != 1:
        msg = f"raman_response.toml names an unknown formula {reference!r}"
        raise DataError(msg)
    numbers: dict[str, float] = {}
    for key, value in _strings(matches[0]["coefficients"]).items():
        try:
            numbers[key] = float(value)
        except ValueError:
            continue
    return numbers


def _shift(entry: Mapping[str, str], temperature: float) -> float:
    """Return a band's line shift per amagat at a temperature.

    Parameters
    ----------
    entry : mapping of str to str
        ``law = "polynomial"``: ``formula``, whose coefficients ``C0``, ``C1``, ``C2``
        and ``T0`` give ``C = C0 + C1 (T - T0) + C2 (T - T0)^2`` in MHz/amagat
        (Bischel & Dyer 1986, Eq. (4)).
    temperature : float
        Temperature, K.

    Returns
    -------
    float
        The shift of the Raman frequency, Hz/amagat.

    Raises
    ------
    DataError
        If the law is unknown.
    """
    if entry.get("law") != "polynomial":
        msg = f"raman_response.toml: unknown shift law {entry.get('law')!r}"
        raise DataError(msg)
    c = _formula_coefficients(entry["formula"])
    dt = temperature - c["T0"]
    return (c["C0"] + c["C1"] * dt + c["C2"] * dt**2) * 1e6


def _doppler_constant() -> float:
    """Return the constant of the forward-scattering Doppler width.

    Herring, Dyer & Bischel 1986, Eq. (4): ``Delta nu_D = k nu_R (T/m)^(1/2)``, FWHM,
    with ``nu_R`` and ``Delta nu_D`` in one unit, ``T`` in K and ``m`` in u
    (``herring_1986/eq4_doppler``).

    Returns
    -------
    float
        ``k``, K^-1/2 u^1/2.
    """
    return _formula_coefficients("herring_1986/eq4_doppler")["k"]


def _cutoff_factor() -> float:
    """Return the factor of Herring et al.'s cutoff density of the diffusion law.

    Herring, Dyer & Bischel 1986, Eq. (3): ``rho_c = 3.33 A / Delta nu_D``, the density
    at which the diffusion (Dicke) law is 10 % above the hard-collision line shape
    (``herring_1986/eq3_cutoff``).

    Returns
    -------
    float
        The factor.
    """
    return _formula_coefficients("herring_1986/eq3_cutoff")["factor"]


def _width_coefficients(
    widths: tuple[Float[np.ndarray, " n"], ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``c``, ``d`` and ``beta`` of ``c + d/n + beta n`` from widths at probes.

    Parameters
    ----------
    widths : tuple of numpy.ndarray
        Each line's FWHM (Hz) at the densities of :data:`_PROBE_DENSITIES`, then at
        :data:`_CHECK_DENSITY`.

    Returns
    -------
    tuple of numpy.ndarray
        ``c`` (Hz), ``d`` (Hz amagat), ``beta`` (Hz/amagat).

    Raises
    ------
    DataError
        If the coefficients do not reproduce the width at the fourth density.
    """
    design = np.array([[1.0, 1.0 / n, n] for n in _PROBE_DENSITIES])
    constant, dicke, broadening = np.linalg.solve(design, np.array(widths[:3]))
    predicted = constant + dicke / _CHECK_DENSITY + broadening * _CHECK_DENSITY
    scale = np.maximum(np.abs(widths[3]), np.finfo(float).tiny)
    if np.any(np.abs(predicted - widths[3]) > _LAW_TOLERANCE * scale):
        msg = "a width law is not of the form c + d/n + beta n"
        raise DataError(msg)
    # A coefficient the law does not have comes out at rounding level; zero it, so
    # that lines of one law stay in one damping group and a Dicke cap is not applied.
    size = np.max(np.abs(np.array(widths[:3])), axis=0)
    for coefficient, weight in ((constant, 1.0), (dicke, 1.0), (broadening, 1.0)):
        coefficient[np.abs(coefficient) * weight < 1e-12 * size] = 0.0
    return constant, dicke, broadening


def _anisotropic_share(part: BandLines) -> Float[np.ndarray, " n"]:
    """Return the anisotropic (rank-2) share of each line of a band (TC-686).

    Parameters
    ----------
    part : BandLines
        One band's lines; their amplitudes are the band's isotropic shares for a
        vibrational band, the anisotropic ones for the rotational band.

    Returns
    -------
    numpy.ndarray
        The anisotropic shares, m^5 V^-2 s^-1.
    """
    lines, band = part.lines, part.band
    if isinstance(band, RotationalBand):
        return lines.amplitude
    isotropic = lines.amplitude
    if band.transition_anisotropy is not None and band.lower is not None:
        # Per level: (4/45) b_JJ gamma'^2 against a'^2, both dispersed alike.
        ratio = band.transition_anisotropy(lines.j) / _mean_element(band, lines.j)
        return isotropic * 4 / 45 * placzek_teller_q(lines.j) * ratio**2
    rho = band.depolarization or 0.0
    band_ratio = 4 * rho / (3 - 4 * rho)  # (4/45) b gamma'^2 / a'^2 of the branch
    if band.lower is None:
        return isotropic * band_ratio
    # Resolved with one ratio for the branch: shared out as b_JJ over its mean.
    b = placzek_teller_q(lines.j)
    mean_b = float(np.sum(isotropic * b) / np.sum(isotropic))
    return isotropic * band_ratio * b / mean_b


def _mean_element(band: VibrationalBand, j: Int[np.ndarray, " n"]) -> np.ndarray:
    """Return a band's curated mean transition element at levels ``j``, C m^2/V.

    Parameters
    ----------
    band : VibrationalBand
        The band.
    j : numpy.ndarray
        Rotational levels.

    Returns
    -------
    numpy.ndarray
        ``<0J| alpha_mean |1J>`` as curated (not dispersed).
    """
    element = band.transition_polarizability
    if isinstance(element, AnisotropyTable):
        return element(j)
    return np.full(np.shape(j), element)


def _band_limits(
    amplitude: np.ndarray,
    frequency: np.ndarray,
    constant: np.ndarray,
    dicke: np.ndarray,
    broadening: np.ndarray,
    doppler: np.ndarray,
) -> tuple[float, float]:
    """Return the densities between which a band's lines are collisional Lorentzians.

    Found for the band's strongest line (module docstring, TC-688).

    Parameters
    ----------
    amplitude : numpy.ndarray
        The band's lines' amplitudes (any common factor).
    frequency : numpy.ndarray
        Their frequencies, rad/s.
    constant, dicke, broadening : numpy.ndarray
        Their width coefficients: Hz, Hz amagat, Hz/amagat.
    doppler : numpy.ndarray
        Their Doppler widths for forward scattering, Hz.

    Returns
    -------
    tuple of float
        The lower and upper density, amagat (0 and inf where there is no limit).
    """
    k = int(np.argmax(np.abs(amplitude)))
    if dicke[k] > 0:
        low = _cutoff_factor() * dicke[k] / doppler[k]
    elif constant[k] >= doppler[k] or broadening[k] <= 0:
        low = 0.0
    else:
        low = (doppler[k] - constant[k]) / broadening[k]
    high = math.inf
    if len(frequency) > 1 and broadening[k] > 0:
        nu = frequency / (2 * math.pi)
        spacing = float(np.min(np.abs(np.delete(nu, k) - nu[k])))
        high = max(spacing - constant[k], 0.0) / broadening[k]
    return float(low), float(high)


def _masses(
    species: Sequence[str], mass: float | Mapping[str, float] | None
) -> dict[str, float]:
    """Return the molecular mass of every species of a response's lines.

    Parameters
    ----------
    species : sequence of str
        The lines' species.
    mass : float or Mapping or None
        A caller's mass for all, or by species; else the curated ones.

    Returns
    -------
    dict
        Species to mass, u.

    Raises
    ------
    ValueError
        If a species has neither a given nor a curated mass.
    """
    curated = {g.value: c.mass for g, c in _config()[0].items()}
    given = (
        dict.fromkeys(species, float(mass))
        if isinstance(mass, (int, float))
        else dict(mass or {})
    )
    out = {}
    for name in dict.fromkeys(species):
        if name not in given and name not in curated:
            msg = f"no molecular mass for {name}: give mass="
            raise ValueError(msg)
        out[name] = float(given.get(name, curated.get(name, math.nan)))
    return out


def _damping_groups(
    species: Sequence[str], bands: Sequence[str], arrays: Mapping[str, np.ndarray]
) -> np.ndarray:
    """Return the damping groups: lines of one band and one width law share one.

    Parameters
    ----------
    species, bands : sequence of str
        Each line's species and band.
    arrays : Mapping
        The lines' ``width_constant``, ``width_dicke``, ``width_broadening``,
        ``doppler_width`` and ``shift``.

    Returns
    -------
    numpy.ndarray
        The group of each line, numbered by first appearance.
    """
    groups: dict[tuple[object, ...], int] = {}
    keys = zip(
        species,
        bands,
        arrays["width_constant"],
        arrays["width_dicke"],
        arrays["width_broadening"],
        arrays["doppler_width"],
        arrays["shift"],
        strict=True,
    )
    return np.array([groups.setdefault(tuple(key), len(groups)) for key in keys])


def _density_range(
    species: Sequence[str],
    bands: Sequence[str],
    frequency: np.ndarray,
    arrays: Mapping[str, np.ndarray],
    temperature: float,
    masses: Mapping[str, float],
) -> tuple[float, float]:
    """Return the densities where every band's lines are collisional Lorentzians.

    Parameters
    ----------
    species, bands : sequence of str
        Each line's species and band.
    frequency : numpy.ndarray
        The lines' frequencies, rad/s.
    arrays : Mapping
        The lines' amplitudes and width coefficients (:meth:`GasRamanResponse
        .from_lines`).
    temperature : float
        K.
    masses : Mapping
        Each species' molecular mass, u.

    Returns
    -------
    tuple of float
        The lower and upper density, m^-3.
    """
    labels = np.array([f"{s}\0{b}" for s, b in zip(species, bands, strict=True)])
    limits = []
    for label in dict.fromkeys(labels):
        chosen = labels == label
        name = str(label).split("\0", 1)[0]
        isotropic = arrays["isotropic"][chosen]
        # A band's lines carry their curated amplitude as the isotropic share
        # (vibration) or the anisotropic one (rotation).
        amplitude = isotropic if np.any(isotropic) else arrays["anisotropic"][chosen]
        doppler = (
            _doppler_constant()
            * frequency[chosen]
            / (2 * math.pi)
            * math.sqrt(temperature / masses[name])
        )
        limits.append(
            _band_limits(
                amplitude,
                frequency[chosen],
                arrays["width_constant"][chosen],
                arrays["width_dicke"][chosen],
                arrays["width_broadening"][chosen],
                doppler,
            )
        )
    return (
        max(low for low, _ in limits) * AMAGAT,
        min(high for _, high in limits) * AMAGAT,
    )


def _composition(
    gas: Gas | str | Mapping[str, float],
    mixtures: Mapping[str, Mapping[str, object]],
) -> tuple[tuple[tuple[str, float], ...], dict[str, str], dict[str, str]]:
    """Return the species and mole fractions of a gas or mixture.

    Parameters
    ----------
    gas : Gas or str or mapping of str to float
        A species, a mixture of ``raman_response.toml`` (``"dry air"``), or mole
        fractions by species.
    mixtures : mapping
        The mixtures' entries.

    Returns
    -------
    tuple
        The composition, and the sources and reasons it adds.

    Raises
    ------
    DataError
        If the fractions are not positive or do not sum to one.
    """
    if isinstance(gas, Mapping):
        pairs = tuple((as_gas(k).value, float(v)) for k, v in gas.items())
        total = sum(v for _, v in pairs)
        if any(v <= 0 for _, v in pairs) or not math.isclose(total, 1.0, rel_tol=1e-9):
            msg = "mole fractions must be positive and sum to one"
            raise DataError(msg)
        return pairs, {}, {}
    if gas == "dry air (CIPM-2007)":
        msg = (
            "the full CIPM-2007 air has no Raman response: its CO2 has no "
            "broadening coefficients by N2 and O2; use 'dry air (N2, O2, Ar)' or "
            "'dry air (N2, O2)'"
        )
        raise DataError(msg)
    name = _AIR_ALIAS.get(gas, gas) if isinstance(gas, str) else gas
    if isinstance(name, str) and name in mixtures:
        entry = mixtures[name]
        keep = [str(s) for s in entry["species"]]  # pyright: ignore[reportGeneralTypeIssues]
        constituents = {c.material: c for c in composition(name)}
        if set(constituents) != set(keep):
            msg = f"raman_response.toml's {name!r} differs from its composition"
            raise DataError(msg)
        pairs = tuple((s, constituents[s].mole_fraction) for s in keep)
        sources = {"composition": ", ".join(constituents[s].source for s in keep)}
        return pairs, sources, {"composition": str(entry["reason"])}
    return ((as_gas(gas).value, 1.0),), {}, {}


@dataclass(frozen=True)
class _BandTerms:
    """One band's lines of one species, ready to be joined into a response.

    Attributes
    ----------
    species, name : str
        The species and the band (``rotation`` or the mode).
    index : int
        The band's place among the species' bands.
    lines : RamanLines
        The lines (their amplitudes as curated, per molecule of the species).
    isotropic, anisotropic : numpy.ndarray
        The two shares of each line, per molecule of the mixture.
    constant, dicke, broadening, doppler : numpy.ndarray
        The width coefficients and the Dicke cap.
    shift : float
        The band's line shift, Hz/amagat.
    temperatures : tuple of float
        The temperatures of its width measurements, K.
    sources : dict of str to str
        Its sources.
    """

    species: str
    name: str
    index: int
    lines: RamanLines
    isotropic: np.ndarray
    anisotropic: np.ndarray
    constant: np.ndarray
    dicke: np.ndarray
    broadening: np.ndarray
    doppler: np.ndarray
    shift: float
    temperatures: tuple[float, float]
    sources: Mapping[str, str]


def _band_terms(
    species: str,
    fraction: float,
    index: int,
    parts: tuple[BandLines, ...],
    choice: _GasChoices,
    temperature: float,
) -> _BandTerms:
    """Return one band's terms from its lines at the probe densities.

    Parameters
    ----------
    species : str
        The species.
    fraction : float
        Its mole fraction in the gas.
    index : int
        The band's place among the species' bands.
    parts : tuple of BandLines
        The band's lines at :data:`_PROBE_DENSITIES` and :data:`_CHECK_DENSITY`.
    choice : _GasChoices
        The species' choices.
    temperature : float
        Temperature, K.

    Returns
    -------
    _BandTerms
        The terms.
    """
    part = parts[0]
    lines = part.lines
    rotational = isinstance(part.band, RotationalBand)
    kind = "rotation" if rotational else "vibration"
    name = "rotation" if isinstance(part.band, RotationalBand) else part.band.key
    constant, dicke, broadening = _width_coefficients(
        tuple(np.asarray(p.lines.fwhm) for p in parts)
    )
    shift = _shift(choice.shifts[kind], temperature) if kind in choice.shifts else 0.0
    doppler = (
        _doppler_constant()
        * lines.frequency
        / (2 * math.pi)
        * math.sqrt(temperature / choice.mass)
    )
    # Only a Dicke term is capped, and with the Doppler width of the band's strongest
    # line: a band's Dicke term is one line's (Q(1) for H2's Q branch, S(1) for its S
    # lines), so its lines keep one width law and one damping group. The other lines'
    # Doppler widths are range criteria only.
    strongest = doppler[int(np.argmax(np.abs(lines.amplitude)))]
    sources = {f"{species} {name}: {k}": v for k, v in part.band.sources.items()}
    sources[f"{species} {name}: widths"] = ", ".join(part.band.widths.sources)
    return _BandTerms(
        species=species,
        name=name,
        index=index,
        lines=lines,
        isotropic=fraction * (np.zeros(len(lines)) if rotational else lines.amplitude),
        anisotropic=fraction * _anisotropic_share(part),
        constant=constant,
        dicke=dicke,
        broadening=broadening,
        doppler=np.where(dicke > 0, strongest, 0.0),
        shift=shift,
        temperatures=choice.temperature_ranges[kind],
        sources=sources,
    )


def _assemble(
    pairs: tuple[tuple[str, float], ...],
    temperature: float,
    wavelength: float | None,
    bands: list[_BandTerms],
    sources: Mapping[str, str],
    reasons: Mapping[str, str],
) -> GasRamanResponse:
    """Join the bands' terms into a response (:meth:`GasRamanResponse.from_lines`).

    Parameters
    ----------
    pairs : tuple of tuple of (str, float)
        The composition.
    temperature : float
        Temperature, K.
    wavelength : float or None
        The pump wavelength, m.
    bands : list of _BandTerms
        The bands.
    sources, reasons : mapping of str to str
        The response's sources and reasons.

    Returns
    -------
    GasRamanResponse
        The response.
    """
    return GasRamanResponse.from_lines(
        tuple(b.species for b in bands for _ in range(len(b.lines))),
        temperature,
        np.concatenate([b.lines.frequency for b in bands]),
        np.concatenate([b.isotropic for b in bands]),
        np.concatenate([b.anisotropic for b in bands]),
        np.concatenate([b.constant for b in bands]),
        np.concatenate([b.dicke for b in bands]),
        np.concatenate([b.broadening for b in bands]),
        np.concatenate([b.doppler for b in bands]),
        np.concatenate([np.full(len(b.lines), b.shift) for b in bands]),
        band=tuple(b.name for b in bands for _ in range(len(b.lines))),
        branch=tuple(x for b in bands for x in b.lines.branch),
        j=np.concatenate([b.lines.j for b in bands]),
        composition=dict(pairs),
        wavelength=wavelength,
        temperature_range=(
            max(b.temperatures[0] for b in bands),
            min(b.temperatures[1] for b in bands),
        ),
        sources=sources,
        reasons=reasons,
    )


def raman_response(
    gas: Gas | str | Mapping[str, float],
    temperature: float,
    *,
    wavelength: float | None = None,
    population_floor: float = 1e-12,
) -> GasRamanResponse:
    """Return the Raman response of a gas or a mixture at a temperature.

    Parameters
    ----------
    gas : Gas or str or mapping of str to float
        A molecular gas (``"N2"``), a mixture named in ``raman_response.toml``
        (``"dry air"``), or mole fractions by species (atoms take part as collision
        partners only).
    temperature : float
        Temperature, K, > 0: the populations and the widths.
    wavelength : float, optional
        Vacuum wavelength of the pump, m, at which the vibrational strengths are
        dispersed (Bischel & Black's law, TC-542); by default as curated.
    population_floor : float, optional
        Lines from levels less populated than this fraction are left out.

    Returns
    -------
    GasRamanResponse
        The response, with its density dependence.

    Raises
    ------
    DataError
        If a species has no Raman parameters or choices, or a collision partner has
        no broadening coefficient.

    Warns
    -----
    ModelRangeWarning
        If the temperature is outside that of the width measurements.

    Examples
    --------
    >>> from optmat.gas.linewidths import AMAGAT
    >>> air = raman_response("dry air", 295.0, wavelength=800e-9)
    >>> sorted(set(air.species))
    ['N2', 'O2']
    >>> f"{air.rotational_nonlinear_index(AMAGAT) * 1e23:.2f}"
    '3.22'
    """
    choices, mixtures = _config()
    pairs, sources, reasons = _composition(gas, mixtures)
    bands: list[_BandTerms] = []
    for species, fraction in pairs:
        gas_id = as_gas(species)
        if info(gas_id).atomic:
            continue
        if gas_id not in choices:
            msg = f"raman_response.toml has no choices for {species}"
            raise DataError(msg)
        choice = choices[gas_id]
        partners = (
            None
            if len(pairs) == 1
            else {("self" if s == species else s): x for s, x in pairs}
        )
        probes = [
            band_lines(
                gas_id,
                temperature,
                n * AMAGAT,
                partners,
                wavelength=wavelength,
                population_floor=population_floor,
            )
            for n in (*_PROBE_DENSITIES, _CHECK_DENSITY)
        ]
        bands += [
            _band_terms(species, fraction, index, parts, choice, temperature)
            for index, parts in enumerate(zip(*probes, strict=True))
        ]
        sources.update({f"{species}: {k}": v for k, v in choice.sources.items()})
        reasons.update({f"{species}: {k}": v for k, v in choice.reasons.items()})
    if not bands:
        msg = f"{gas!r} has no Raman-active constituent"
        raise DataError(msg)
    outside = [
        f"{b.species} {b.name} ({b.temperatures[0]:g}-{b.temperatures[1]:g} K)"
        for b in bands
        if not b.temperatures[0] <= temperature <= b.temperatures[1]
    ]
    if outside:
        msg = (
            f"{temperature} K is outside the temperatures of the width measurements "
            f"of {', '.join(outside)} (TC-688)"
        )
        warnings.warn(msg, ModelRangeWarning, stacklevel=2)
    for b in bands:
        sources.update(b.sources)
    status, reason = _status(pairs, choices)
    return replace(
        _assemble(pairs, temperature, wavelength, bands, sources, reasons),
        material=_material(gas, mixtures),
        population_floor=population_floor,
        status=status,
        status_reason=reason,
    )


def raman_response_status(
    gas: Gas | str | Mapping[str, float],
) -> tuple[Status, str] | None:
    """Return the status a gas's Raman response would have, without building it.

    The status and its reason are those of :func:`raman_response`'s result at any
    temperature and wavelength, read from the curated choices alone: a propagation
    code can ask whether a gas's response is validated or provisional, or whether it
    has none, without computing the lines. The collision partners are not checked: a
    mixture whose constituents lack broadening coefficients by one another (CO2 with
    N2) has a status here, and :func:`raman_response` refuses it.

    Parameters
    ----------
    gas : Gas or str or mapping of str to float
        As for :func:`raman_response`.

    Returns
    -------
    tuple of (str, str) or None
        The status (``validated``, ``provisional``, ...) and its reason (empty for
        ``validated``); None where no constituent is Raman-active, as for the noble
        gases, which have no Raman response at all.

    Raises
    ------
    DataError
        As :func:`raman_response` does for a molecule without curated choices or a
        composition it refuses.

    Examples
    --------
    >>> raman_response_status("N2")
    ('validated', '')
    >>> raman_response_status("CO2")[0]
    'provisional'
    >>> raman_response_status("Ar") is None
    True
    """
    choices, mixtures = _config()
    pairs, _, _ = _composition(gas, mixtures)
    for species, _ in pairs:
        gas_id = as_gas(species)
        if not info(gas_id).atomic and gas_id not in choices:
            msg = f"raman_response.toml has no choices for {species}"
            raise DataError(msg)
    if all(info(as_gas(s)).atomic for s, _ in pairs):
        return None
    return _status(pairs, choices)


def _status(
    pairs: tuple[tuple[str, float], ...], choices: Mapping[Gas, _GasChoices]
) -> tuple[Status, str]:
    """Return a response's status: that of its least trustworthy active constituent.

    Parameters
    ----------
    pairs : tuple
        The composition, species and mole fraction.
    choices : mapping
        The gases' choices.

    Returns
    -------
    tuple of (str, str)
        The status and its reason, the reasons of every constituent not validated.
    """
    active = [choices[as_gas(s)] for s, _ in pairs if not info(as_gas(s)).atomic]
    stopgaps = [c for c in active if c.status != "validated"]
    status: Status = max(  # pyright: ignore[reportAssignmentType] - statuses only
        (c.status for c in stopgaps), key=_RANK.index, default="validated"
    )
    return status, "; ".join(c.status_reason for c in stopgaps)


_RANK: Final = ("validated", "comparison", "provisional", "untested")
_AIR_ALIAS: Final = {"dry air": "dry air (N2, O2, Ar)"}


def _material(
    gas: Gas | str | Mapping[str, float], mixtures: Mapping[str, object]
) -> str:
    """Return the name of what a response is of.

    Parameters
    ----------
    gas : Gas or str or mapping of str to float
        What :func:`raman_response` was given.
    mixtures : mapping
        The mixtures of ``raman_response.toml``.

    Returns
    -------
    str
        A species, a named mixture, or the species of mole fractions joined by ``+``.
    """
    if isinstance(gas, Mapping):
        return "+".join(as_gas(k).value for k in gas)
    name = _AIR_ALIAS.get(gas, gas) if isinstance(gas, str) else gas
    if isinstance(name, str) and name in mixtures:
        return name
    return as_gas(name).value

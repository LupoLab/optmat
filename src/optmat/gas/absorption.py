# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The complex polarizability of a measured absorption spectrum (Kramers-Kronig).

A molecule that absorbs within the range of interest (ozone's Hartley, Huggins and
Chappuis bands) contributes to its polarizability a term fixed by the measured
photoabsorption cross section :math:`\sigma(\omega)`. For one molecule,

.. math::

   \mathrm{Im}\,\alpha(\omega) = \frac{\varepsilon_0 c\,\sigma(\omega)}{\omega},
   \qquad
   \mathrm{Re}\,\alpha(\omega) = \frac{2\varepsilon_0 c}{\pi}\,
   \mathcal{P}\!\int_0^\infty \frac{\sigma(\Omega)}{\Omega^2-\omega^2}\,\mathrm{d}\Omega
   = \frac{\varepsilon_0}{\pi^2}\,
   \mathcal{P}\!\int_0^\infty \frac{\sigma(\tilde\nu')}{\tilde\nu'^2-\tilde\nu^2}\,
   \mathrm{d}\tilde\nu',

in SI, with :math:`\tilde\nu` the wavenumber (the dilute-gas relation
:math:`n - 1 = N\alpha/2\varepsilon_0` gives the familiar form
:math:`(n-1)/N = (c/\pi)\,\mathcal{P}\!\int
\sigma/(\Omega^2-\omega^2)\,\mathrm{d}\Omega`).

The spectrum is assembled from pieces (``data/gas/absorption/<gas>.toml``), each a
measured data set over a wavelength window, on a uniform wavenumber grid, and smoothed
with a Gaussian of a stated width in wavenumber. Smoothing in frequency before the
transform is exactly causal: it is the same as smoothing the complex polarizability with
that kernel, or damping the molecular response in time beyond about ``0.88 / (c FWHM)``.
The principal value is a discrete convolution on the grid
(:math:`1/(\tilde\nu'^2-\tilde\nu^2) = [1/(\tilde\nu'-\tilde\nu) -
1/(\tilde\nu'+\tilde\nu)]/2\tilde\nu`, each part by FFT), accurate to about 1e-5 of the
term at the default grid and width.

The same spectrum gives its dipole sums :math:`S(k) = \int (\mathrm{d}f/\mathrm{d}E)
E^k\,\mathrm{d}E` (the oscillator-strength density is
:math:`\sigma / (\pi e^2\hbar/2\varepsilon_0 m_e c)`), so that the absorption term can
be checked against, and joined to, a dipole oscillator-strength distribution above the
measured range.
"""

import gzip
import itertools
import tomllib
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Final

import numpy as np
from jaxtyping import Float

from optmat import _catalogue
from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_ENERGY,
    AU_POLARIZABILITY,
    E_CHARGE,
    EPS0,
    HBAR,
    M_E,
    C,
)
from optmat.datasets import DATA_DIRECTORY
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.identity import Identified

if TYPE_CHECKING:
    from scipy.interpolate import CubicSpline

__all__ = [
    "ABSORPTION_DIRECTORY",
    "MeasuredAbsorption",
    "Piece",
    "assemble",
    "measured_absorption",
]

ABSORPTION_DIRECTORY: Final = GAS_DATA_DIRECTORY / "absorption"
"""Configurations of measured absorption spectra, ``<gas>.toml``."""

_EXTERNAL: Final = DATA_DIRECTORY / "external"
_FWHM_TO_SIGMA: Final = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
# Oscillator-strength density per unit energy (1/J) per unit cross section (m^2).
_DF_DE_PER_SIGMA: Final = 2.0 * EPS0 * M_E * C / (np.pi * E_CHARGE**2 * HBAR)


@dataclass(frozen=True)
class Piece:
    """One measured data set over a wavelength window.

    Attributes
    ----------
    source : str
        A file under ``data/external`` (two columns: vacuum wavelength in nm, cross
        section in cm^2) or a curated dataset ``source/dataset`` whose values are
        ``absorption_cross_section``.
    from_nm : float
        The window's short-wavelength end, nm.
    to_nm : float
        Its long-wavelength end, nm.
    note : str
        Why this data set covers this window.
    """

    source: str
    from_nm: float
    to_nm: float
    note: str = ""


def _read_file(relative: str) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column external spectrum.

    Parameters
    ----------
    relative : str
        Path under ``data/external`` (``.gz`` read transparently).

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelength (nm) and cross section (m^2), sorted by wavelength, duplicate
        wavelengths averaged.
    """
    path = _EXTERNAL / relative
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        rows = [
            line.split()[:2]
            for line in handle
            if line.strip() and not line.lstrip()[0].isalpha()
        ]
    table = np.array(rows, dtype=float)
    wavelength, inverse = np.unique(table[:, 0], return_inverse=True)
    sigma = np.bincount(inverse, weights=table[:, 1]) / np.bincount(inverse)
    return wavelength, sigma * 1e-4


def _read_dataset(key: str) -> tuple[np.ndarray, np.ndarray]:
    """Read a curated absorption dataset.

    Parameters
    ----------
    key : str
        ``source/dataset``.

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelength (nm) and cross section (m^2), sorted.

    Raises
    ------
    DataError
        If the dataset is unknown or not an absorption cross section in cm2 or Mb.
    """
    dataset = _catalogue.by_key().get(key)
    if dataset is None:
        msg = f"unknown absorption dataset {key!r}"
        raise DataError(msg)
    factor = {"cm2": 1e-4, "Mb": 1e-22}.get(dataset.unit)
    if dataset.quantity != "absorption_cross_section" or factor is None:
        msg = f"{key}: not an absorption_cross_section in cm2 or Mb"
        raise DataError(msg)
    wavelength = np.asarray(dataset.vacuum_wavelength(), dtype=float) * 1e9
    sigma = np.asarray(dataset.scaled_values(), dtype=float) * factor
    order = np.argsort(wavelength)
    return wavelength[order], sigma[order]


def _read(source: str) -> tuple[np.ndarray, np.ndarray]:
    """Read a piece's data set, from a file or a curated dataset.

    Parameters
    ----------
    source : str
        See :attr:`Piece.source`.

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelength (nm) and cross section (m^2).
    """
    if (_EXTERNAL / source).is_file():
        return _read_file(source)
    return _read_dataset(source)


@dataclass(frozen=True, eq=False)
class MeasuredAbsorption(Identified):
    """A measured absorption spectrum and its complex polarizability.

    It is compared by identity, not by value (``eq=False``), since its spline is cached
    per object.

    Attributes
    ----------
    gas : Gas
        The species.
    temperature_k : float
        The temperature of the measurements, K.
    pieces : tuple of Piece
        The data sets and their windows, from long to short wavelength or in any order;
        the windows must not overlap.
    fwhm_cm : float
        The Gaussian smoothing, full width at half maximum in wavenumber, cm^-1 (0 for
        none beyond the grid).
    step_cm : float
        The wavenumber grid step, cm^-1.
    wavenumber_cm : numpy.ndarray
        The grid, cm^-1, from 0 to the largest wavenumber of the pieces.
    sigma_m2 : numpy.ndarray
        The smoothed cross section on the grid, m^2.
    real_au : numpy.ndarray
        ``Re alpha`` on the grid, atomic units.
    """

    gas: Gas
    temperature_k: float
    pieces: tuple[Piece, ...]
    fwhm_cm: float
    step_cm: float
    wavenumber_cm: Float[np.ndarray, " n"] = field(repr=False)
    sigma_m2: Float[np.ndarray, " n"] = field(repr=False)
    real_au: Float[np.ndarray, " n"] = field(repr=False)

    IDENTITY_KIND: ClassVar[str] = "gas_absorption"
    EVALUATION_VERSION: ClassVar[int] = 1

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the spectrum's name, gas and data sets.

        Returns
        -------
        tuple
            ``("absorption/<gas>", gas, sources)``, the sources being the pieces'
            curated datasets' records (an external file is named by its path).
        """
        sources = tuple(sorted({p.source.split("/", 1)[0] for p in self.pieces}))
        return f"absorption/{self.gas.value}", self.gas.value, sources

    @property
    def max_wavenumber_cm(self) -> float:
        """The top of the measured range, cm^-1.

        Returns
        -------
        float
            The largest grid wavenumber.
        """
        return float(self.wavenumber_cm[-1])

    def complex_polarizability_au(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> tuple[Float[np.ndarray, "*shape"], Float[np.ndarray, "*shape"]]:
        """Return the real and imaginary polarizability of the measured term.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree; above the measured range the real part is the
            regular integral and the imaginary part zero.

        Returns
        -------
        tuple of numpy.ndarray
            ``Re alpha`` and ``Im alpha``, atomic units.
        """
        omega = np.asarray(omega, dtype=float)
        nu = omega * AU_ANGULAR_FREQUENCY / (2 * np.pi * C) / 100  # cm^-1
        spline = _spline(self)
        inside = nu <= self.max_wavenumber_cm
        real = np.where(inside, spline(np.minimum(nu, self.max_wavenumber_cm)), 0.0)
        if np.any(~inside):
            outside = nu[~inside]
            grid = self.wavenumber_cm * 100
            weights = self.sigma_m2[:, None] / (
                grid[:, None] ** 2 - (outside * 100) ** 2
            )
            integral = np.trapezoid(weights, grid, axis=0)
            real = real.copy()
            real[~inside] = EPS0 / np.pi**2 * integral / AU_POLARIZABILITY
        sigma = np.interp(nu, self.wavenumber_cm, self.sigma_m2, right=0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            imaginary = np.where(
                omega > 0,
                EPS0 * C * sigma / (omega * AU_ANGULAR_FREQUENCY) / AU_POLARIZABILITY,
                0.0,
            )
        return np.asarray(real), np.asarray(imaginary)

    def complex_derivative_au(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> tuple[Float[np.ndarray, "*shape"], Float[np.ndarray, "*shape"]]:
        """Return the frequency derivative of the measured term, branch by branch.

        The derivative is that of the evaluator :meth:`complex_polarizability_au` uses
        (OEP 1, section 4.1):

        - inside the grid, the real part is a cubic spline, and its derivative is the
          spline's own;
        - above the grid, the real part is a trapezoid sum over the grid, and its
          derivative is the same sum of the differentiated integrand;
        - the imaginary part is ``eps0 c sigma(nu) / omega`` with ``sigma`` linearly
          interpolated, so its derivative takes the slope of the segment: the mean of
          the two slopes at a knot, the segment below at the grid's upper end, zero
          above it.

        Both are 0 at ``omega = 0``, by definition: the real part is even in
        ``omega`` and the imaginary part vanishes there.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.

        Returns
        -------
        tuple of numpy.ndarray
            ``d Re alpha / d omega`` and ``d Im alpha / d omega``, atomic units per
            hartree.
        """
        omega = np.asarray(omega, dtype=float)
        per_hartree = AU_ANGULAR_FREQUENCY / (2 * np.pi * C) / 100  # cm^-1 per hartree
        nu = omega * per_hartree
        top = self.max_wavenumber_cm
        inside = nu <= top
        real = np.where(
            inside, _spline(self)(np.minimum(nu, top), 1) * per_hartree, 0.0
        )
        if np.any(~inside):
            outside = nu[~inside] * 100  # m^-1
            grid = self.wavenumber_cm * 100
            weights = (
                self.sigma_m2[:, None]
                * 2
                * outside
                / (grid[:, None] ** 2 - outside**2) ** 2
            )
            integral = np.trapezoid(weights, grid, axis=0)
            real = real.copy()
            real[~inside] = (
                EPS0 / np.pi**2 * integral / AU_POLARIZABILITY * per_hartree * 100
            )
        grid, sigma = self.wavenumber_cm, self.sigma_m2
        slopes = np.diff(sigma) / np.diff(grid)
        k = np.clip(np.searchsorted(grid, nu, side="right") - 1, 0, len(slopes) - 1)
        at_knot = (nu == grid[k]) & (k > 0)
        slope = np.where(at_knot, 0.5 * (slopes[k - 1] + slopes[k]), slopes[k])
        slope = np.where(nu == top, slopes[-1], slope)
        slope = np.where(nu > top, 0.0, slope)
        value = np.interp(nu, grid, sigma, right=0.0)
        scale = EPS0 * C / (AU_ANGULAR_FREQUENCY * AU_POLARIZABILITY)
        with np.errstate(divide="ignore", invalid="ignore"):
            imaginary = np.where(
                omega > 0,
                scale * (slope * per_hartree / omega - value / omega**2),
                0.0,
            )
        return np.asarray(np.where(omega > 0, real, 0.0)), np.asarray(imaginary)

    def real_derivative_au(
        self, omega: Float[np.ndarray, "*shape"], order: int
    ) -> Float[np.ndarray, "*shape"]:
        """Return a derivative of ``Re alpha`` with respect to the photon energy.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree, within the grid.
        order : int
            The order, 0 to 4 (finite differences of the smoothed grid values, then
            linear interpolation).

        Returns
        -------
        numpy.ndarray
            ``d^order Re alpha / d omega^order``, atomic units per hartree^order.
        """
        omega = np.asarray(omega, dtype=float)
        per_hartree = AU_ANGULAR_FREQUENCY / (2 * np.pi * C) / 100  # cm^-1 per hartree
        nu = omega * per_hartree
        values = _derivatives(self)[order]
        return np.asarray(
            np.interp(nu, self.wavenumber_cm, values) * per_hartree**order
        )

    def sum_rule(self, k: float) -> float:
        """Return the dipole sum ``S(k)`` of the measured spectrum.

        Parameters
        ----------
        k : float
            The order.

        Returns
        -------
        float
            ``int (df/dE) E^k dE`` over the grid, atomic units.
        """
        energy = self.wavenumber_cm * 100 * 2 * np.pi * C * HBAR  # J
        density = self.sigma_m2 * _DF_DE_PER_SIGMA * AU_ENERGY  # per hartree
        energy_au = energy / AU_ENERGY
        positive = energy_au > 0
        return float(
            np.trapezoid(
                density[positive] * energy_au[positive] ** k, energy_au[positive]
            )
        )

    def strength(self, low_ev: float, high_ev: float) -> float:
        """Return the oscillator strength between two photon energies.

        Parameters
        ----------
        low_ev, high_ev : float
            The interval, eV.

        Returns
        -------
        float
            ``int df/dE dE`` over the interval.
        """
        energy_ev = self.wavenumber_cm * 100 * 2 * np.pi * C * HBAR / E_CHARGE
        inside = (energy_ev >= low_ev) & (energy_ev <= high_ev)
        density = self.sigma_m2 * _DF_DE_PER_SIGMA * E_CHARGE  # per eV
        return float(np.trapezoid(density[inside], energy_ev[inside]))


def assemble(
    gas: Gas | str,
    pieces: tuple[Piece, ...],
    *,
    temperature_k: float,
    fwhm_cm: float = 50.0,
    step_cm: float = 0.5,
) -> MeasuredAbsorption:
    """Assemble, smooth and transform a spectrum.

    Parameters
    ----------
    gas : Gas or str
        The species.
    pieces : tuple of Piece
        The data sets and their windows.
    temperature_k : float
        The temperature of the measurements, K.
    fwhm_cm : float, optional
        Gaussian smoothing FWHM in wavenumber, cm^-1; 50 by default.
    step_cm : float, optional
        The grid step, cm^-1; 0.5 by default.

    Returns
    -------
    MeasuredAbsorption
        The spectrum with its real polarizability on the grid.

    Raises
    ------
    DataError
        If two windows overlap.
    """
    windows = sorted((p.from_nm, p.to_nm) for p in pieces)
    for (_, high), (low, _) in itertools.pairwise(windows):
        if low < high:
            msg = f"absorption pieces overlap near {low} nm"
            raise DataError(msg)
    top = 1e7 / min(p.from_nm for p in pieces)
    grid = np.arange(0.0, top + step_cm, step_cm)
    sigma = np.zeros_like(grid)
    for piece in pieces:
        wavelength, values = _read(piece.source)
        nu = 1e7 / wavelength[::-1]
        inside = (grid >= 1e7 / piece.to_nm) & (grid <= 1e7 / piece.from_nm)
        sigma[inside] = np.interp(grid[inside], nu, values[::-1])
    if fwhm_cm > 0:
        # Imported on first use: scipy.ndimage is slow to import.
        from scipy.ndimage import gaussian_filter1d  # noqa: PLC0415

        sigma = gaussian_filter1d(
            sigma, fwhm_cm * _FWHM_TO_SIGMA / step_cm, mode="constant"
        )
    sigma = np.maximum(sigma, 0.0)
    return MeasuredAbsorption(
        gas=as_gas(gas),
        temperature_k=temperature_k,
        pieces=tuple(pieces),
        fwhm_cm=fwhm_cm,
        step_cm=step_cm,
        wavenumber_cm=grid,
        sigma_m2=sigma,
        real_au=_kramers_kronig(grid, sigma, step_cm),
    )


def _kramers_kronig(
    grid_cm: Float[np.ndarray, " n"], sigma: Float[np.ndarray, " n"], step_cm: float
) -> Float[np.ndarray, " n"]:
    """Return ``Re alpha`` (atomic units) on a uniform wavenumber grid from 0.

    Parameters
    ----------
    grid_cm : numpy.ndarray
        Wavenumbers ``m * step``, cm^-1, starting at 0.
    sigma : numpy.ndarray
        Cross section on the grid, m^2.
    step_cm : float
        The grid step, cm^-1.

    Returns
    -------
    numpy.ndarray
        ``(eps0 / pi^2) PV int sigma / (nu'^2 - nu^2) dnu'`` in atomic units; at
        ``nu = 0`` the static value ``int sigma / nu'^2``.
    """
    # Imported on first use: scipy.signal is slow to import (it brings scipy.stats).
    from scipy.signal import fftconvolve  # noqa: PLC0415

    size = len(grid_cm)
    step = step_cm * 100.0  # m^-1
    m = np.arange(-(size - 1), size)
    safe = np.where(m == 0, 1, m)
    difference = np.where(m == 0, 0.0, 1.0 / (safe * step))
    # sum_j sigma_j / (nu_j - nu_i): the principal value, the singular term dropped.
    minus = fftconvolve(sigma, difference[::-1])[size - 1 : 2 * size - 1] * step
    index = np.arange(2 * size)
    total = np.where(index == 0, 0.0, 1.0 / (np.where(index == 0, 1, index) * step))
    # sum_j sigma_j / (nu_j + nu_i) = (reversed sigma * 1/((i + j) step))[i + size - 1]
    plus = fftconvolve(sigma[::-1], total)[size - 1 : 2 * size - 1] * step
    nu = grid_cm * 100.0
    with np.errstate(divide="ignore", invalid="ignore"):
        integral = (minus - plus) / (2.0 * nu)
    positive = nu > 0
    static = float(np.sum(sigma[positive] / nu[positive] ** 2) * step)
    integral = np.where(positive, integral, static)
    return EPS0 / np.pi**2 * integral / AU_POLARIZABILITY


@cache
def _derivatives(absorption: MeasuredAbsorption) -> tuple[np.ndarray, ...]:
    """Return ``Re alpha`` and its first four wavenumber derivatives on the grid.

    Parameters
    ----------
    absorption : MeasuredAbsorption
        The spectrum.

    Returns
    -------
    tuple of numpy.ndarray
        Derivatives of order 0 to 4, atomic units per (cm^-1)^order.
    """
    values = [absorption.real_au]
    for _ in range(4):
        values.append(np.gradient(values[-1], absorption.step_cm))
    return tuple(values)


@cache
def _spline(absorption: MeasuredAbsorption) -> "CubicSpline":
    """Return a cubic spline of the real part on the grid (cached per spectrum).

    Parameters
    ----------
    absorption : MeasuredAbsorption
        The spectrum.

    Returns
    -------
    scipy.interpolate.CubicSpline
        ``Re alpha`` against wavenumber, cm^-1.
    """
    # Imported on first use: scipy.interpolate is slow to import.
    from scipy.interpolate import CubicSpline  # noqa: PLC0415

    return CubicSpline(absorption.wavenumber_cm, absorption.real_au)


def _pieces(config: dict[str, object], path: Path) -> tuple[Piece, ...]:
    """Read the ``[[pieces]]`` of a configuration.

    Parameters
    ----------
    config : dict
        The parsed TOML.
    path : pathlib.Path
        For error messages.

    Returns
    -------
    tuple of Piece
        The pieces.

    Raises
    ------
    DataError
        If the pieces are not an array of tables or a piece lacks a key.
    """
    entries = config.get("pieces", [])
    if not isinstance(entries, list):
        msg = f"{path.name}: pieces must be an array of tables"
        raise DataError(msg)
    pieces = []
    for entry in entries:
        try:
            pieces.append(
                Piece(
                    source=str(entry["source"]),
                    from_nm=float(entry["from_nm"]),
                    to_nm=float(entry["to_nm"]),
                    note=str(entry.get("note", "")),
                )
            )
        except KeyError as error:
            msg = f"{path.name}: a piece lacks {error}"
            raise DataError(msg) from None
    return tuple(pieces)


@cache
def measured_absorption(
    gas: Gas | str,
    *,
    fwhm_cm: float | None = None,
    variant: str | None = None,
) -> MeasuredAbsorption:
    """Return a gas's measured absorption spectrum (``data/gas/absorption/<gas>.toml``).

    Parameters
    ----------
    gas : Gas or str
        The species.
    fwhm_cm : float, optional
        The smoothing FWHM, cm^-1; the configuration's by default.
    variant : str, optional
        A ``[variants.<name>]`` table of the configuration whose pieces replace those of
        the same window (for instance another temperature).

    Returns
    -------
    MeasuredAbsorption
        The spectrum.

    Raises
    ------
    DataError
        If the gas has no configuration or the variant is unknown.
    """
    gas = as_gas(gas)
    path = ABSORPTION_DIRECTORY / f"{gas.value}.toml"
    if not path.is_file():
        msg = f"no measured absorption spectrum for {gas.value}"
        raise DataError(msg)
    config = tomllib.loads(path.read_text(encoding="utf-8"))
    pieces = _pieces(config, path)
    temperature = float(config["temperature_k"])
    if variant is not None:
        variants = config.get("variants", {})
        if variant not in variants:
            msg = f"{path.name}: no variant {variant!r}"
            raise DataError(msg)
        replacement = _pieces(variants[variant], path)
        replaced = {(p.from_nm, p.to_nm) for p in replacement}
        pieces = (
            tuple(p for p in pieces if (p.from_nm, p.to_nm) not in replaced)
            + replacement
        )
        temperature = float(variants[variant].get("temperature_k", temperature))
    return assemble(
        gas,
        pieces,
        temperature_k=temperature,
        fwhm_cm=float(config.get("fwhm_cm", 50.0)) if fwhm_cm is None else fwhm_cm,
        step_cm=float(config.get("step_cm", 0.5)),
    )

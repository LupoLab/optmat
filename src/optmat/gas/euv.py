# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Build the complex polarizability from photoabsorption (the extreme ultraviolet).

Above its first resonance an atom absorbs, and its polarizability is complex. Both parts
follow from the dipole oscillator-strength distribution, the discrete lines
:math:`(E_i, f_i)` and the continuum density
:math:`g(E) = \mathrm{d}f/\mathrm{d}E` above the ionisation threshold :math:`I`:

.. math::

   \mathrm{Re}\,\alpha(\omega) = \sum_i \frac{f_i}{E_i^2 - \omega^2}
   + \mathcal{P}\!\int_I^\infty \frac{g(E)}{E^2 - \omega^2}\,\mathrm{d}E,
   \qquad
   \mathrm{Im}\,\alpha(\omega) = \frac{\pi}{2\omega}\,g(\omega),

in atomic units (the Kramers-Kronig pair of the causal response; the lines contribute
delta functions to the imaginary part, which are left out). The continuum is measured
as a photoionisation cross section, :math:`\sigma = (\pi e^2\hbar/2\varepsilon_0 m_e c)
\,g(E)`, 109.76 Mb per eV\ :sup:`-1`. The same distribution gives the dipole sums
:math:`S(k) = \sum_i f_i E_i^k + \int g(E) E^k\,\mathrm{d}E`: the Thomas-Reiche-Kuhn sum
:math:`S(0)` is the number of electrons and :math:`S(-2)` the static polarizability,
which test the data.

The principal value is taken by subtracting :math:`g(\omega)` from the integrand, whose
remainder is regular at :math:`E = \omega`, and adding
:math:`g(\omega)\,\mathcal{P}\!\int_I^{E_\mathrm{max}} \mathrm{d}E/(E^2-\omega^2)
= (g(\omega)/2\omega) \ln[(E_\mathrm{max}-\omega)(I+\omega)/((E_\mathrm{max}+\omega)
(\omega-I))]` analytically. Beyond the last tabulated energy the continuum falls as
:math:`E^{-7/2}`, the nonrelativistic asymptote.

The refractive index of the dilute gas is :math:`n = 1 + N\alpha/2\varepsilon_0`: above
threshold :math:`\mathrm{Re}\,n - 1` becomes negative and :math:`\mathrm{Im}\,n` is the
absorption index. For comparison, Henke's atomic scattering factors give
:math:`\alpha = -(f_1 - i f_2)/\omega^2` in atomic units: Henke et al. write
:math:`n = 1 - \delta - i\beta` with :math:`\beta \propto f_2 > 0`, the opposite sign
convention to the package's :math:`\mathrm{Im}\,n > 0` for absorption.
"""

import csv
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from jaxtyping import Complex, Float

from optmat import _catalogue
from optmat.conditions import ONE_ATM_0C, State
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
from optmat.gas.eos import number_density
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.gases import Gas, as_gas

__all__ = [
    "EV_PER_HARTREE",
    "MB_PER_INVERSE_EV",
    "ComplexModel",
    "OscillatorDistribution",
    "constrained",
    "euv_model",
    "nist_resonance_lines",
]

EV_PER_HARTREE: Final = AU_ENERGY / E_CHARGE
"""Electronvolts per hartree."""

MB_PER_INVERSE_EV: Final = (
    np.pi * E_CHARGE**2 * HBAR / (2 * EPS0 * M_E * C) / E_CHARGE / 1e-22
)
"""Photoionisation cross section (Mb) per unit ``df/dE`` (1/eV), 109.761:
``pi e^2 hbar / (2 eps0 m_e c)``, per eV and in units of 1e-22 m^2."""

_TAIL_EXPONENT: Final = 3.5


@dataclass(frozen=True)
class OscillatorDistribution:
    """A dipole oscillator-strength distribution: discrete lines and a continuum.

    Attributes
    ----------
    gas : Gas
        The species.
    name : str
        What the distribution is built from.
    line_energies : numpy.ndarray
        Discrete line energies, hartree.
    line_strengths : numpy.ndarray
        Their oscillator strengths.
    continuum_energies : numpy.ndarray
        Photon energies of the continuum samples, hartree, increasing from the
        threshold.
    continuum_density : numpy.ndarray
        ``df/dE`` at those energies, per hartree.
    """

    gas: Gas
    name: str
    line_energies: Float[np.ndarray, " n_lines"] = field(repr=False)
    line_strengths: Float[np.ndarray, " n_lines"] = field(repr=False)
    continuum_energies: Float[np.ndarray, " n_continuum"] = field(repr=False)
    continuum_density: Float[np.ndarray, " n_continuum"] = field(repr=False)

    def __post_init__(self) -> None:
        """Check the tables.

        Raises
        ------
        ValueError
            If the continuum energies do not increase or the tables disagree in length.
        """
        if np.any(np.diff(self.continuum_energies) <= 0):
            msg = f"{self.name}: continuum energies must increase"
            raise ValueError(msg)
        if len(self.line_energies) != len(self.line_strengths) or len(
            self.continuum_energies
        ) != len(self.continuum_density):
            msg = f"{self.name}: tables of unequal length"
            raise ValueError(msg)

    @property
    def threshold(self) -> float:
        """The ionisation threshold, the first continuum energy, hartree.

        Returns
        -------
        float
            The threshold.
        """
        return float(self.continuum_energies[0])

    def _grid(self, points: int = 20000) -> tuple[np.ndarray, np.ndarray]:
        """Return a fine logarithmic grid over the continuum and its tail, with ``g``.

        Parameters
        ----------
        points : int, optional
            The number of grid points.

        Returns
        -------
        tuple of numpy.ndarray
            Energies (hartree) and ``df/dE`` there: interpolated linearly in log-log
            inside the table, the ``E^-3.5`` tail beyond it (to 1e4 times the last
            energy).
        """
        energies = np.geomspace(
            self.threshold, self.continuum_energies[-1] * 1e4, points
        )
        inside = energies <= self.continuum_energies[-1]
        density = np.empty_like(energies)
        density[inside] = np.exp(
            np.interp(
                np.log(energies[inside]),
                np.log(self.continuum_energies),
                np.log(self.continuum_density),
            )
        )
        last_e, last_g = self.continuum_energies[-1], self.continuum_density[-1]
        density[~inside] = last_g * (energies[~inside] / last_e) ** -_TAIL_EXPONENT
        return energies, density

    def sum_rule(self, k: float) -> float:
        """Return the dipole sum ``S(k)``.

        Parameters
        ----------
        k : float
            The order; the integral converges for ``k < 2.5`` with the ``E^-3.5`` tail.

        Returns
        -------
        float
            ``sum_i f_i E_i^k + int g(E) E^k dE``, atomic units.
        """
        energies, density = self._grid()
        tail_start = energies[-1]
        # The part of the tail beyond the grid, analytically.
        tail = density[-1] * tail_start ** (k + 1) / (_TAIL_EXPONENT - k - 1)
        continuum = float(np.trapezoid(density * energies**k, energies)) + tail
        return float(np.sum(self.line_strengths * self.line_energies**k)) + continuum

    def complex_polarizability_au(
        self, omega: Float[np.ndarray, " n"]
    ) -> tuple[Float[np.ndarray, " n"], Float[np.ndarray, " n"]]:
        """Return the real and imaginary parts of the polarizability.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree; exactly at a discrete line the real part is
            infinite.

        Returns
        -------
        tuple of numpy.ndarray
            ``Re alpha`` and ``Im alpha``, atomic units.
        """
        omega = np.atleast_1d(np.asarray(omega, dtype=float))
        energies, density = self._grid()
        e_max, threshold = energies[-1], self.threshold
        lines = np.sum(
            self.line_strengths[:, None]
            / (self.line_energies[:, None] ** 2 - omega[None, :] ** 2),
            axis=0,
        )
        real = np.empty_like(omega)
        imaginary = np.zeros_like(omega)
        for index, w in enumerate(omega):
            above = threshold < w < e_max
            g_w = float(np.interp(w, energies, density)) if above else 0.0
            numerator = density - g_w
            denominator = energies**2 - w**2
            with np.errstate(divide="ignore", invalid="ignore"):
                integrand = np.where(denominator != 0, numerator / denominator, 0.0)
            integral = float(np.trapezoid(integrand, energies))
            if above:
                integral += (
                    g_w
                    / (2 * w)
                    * np.log(
                        (e_max - w) * (threshold + w) / ((e_max + w) * (w - threshold))
                    )
                )
                imaginary[index] = np.pi / (2 * w) * g_w
            real[index] = lines[index] + integral
        return real, imaginary


def constrained(
    distribution: OscillatorDistribution,
    region_max: float,
    target: tuple[float, float],
) -> tuple[OscillatorDistribution, float]:
    """Scale the continuum below an energy so that one dipole sum takes a target value.

    Typically the measured continuum near threshold is scaled, within its stated
    uncertainty, to the static polarizability ``S(-2)`` of an ab initio calculation, the
    simplest of Kumar and Meath's constrained distributions; the higher-energy data are
    left as they are, and the Thomas-Reiche-Kuhn sum ``S(0)`` is kept as a check.

    Parameters
    ----------
    distribution : OscillatorDistribution
        The measured distribution.
    region_max : float
        The continuum below this energy (hartree) is scaled.
    target : tuple of float
        ``(k, value)``: the sum ``S(k)`` to meet.

    Returns
    -------
    tuple
        The scaled distribution and the factor applied.

    Raises
    ------
    DataError
        If the target does not respond to scaling, or a secant step requires a
        non-positive or non-finite continuum scale. The local iteration is intended
        for adjustments within the measured uncertainty, not arbitrary targets.
    """
    k, value = target
    below = distribution.continuum_energies < region_max

    def scaled_by(factor: float) -> OscillatorDistribution:
        """Return the distribution with the region's continuum times ``factor``.

        Parameters
        ----------
        factor : float
            The scale.

        Returns
        -------
        OscillatorDistribution
            The scaled distribution.
        """
        return OscillatorDistribution(
            gas=distribution.gas,
            name=f"{distribution.name}, S({k:g}) constrained",
            line_energies=distribution.line_energies,
            line_strengths=distribution.line_strengths,
            continuum_energies=distribution.continuum_energies,
            continuum_density=distribution.continuum_density
            * np.where(below, factor, 1.0),
        )

    # S(k) is nearly linear in the factor (exactly, but for the log-log interpolation
    # across the split): a few secant steps meet the target to round-off.
    factors = [1.0, 1.01]
    residuals = [scaled_by(f).sum_rule(k) - value for f in factors]
    for _ in range(6):
        if abs(residuals[-1]) <= 1e-14 * abs(value) or factors[-1] == factors[-2]:
            break
        slope = (residuals[-1] - residuals[-2]) / (factors[-1] - factors[-2])
        if slope == 0:
            msg = "the constraint target does not respond to continuum scaling"
            raise DataError(msg)
        factor = factors[-1] - residuals[-1] / slope
        # Log-log interpolation requires a positive density at every sample (TC-736).
        if not np.isfinite(factor) or factor <= 0:
            msg = "the constraint iteration proposed an invalid continuum scale"
            raise DataError(msg)
        factors.append(factor)
        residuals.append(scaled_by(factors[-1]).sum_rule(k) - value)
    return scaled_by(factors[-1]), float(factors[-1])


def nist_resonance_lines(gas: Gas) -> list[tuple[float, float, str]]:
    """Return the ground-state J = 1 lines with oscillator strengths from the NIST file.

    Parameters
    ----------
    gas : Gas
        A noble gas with a ``data/external/nist_asd/<gas>_I_lines.tsv`` file.

    Returns
    -------
    list of tuple
        ``(energy in hartree, f, upper configuration)``, by increasing energy.
    """
    path = DATA_DIRECTORY / "external" / "nist_asd" / f"{gas.value}_I_lines.tsv"
    lines = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            clean = {k: (v or "").strip('"[]+ ') for k, v in row.items() if k}
            if not clean.get("fik") or clean.get("J_k") != "1":
                continue
            if float(clean["Ei(cm-1)"] or "nan") != 0.0:
                continue
            wavelength = clean["ritz_wl_vac(nm)"] or clean["obs_wl_vac(nm)"]
            lines.append(
                (HARTREE_NM / float(wavelength), float(clean["fik"]), clean["conf_k"])
            )
    return sorted(lines)


def _rydberg_continuation(
    last: tuple[float, float], threshold: float, n_last: int, n_max: int
) -> list[tuple[float, float]]:
    """Continue a Rydberg series from its last listed member to ``n_max``.

    Parameters
    ----------
    last : tuple of float
        The last listed line's energy (hartree) and oscillator strength.
    threshold : float
        The series limit, hartree.
    n_last : int
        The last listed line's principal quantum number.
    n_max : int
        The last member to add.

    Returns
    -------
    list of tuple
        ``(energy, f)`` of the added members: constant quantum defect, ``f`` scaling as
        the inverse cube of the effective quantum number.
    """
    energy, strength = last
    n_star = 1.0 / np.sqrt(2.0 * (threshold - energy))
    added = []
    for n in range(n_last + 1, n_max + 1):
        n_eff = n_star + (n - n_last)
        added.append((threshold - 0.5 / n_eff**2, strength * (n_star / n_eff) ** 3))
    return added


def _continuum_part(
    entry: Mapping[str, object], gas: Gas
) -> tuple[np.ndarray, np.ndarray]:
    """Return one continuum source's energies (hartree) and ``df/dE`` (per hartree).

    Parameters
    ----------
    entry : mapping
        ``{dataset = key}`` (a photoionisation cross section in Mb or a
        ``differential_oscillator_strength`` in 1/eV, photon energies in eV) or ``{henke
        = element}``, with optional ``min_ev`` and ``max_ev``, and ``extend_to_ev`` to
        carry the first tabulated density down to a lower ionisation limit.
    gas : Gas
        The species.

    Returns
    -------
    tuple of numpy.ndarray
        Energies and densities.

    Raises
    ------
    DataError
        If the dataset's quantity or units are not a continuum.
    """
    low = float(entry.get("min_ev", 0.0))  # pyright: ignore[reportArgumentType]
    high = float(entry.get("max_ev", np.inf))  # pyright: ignore[reportArgumentType]
    if "henke" in entry:
        table = np.loadtxt(
            DATA_DIRECTORY / "external" / "henke" / f"{entry['henke']}.nff", skiprows=1
        )
        energy_ev, f2 = table[:, 0], table[:, 2]
        energies = energy_ev / EV_PER_HARTREE
        density = (
            2.0 * f2 / (np.pi * energies)
        )  # Im alpha = f2/omega^2 = pi g/(2 omega)
    else:
        dataset = _catalogue.by_key(gas)[str(entry["dataset"])]
        if dataset.abscissa_unit != "eV":
            msg = f"{dataset.key}: a continuum needs photon energies in eV"
            raise DataError(msg)
        energy_ev = dataset.abscissa_values()
        values = dataset.scaled_values()
        if dataset.quantity == "photoionization_cross_section" and dataset.unit == "Mb":
            per_ev = values / MB_PER_INVERSE_EV
        elif (
            dataset.quantity == "differential_oscillator_strength"
            and dataset.unit == "1/eV"
        ):
            per_ev = values
        else:
            msg = f"{dataset.key}: not a continuum in Mb or 1/eV"
            raise DataError(msg)
        energies = energy_ev / EV_PER_HARTREE
        density = per_ev * EV_PER_HARTREE
    keep = (energy_ev >= low) & (energy_ev <= high)
    energies, density = energies[keep], density[keep]
    if "extend_to_ev" in entry:
        # The autoionising region between the two ionisation limits of a heavier noble
        # gas is not tabulated: its density is taken as the first tabulated value, down
        # to the lower limit (a stated approximation, checked by the sum rules).
        start = float(entry["extend_to_ev"]) / EV_PER_HARTREE  # pyright: ignore[reportArgumentType]
        energies = np.concatenate([[start], energies])
        density = np.concatenate([[density[0]], density])
    return energies, density


@dataclass(frozen=True)
class ComplexModel:
    """A constrained oscillator-strength distribution as a complex-index model.

    Attributes
    ----------
    gas : Gas
        The species.
    name : str
        ``euv/<gas>``.
    distribution : OscillatorDistribution
        The constrained distribution.
    chi0_au : float
        The magnetic susceptibility added to the electric polarizability.
    factor : float
        The scale applied to the measured continuum by the constraint.
    """

    gas: Gas
    name: str
    distribution: OscillatorDistribution = field(repr=False)
    chi0_au: float = 0.0
    factor: float = 1.0

    def polarizability_au(
        self, omega: Float[np.ndarray, " n"]
    ) -> Complex[np.ndarray, " n"]:
        """Return the complex refractive polarizability, atomic units.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.

        Returns
        -------
        numpy.ndarray
            ``Re alpha + chi_0 + i Im alpha`` (complex).
        """
        real, imaginary = self.distribution.complex_polarizability_au(omega)
        return (real + self.chi0_au) + 1j * imaginary

    def refractive_index(
        self,
        wavelength: Float[np.ndarray, "*shape"] | float,
        state: State = ONE_ATM_0C,
    ) -> Complex[np.ndarray, "*shape"]:
        """Return the complex refractive index at a state (Lorentz-Lorenz).

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State, optional
            Pressure and temperature; 0 degrees C and 1 atm by default.

        Returns
        -------
        numpy.ndarray
            ``n`` (complex; ``Im n > 0`` is absorption), the shape of ``wavelength``.
        """
        shape = np.shape(wavelength)
        flat = np.atleast_1d(np.asarray(wavelength, dtype=float)).ravel()
        omega = 2 * np.pi * C / flat / AU_ANGULAR_FREQUENCY
        alpha = self.polarizability_au(omega) * AU_POLARIZABILITY
        y = number_density(self.gas, state) * alpha / (3 * EPS0)
        return np.sqrt(1 + 3 * y / (1 - y)).reshape(shape)


def _lines(config: Mapping[str, object], gas: Gas) -> list[tuple[float, float]]:
    """Return the discrete lines a configuration names.

    Parameters
    ----------
    config : mapping
        The ``[lines]`` table: ``dataset`` (oscillator strengths against photon
        energies in eV) or, without it, the NIST lines; either optionally with ``total``
        (a dataset whose value is the summed strength to the ionisation limit) and
        ``pseudo_line_ev`` (where the unlisted remainder is placed).
    gas : Gas
        The species.

    Returns
    -------
    list of tuple
        ``(energy in hartree, f)``.
    """
    by_key = _catalogue.by_key(gas)
    if "dataset" in config:
        listed = by_key[str(config["dataset"])]
        energies = listed.abscissa_values() / EV_PER_HARTREE
        lines = list(zip(energies, listed.scaled_values(), strict=True))
    else:
        lines = [(e, f) for e, f, _ in nist_resonance_lines(gas)]
    if "total" in config:
        total = float(by_key[str(config["total"])].scaled_values()[0])
        remainder = total - sum(f for _, f in lines)
        if remainder > 0:
            pseudo = float(config["pseudo_line_ev"]) / EV_PER_HARTREE  # pyright: ignore[reportArgumentType]
            lines.append((pseudo, remainder))
    return sorted((float(e), float(f)) for e, f in lines)


def euv_model(gas: Gas | str) -> ComplexModel:
    """Build a gas's complex-index model from ``data/gas/euv/<gas>.toml``.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    ComplexModel
        The constrained distribution as a model.

    Raises
    ------
    DataError
        If the configuration is malformed or refers to missing data.
    """
    gas = as_gas(gas)
    path = GAS_DATA_DIRECTORY / "euv" / f"{gas.value}.toml"
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    line_list = _lines(config.get("lines", {}), gas)
    parts = [_continuum_part(entry, gas) for entry in config["continuum"]]
    energies = np.concatenate([p[0] for p in parts])
    density = np.concatenate([p[1] for p in parts])
    order = np.argsort(energies)
    energies, density = energies[order], density[order]
    rydberg = config.get("lines", {}).get("rydberg_to_n")
    if rydberg and "dataset" not in config.get("lines", {}):
        n_last = int(config["lines"]["last_listed_n"])
        line_list += _rydberg_continuation(line_list[-1], energies[0], n_last, rydberg)
    distribution = OscillatorDistribution(
        gas=gas,
        name=f"euv/{gas.value}",
        line_energies=np.array([e for e, _ in line_list]),
        line_strengths=np.array([f for _, f in line_list]),
        continuum_energies=energies,
        continuum_density=density,
    )
    constraint = config["constraint"]
    target_data = _catalogue.by_key(gas)[str(constraint["dataset"])]
    row = target_data.rows.index(str(constraint["row"]))
    target = float(target_data.scaled_values()[row])
    order_k = float(constraint.get("k", constraint["row"]))
    constrained_distribution, factor = constrained(
        distribution,
        float(constraint["region_max_ev"]) / EV_PER_HARTREE,
        (order_k, target),
    )
    return ComplexModel(
        gas=gas,
        name=f"euv/{gas.value}",
        distribution=constrained_distribution,
        chi0_au=float(config.get("chi0_au", 0.0)),
        factor=factor,
    )

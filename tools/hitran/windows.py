# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The infrared windows of the molecular gases, their band rules and bound (OEP 3).

    uv run python tools/hitran/windows.py CO CO2 CH4 N2O

reads the cached HITRAN line lists of each gas (``~/.cache/optmat/hitran``, fetched by
``band_strengths.py``) and writes, to ``src/optmat/data/external/hitran/``:

- ``<gas>_band_rules.csv``: every band system (a band split at the windows, the
  rotational lines, and the weak bands inside each window folded onto two poles), its
  totals and the nodes and strengths of its Gaussian rule;
- ``<gas>_windows.csv``: the qualified intervals from 12.5 um to the fit's shortest
  wavelength and the gaps between them, with each gap's reason and the criteria's
  values;
- ``<gas>_bound.csv``: the representation bound B_0, B_2, B_3 per cell of the windows;
- ``<gas>_absorption.csv``: the smoothed absorption a_sigma of the pure gas at the
  reference state, its maximum per 1 cm^-1 cell where it exceeds 1e-4 m^-1 (for the
  windows of a mixture, OEP 3 Q13).

Each file starts with ``#`` lines of provenance: the SHA-256 of the raw inputs, the
query record, every setting, the parent fit whose electronic terms and statistical
sigma set the tolerance, and the tool's commit.

The method (OEP 3, sections 4 and 7):

- **The reference** is the gas's refractivity smoothed by a Gaussian of standard
  deviation sigma in frequency, over every line: each line's resonant term
  1/(E - omega) becomes (2^1/2/sigma) F((E - omega)/(2^1/2 sigma)), F Dawson's
  function. It is evaluated with its first three frequency derivatives by FFT
  convolution on a uniform wavenumber grid, the lines deposited linearly on it, and by
  direct sums above the grid. D_0 = alpha, D_2 = d^2(omega alpha)/d omega^2 and D_3.
- **The windows** start where the smoothed absorption of the pure gas at 101325 Pa and
  296 K is below the threshold; a band's lines are split at them.
- **The model** is the parent's electronic terms plus, per band system, the smallest
  Gaussian rule (at most ``node_cap`` nodes) whose error is below a tenth of the
  tolerance throughout the windows; the in-window bands are kept (a band system whose
  span splits the window) or folded onto the nearest nodes below and above, chosen per
  window for the widest qualified span.
- **The tolerance**: |D_0 - D_0 ref| below ``alpha_fraction`` of the parent's
  statistical sigma of alpha, and |D_2 - D_2 ref| below ``gvd_fraction`` of the larger
  of the reference's and the electronic |D_2|. A window keeps the runs where both hold,
  at least ``min_width`` wide.
- **The bound** is the largest |D_m - D_m ref| in each cell, from a grid of spacing at
  most sigma/4 and a tenth of the distance to the nearest pole, checked on a second
  grid offset by half a step; a cell whose check exceeds it is resampled on the fine
  grid and its bound raised.
"""

import csv
import dataclasses
import datetime
import hashlib
import itertools
import json
import math
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.signal import fftconvolve
from scipy.special import dawsn, wofz

from optmat.gas.infrared import fold_poles, gauss_rule

ROOT = Path(__file__).resolve().parents[2]
CACHE = Path("~/.cache/optmat/hitran").expanduser()
OUT = ROOT / "src" / "optmat" / "data" / "external" / "hitran"
VERSION = 1
GASES = ("CO", "CO2", "CH4", "N2O")
HCM = constants.physical_constants["hartree-inverse meter relationship"][0] / 100
RE_CM = constants.physical_constants["classical electron radius"][0] * 100
K_B = constants.k


@dataclass(frozen=True)
class Settings:
    """Every setting of the derivation (OEP 3, section 7)."""

    sigma_cm: float = 1.0  # the Gaussian's standard deviation, cm^-1
    pressure_pa: float = 101325.0  # the reference state of the absorption
    temperature_k: float = 296.0
    absorption_threshold_m: float = 0.1  # transparency: a_sigma below this, m^-1
    alpha_fraction: float = 0.1  # of the statistical sigma of alpha
    gvd_fraction: float = 0.01  # of max(|D_2 ref|, |D_2 electronic|)
    node_cap: int = 8
    share: float = 0.1  # a system's share of the tolerance
    min_width_cm: float = 50.0
    # Kept in-window bands: those above this fraction of the molecule's M_0; None folds
    # every one. The widest qualified span decides per window.
    keep_fractions: tuple[float | None, ...] = (None, 1e-5, 1e-6, 1e-7, 1e-8)
    low_cm: float = 800.0  # 12.5 um, the long end of the analysis
    fine_cm: float = 0.005  # the FFT grid
    top_cm: float = 11100.0  # the FFT grid's top; direct sums above
    grid_cm: float = 0.25  # the bound's grid, at most sigma / 3
    decision_cm: float = 2.0  # the decision grid of the adaptive rules
    direct_ratio: float = 1.002  # the geometric grid above the FFT grid
    cell_ratio: float = 1.001  # the stored cells of the bound
    wing_cutoff_cm: float = 25.0  # the far wings: beyond this from a line
    near_cm: float = 30.0  # lines nearer than this are summed directly
    absorption_floor_m: float = 1e-4  # absorption cells stored above this


@dataclass
class Lines:
    """A gas's lines: wavenumber, intensity, self-broadened width, band index."""

    nu: np.ndarray  # cm^-1
    strength: np.ndarray  # cm/molecule at 296 K
    gamma_self: np.ndarray  # cm^-1/atm, HWHM
    band: np.ndarray  # index into keys
    keys: list[tuple[str, str, str]]  # (global isotopologue id, upper, lower)

    @property
    def rotational(self) -> np.ndarray:
        """Whether each band is pure rotational (the same state above and below)."""
        return np.array([bool(lo) and up == lo for _, up, lo in self.keys])


def load_lines(gas: str) -> Lines:
    """Read a gas's cached line list."""
    path = CACHE / f"{gas}_0_11000.par"
    index: dict[tuple[str, str, str], int] = {}
    nu, strength, gamma, band = [], [], [], []
    with path.open() as handle:
        for line in handle:
            record, *extra = line.rstrip("\r\n").split(",")
            if len(record) < 160:
                continue
            iso = extra[-1] if extra else record[2]
            key = (iso, record[67:82].strip(), record[82:97].strip())
            band.append(index.setdefault(key, len(index)))
            nu.append(float(record[3:15]))
            strength.append(float(record[15:25]))
            gamma.append(float(record[40:45]))
    order = np.argsort(nu, kind="stable")
    return Lines(
        np.array(nu)[order],
        np.array(strength)[order],
        np.array(gamma)[order],
        np.array(band)[order],
        list(index),
    )


def oscillator(lines_nu: np.ndarray, strength: np.ndarray) -> tuple:
    """Return the lines' oscillator strengths and energies (hartree)."""
    return strength / (math.pi * RE_CM), lines_nu / HCM


# ---------------------------------------------------------------- the reference


def _deposit(nu: np.ndarray, weight: np.ndarray, h: float, size: int) -> np.ndarray:
    """Deposit weights at wavenumbers onto a uniform grid, linearly (two moments)."""
    x = nu / h
    i = np.floor(x).astype(int)
    t = x - i
    grid = np.zeros(size)
    np.add.at(grid, i, weight * (1 - t))
    np.add.at(grid, i + 1, weight * t)
    return grid


def _resonant_kernel(delta: np.ndarray, sigma: float, order: int) -> np.ndarray:
    """Return the order-th derivative in delta of the smoothed 1/delta (hartree)."""
    c = 1 / (math.sqrt(2) * sigma)
    y = delta * c
    # The derivative recurrence loses digits to cancellation at large |y|; there the
    # asymptotic series R = sum_k (2k - 1)!! sigma^2k / delta^(2k+1) (six terms)
    # takes over at |y| = 15, where the two agree to 1e-10 in every order.
    far = np.abs(y) > 15
    with np.errstate(divide="ignore", invalid="ignore"):
        f = [dawsn(np.where(far, 0.0, y))]
        f.append(1 - 2 * y * f[0])
        f.append(-2 * f[0] - 2 * y * f[1])
        f.append(-4 * f[1] - 2 * y * f[2])
        near = math.sqrt(2) / sigma * c**order * f[order]
        d = np.where(far, delta, 1.0)
        series = sum(
            double_factorial
            * sigma ** (2 * k)
            * (-1) ** order
            * math.prod(range(2 * k + 1, 2 * k + 1 + order))
            * d ** -(2 * k + 1 + order)
            for k, double_factorial in enumerate((1, 1, 3, 15, 105, 945))
        )
    return np.where(far, series, near)


def far_derivatives(
    lines_nu: np.ndarray,
    strength: np.ndarray,
    settings: Settings,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the grid and alpha^(j), j = 0..3, of the lines farther than near_cm (au).

    alpha = sum_l f_l / (2 E_l) [R_sigma(E_l - omega) + 1 / (E_l + omega)]. The
    resonant term of the lines beyond ``near_cm`` and the antiresonant term of every
    line are FFT convolutions of the deposited amplitudes with smooth kernels, so the
    roundoff stays below the far values; :func:`near_derivatives` adds the near lines
    exactly. The antiresonant kernel is zero below ``low_cm``, where no output is
    used, so that the lines near zero frequency add no roundoff.
    """
    h, top = settings.fine_cm, settings.top_cm
    f, e = oscillator(lines_nu, strength)
    size = math.ceil(top / h) + 2
    amplitude = _deposit(lines_nu, f / (2 * e), h, size)
    grid = np.arange(size) * h
    sigma = settings.sigma_cm / HCM
    lags = (size - 1 - np.arange(2 * size - 1)) * h  # cm^-1: R((K-1-n) h)
    sums = np.arange(2 * size - 1) * h  # cm^-1: E + omega at index m = i + k
    far = np.abs(lags) > settings.near_cm
    used = sums >= settings.low_cm
    out = np.empty((4, size))
    for j in range(4):
        kernel = np.where(far, (-1) ** j * _resonant_kernel(lags / HCM, sigma, j), 0.0)
        resonant = fftconvolve(amplitude, kernel)[size - 1 : 2 * size - 1]
        with np.errstate(divide="ignore"):
            anti_kernel = np.where(
                used, (-1) ** j * math.factorial(j) / (sums / HCM) ** (j + 1), 0.0
            )
        anti = fftconvolve(amplitude[::-1], anti_kernel)[size - 1 : 2 * size - 1]
        out[j] = resonant + anti
    return grid, out


def near_derivatives(
    lines_nu: np.ndarray, strength: np.ndarray, settings: Settings, nu: np.ndarray
) -> np.ndarray:
    """Return alpha^(j) of the resonant terms of the lines within near_cm (au).

    ``lines_nu`` must be sorted.
    """
    f, e = oscillator(lines_nu, strength)
    sigma = settings.sigma_cm / HCM
    low = np.searchsorted(lines_nu, nu - settings.near_cm, side="left")
    high = np.searchsorted(lines_nu, nu + settings.near_cm, side="right")
    out = np.zeros((4, len(nu)))
    # Points in blocks, each with the lines any of its points needs.
    for start in range(0, len(nu), 256):
        stop = min(start + 256, len(nu))
        a, b = low[start:stop].min(), high[start:stop].max()
        if b <= a:
            continue
        ee, ff = e[a:b][None, :], f[a:b][None, :]
        delta_cm = lines_nu[a:b][None, :] - nu[start:stop, None]
        near = np.abs(delta_cm) <= settings.near_cm
        w = nu[start:stop, None] / HCM
        for j in range(4):
            term = (-1) ** j * _resonant_kernel(ee - w, sigma, j)
            out[j, start:stop] = np.sum(
                np.where(near, ff / (2 * ee) * term, 0.0), axis=1
            )
    return out


def direct_derivatives(
    lines_nu: np.ndarray, strength: np.ndarray, sigma_cm: float, nu: np.ndarray
) -> np.ndarray:
    """Return alpha^(j), j = 0..3, of the smoothed line sum by direct summation (au)."""
    f, e = oscillator(lines_nu, strength)
    omega = nu / HCM
    sigma = sigma_cm / HCM
    out = np.zeros((4, len(nu)))
    for start in range(0, len(f), 4000):
        ff = f[start : start + 4000][None, :]
        ee = e[start : start + 4000][None, :]
        w = omega[:, None]
        for j in range(4):
            resonant = (-1) ** j * _resonant_kernel(ee - w, sigma, j)
            anti = (-1) ** j * math.factorial(j) / (ee + w) ** (j + 1)
            out[j] += np.sum(ff / (2 * ee) * (resonant + anti), axis=1)
    return out


def dispersion(derivatives: np.ndarray, omega: np.ndarray) -> dict[int, np.ndarray]:
    """Return D_0 = alpha, D_2 and D_3 of d^m(omega alpha)/d omega^m."""
    a0, a1, a2, a3 = derivatives
    return {0: a0, 2: omega * a2 + 2 * a1, 3: omega * a3 + 3 * a2}


def pole_dispersion(
    f: np.ndarray, e: np.ndarray, omega: np.ndarray
) -> dict[int, np.ndarray]:
    """Return D_0, D_2, D_3 of undamped poles sum f / (E^2 - omega^2) (au)."""
    out = {m: np.zeros(len(omega)) for m in (0, 2, 3)}
    for start in range(0, len(f), 2000):
        ff = f[start : start + 2000][None, :]
        ee = e[start : start + 2000][None, :]
        w = omega[:, None]
        out[0] += np.sum(ff / (ee**2 - w**2), axis=1)
        for m in (2, 3):
            fac = 0.5 * math.factorial(m)
            t = fac * ((ee - w) ** -(m + 1) - (-1) ** m * (ee + w) ** -(m + 1))
            out[m] += np.sum(ff * t, axis=1)
    return out


def smoothed_absorption(
    lines_nu: np.ndarray, strength: np.ndarray, settings: Settings, size: int
) -> np.ndarray:
    """Return a_sigma = N sum S g_sigma on the fine grid, m^-1."""
    h = settings.fine_cm
    deposited = _deposit(lines_nu, strength, h, size)
    x = np.arange(-10 * settings.sigma_cm, 10 * settings.sigma_cm + h / 2, h)
    gauss = np.exp(-(x**2) / (2 * settings.sigma_cm**2)) / (
        settings.sigma_cm * math.sqrt(2 * math.pi)
    )
    smoothed = fftconvolve(deposited, gauss, mode="same")
    # The FFT leaves roundoff of either sign where there are no lines.
    return np.maximum(smoothed, 0.0) * _number_density(settings) * 100


def far_wing_absorption(lines: Lines, settings: Settings, size: int) -> np.ndarray:
    """Return the Lorentz far wings, beyond the cutoff of every line, m^-1, at 1 atm.

    Real wings fall below the Lorentzian far from the line (the chi factor); this is
    an estimate from above of the self-broadened wings at the reference pressure.
    """
    h = settings.fine_cm
    atm = settings.pressure_pa / 101325.0
    deposited = _deposit(
        lines.nu, lines.strength * lines.gamma_self * atm / math.pi, h, size
    )
    x = np.arange(-2000.0, 2000.0 + h / 2, h)
    kernel = np.where(
        np.abs(x) > settings.wing_cutoff_cm, 1 / np.maximum(x**2, 1e-30), 0
    )
    wings = fftconvolve(deposited, kernel, mode="same")
    return np.maximum(wings, 0.0) * _number_density(settings) * 100


def _number_density(settings: Settings) -> float:
    """Return the reference state's number density, cm^-3."""
    return settings.pressure_pa / (K_B * settings.temperature_k) * 1e-6


# ------------------------------------------------------------------ the windows


def runs(mask: np.ndarray, nu: np.ndarray) -> list[tuple[float, float]]:
    """Return the (low, high) wavenumbers of each run of True in a mask on a grid."""
    edges = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    starts, stops = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1) - 1
    return [(float(nu[a]), float(nu[b])) for a, b in zip(starts, stops, strict=True)]


@dataclass
class System:
    """A band system: a group of lines represented by one rule, or a fold."""

    treatment: str  # "rule", "rotational" or "folded"
    key: tuple[str, str, str]
    where: str  # the segment the lines lie in
    index: np.ndarray  # line indices
    strengths: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros(0))
    energies: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros(0))
    bands: int = 1


@dataclass
class Parent:
    """The parent fit's electronic terms and statistical sigma."""

    name: str
    config_sha256: str
    strengths: np.ndarray
    energies: np.ndarray
    offset: float
    upper_cm: float  # the fit's shortest wavelength, as a wavenumber
    optical_cm: float  # the fit's longest wavelength (range_nm), as a wavenumber
    model: object

    def electronic(self, omega: np.ndarray) -> dict[int, np.ndarray]:
        """Return D_0, D_2, D_3 of the electronic terms (and the offset)."""
        out = pole_dispersion(self.strengths, self.energies, omega)
        out[0] = out[0] + self.offset
        return out

    def sigma(self, omega: np.ndarray) -> np.ndarray:
        """Return the fit's statistical sigma of alpha (au)."""
        return np.asarray(self.model.polarizability_sigma_au(omega))  # pyright: ignore[reportAttributeAccessIssue]


def load_parent(gas: str) -> Parent:
    """Return the stored fit of a gas."""
    from optmat.gas.models.registry import model  # noqa: PLC0415

    fit = model(f"fit/{gas}")
    n_fit = len(fit.strengths) - fit.n_fixed  # pyright: ignore[reportAttributeAccessIssue]
    return Parent(
        name=fit.name,
        config_sha256=str(fit.record["config_sha256"]),  # pyright: ignore[reportAttributeAccessIssue]
        strengths=np.asarray(fit.strengths[:n_fit]),  # pyright: ignore[reportAttributeAccessIssue]
        energies=np.asarray(fit.energies[:n_fit]),  # pyright: ignore[reportAttributeAccessIssue]
        offset=float(fit.offset_au),  # pyright: ignore[reportAttributeAccessIssue]
        upper_cm=1e-2 / fit.wavelength_range[0],
        optical_cm=1e-2 / fit.wavelength_range[1],
        model=fit,
    )


@dataclass
class Reference:
    """The reference and the tolerance, on the fine grid and above it.

    The fine arrays hold every point of the FFT grid from ``low_cm``; the analysis grid
    is every ``grid_cm`` of them, followed by the direct points above the grid.
    """

    fine_nu: np.ndarray
    fine_d: np.ndarray  # alpha^(j) of the far lines, au, j = 0..3
    fine_absorption: np.ndarray  # a_sigma, m^-1
    fine_wing: np.ndarray  # far-wing absorption, m^-1
    direct_nu: np.ndarray
    direct_d: dict[int, np.ndarray]
    step: int  # fine points per grid step
    parent: Parent
    settings: Settings
    lines: Lines

    def sample(self, nu: np.ndarray) -> dict[str, np.ndarray]:
        """Return the reference at points of the fine grid or of the direct grid.

        Returns D_0, D_2, D_3 (electronic plus every line), the absorption, the far-wing
        absorption and the two tolerances.
        """
        h = self.settings.fine_cm
        on_fine = nu <= self.fine_nu[-1] + h / 2
        index = np.clip(
            np.rint((nu - self.fine_nu[0]) / h).astype(int), 0, len(self.fine_nu) - 1
        )
        omega = nu / HCM
        fine_points = nu[on_fine]
        alpha = self.fine_d[:, index[on_fine]] + near_derivatives(
            self.lines.nu, self.lines.strength, self.settings, fine_points
        )
        lines_d = dispersion(alpha, fine_points / HCM)
        electronic = self.parent.electronic(omega)
        where = np.searchsorted(self.direct_nu, nu[~on_fine]).clip(
            0, max(len(self.direct_nu) - 1, 0)
        )
        out: dict[str, np.ndarray] = {}
        for m in (0, 2, 3):
            values = np.empty(len(nu))
            values[on_fine] = lines_d[m]
            values[~on_fine] = self.direct_d[m][where] if (~on_fine).any() else 0.0
            out[f"d{m}"] = electronic[m] + values
        out["absorption"] = np.where(on_fine, self.fine_absorption[index], 0.0)
        out["wing"] = np.where(on_fine, self.fine_wing[index], 0.0)
        out["tol0"] = self.settings.alpha_fraction * self.parent.sigma(omega)
        out["tol2"] = self.settings.gvd_fraction * np.maximum(
            np.abs(out["d2"]), np.abs(electronic[2])
        )
        return out

    @property
    def grid(self) -> np.ndarray:
        """The analysis grid: every grid step of the fine grid, then the direct grid."""
        return np.concatenate([self.fine_nu[:: self.step], self.direct_nu])


def build_reference(lines: Lines, parent: Parent, settings: Settings) -> Reference:
    """Return the reference from 12.5 um to the fit's shortest wavelength."""
    grid, alpha = far_derivatives(lines.nu, lines.strength, settings)
    size = len(grid)
    absorption = smoothed_absorption(lines.nu, lines.strength, settings, size)
    wing = far_wing_absorption(lines, settings, size)
    keep = (grid >= settings.low_cm - settings.fine_cm / 2) & (
        grid <= settings.top_cm - 50
    )
    nu_fine = grid[keep]
    fine_d = alpha[:, keep]
    count = int(
        math.log(max(parent.upper_cm, nu_fine[-1]) / nu_fine[-1])
        / math.log(settings.direct_ratio)
    )
    nu_direct = nu_fine[-1] * settings.direct_ratio ** np.arange(1, count + 1)
    nu_direct = nu_direct[nu_direct <= parent.upper_cm]
    omega_direct = nu_direct / HCM
    direct_lines = dispersion(
        direct_derivatives(lines.nu, lines.strength, settings.sigma_cm, nu_direct),
        omega_direct,
    )

    return Reference(
        fine_nu=nu_fine,
        fine_d=fine_d,
        fine_absorption=absorption[keep],
        fine_wing=wing[keep],
        direct_nu=nu_direct,
        direct_d=direct_lines,
        step=round(settings.grid_cm / settings.fine_cm),
        parent=parent,
        settings=settings,
        lines=lines,
    )


def _segments(windows: list[tuple[float, float]]) -> list[tuple[str, float, float]]:
    """Return the windows and the gaps around them, as (name, low, high)."""
    out, low = [], 0.0
    for k, (a, b) in enumerate(windows):
        out += [(f"gap {k}", low, a), (f"window {k}", a, b)]
        low = b
    return [*out, (f"gap {len(windows)}", low, math.inf)]


def _pieces(
    window: tuple[float, float], excluded: list[tuple[float, float]]
) -> list[tuple[float, float]]:
    """Return a window minus excluded spans."""
    pieces = [window]
    for a, b in sorted(excluded):
        pieces = [
            part
            for low, high in pieces
            for part in ((low, min(high, a)), (max(low, b), high))
            if part[1] > part[0]
        ]
    return pieces


def partition(
    lines: Lines, windows: list[tuple[float, float]]
) -> tuple[list[System], dict[int, list[System]]]:
    """Split every band at the windows: systems in the gaps, parts in the windows."""
    rotational = lines.rotational
    segments = _segments(windows)
    edges = np.array([low for _, low, _ in segments[1:]])
    where = np.searchsorted(edges, lines.nu, side="right")
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i, (b, s) in enumerate(zip(lines.band, where, strict=True)):
        groups[(int(b), int(s))].append(i)
    rot_index = np.flatnonzero(rotational[lines.band])
    systems = (
        [System("rotational", ("*", "rotational", ""), "gap 0", rot_index)]
        if len(rot_index)
        else []
    )
    inside: dict[int, list[System]] = defaultdict(list)
    for (b, s), index in groups.items():
        if rotational[b]:
            continue
        name = segments[s][0]
        system = System("rule", lines.keys[b], name, np.array(index))
        if name.startswith("window"):
            inside[int(name.split()[1])].append(system)
        else:
            systems.append(system)
    return systems, inside


_RULES: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}


def adaptive_rule(
    system: System,
    f: np.ndarray,
    e: np.ndarray,
    lines_nu: np.ndarray,
    decision: dict[str, np.ndarray],
    settings: Settings,
) -> None:
    """Give a system the smallest rule within its share of the tolerance.

    The share is checked at the decision points of the candidate windows outside the
    system's own span; the result is cached by the system's lines.
    """
    key = (system.key, system.index.tobytes())
    if key not in _RULES:
        fi, ei = f[system.index], e[system.index]
        span = (lines_nu[system.index].min(), lines_nu[system.index].max())
        away = (decision["nu"] < span[0] - 1.0) | (decision["nu"] > span[1] + 1.0)
        omega = decision["nu"][away] / HCM
        tol0, tol2 = decision["tol0"][away], decision["tol2"][away]
        exact = pole_dispersion(fi, ei, omega)
        strengths, energies = gauss_rule(fi, ei, 1)
        for n in range(1, settings.node_cap + 1):
            strengths, energies = gauss_rule(fi, ei, n)
            rule = pole_dispersion(strengths, energies, omega)
            good = np.all(
                np.abs(rule[0] - exact[0]) <= settings.share * tol0
            ) and np.all(np.abs(rule[2] - exact[2]) <= settings.share * tol2)
            if good or n == settings.node_cap or len(energies) < n:
                break
        _RULES[key] = (strengths, energies)
    system.strengths, system.energies = _RULES[key]


@dataclass
class Derivation:
    """The result of one derivation: the systems, windows and bounds."""

    systems: list[System]
    windows: list[tuple[float, float]]  # the qualified intervals, cm^-1
    candidates: list[tuple[float, float]]  # the transparent candidates
    pieces: list[tuple[float, float]]  # the candidates minus the kept bands' spans
    choice: dict[int, float | None]  # the kept fraction per candidate
    kept: dict[int, list[System]]
    folded: dict[int, int]  # per candidate, the number of folded band parts
    poles: tuple[np.ndarray, np.ndarray]  # every infrared pole: strength, energy


@dataclass
class Analysis:
    """A gas's lines, parent, reference and its values on the analysis grid."""

    gas: str
    lines: Lines
    parent: Parent
    reference: Reference
    settings: Settings
    grid: np.ndarray
    at: dict[str, np.ndarray]  # the reference sampled on the grid
    candidates: list[tuple[float, float]]
    decision: dict[str, np.ndarray]


def analyse(gas: str, settings: Settings) -> Analysis:
    """Load a gas and compute its reference and candidate windows."""
    lines = load_lines(gas)
    parent = load_parent(gas)
    reference = build_reference(lines, parent, settings)
    grid = reference.grid
    at = reference.sample(grid)
    transparent = at["absorption"] < settings.absorption_threshold_m
    candidates = [
        (a, b) for a, b in runs(transparent, grid) if b - a >= settings.min_width_cm
    ]
    inside = np.zeros(len(grid), bool)
    for a, b in candidates:
        inside |= (grid >= a) & (grid <= b)
    on_step = np.isclose(
        np.round(grid / settings.decision_cm) * settings.decision_cm, grid
    )
    direct = grid > reference.fine_nu[-1]
    pick = inside & (on_step | _near_edges(grid, candidates, 20.0) | direct)
    decision = {"nu": grid[pick], "tol0": at["tol0"][pick], "tol2": at["tol2"][pick]}
    return Analysis(
        gas, lines, parent, reference, settings, grid, at, candidates, decision
    )


def derive(analysis: Analysis, choice: dict[int, float | None]) -> Derivation:
    """Derive the systems and pieces for a choice of kept bands per candidate window."""
    lines, settings = analysis.lines, analysis.settings
    f, e = oscillator(lines.nu, lines.strength)
    _, inside = partition(lines, analysis.candidates)
    vibrational = ~lines.rotational[lines.band]
    total_m0 = float(f[vibrational].sum())
    kept_all: dict[int, list[System]] = {}
    pieces: list[tuple[float, float]] = []
    folded_count: dict[int, int] = {}
    for k, window in enumerate(analysis.candidates):
        fraction = choice.get(k)
        parts = inside.get(k, [])
        # A band inside the optical range is always folded: no pole may lie in a
        # qualified interval, and there the tolerance does not split the window.
        kept = [
            p
            for p in parts
            if fraction is not None
            and f[p.index].sum() / total_m0 > fraction
            and lines.nu[p.index].max() + settings.sigma_cm < analysis.parent.optical_cm
        ]
        kept_all[k] = kept
        spans = [
            (
                float(lines.nu[p.index].min()) - settings.sigma_cm,
                float(lines.nu[p.index].max()) + settings.sigma_cm,
            )
            for p in kept
        ]
        pieces += [
            piece
            for piece in _pieces(window, spans)
            if piece[1] - piece[0] >= settings.min_width_cm
        ]
        folded_count[k] = len(parts) - len(kept)
    # Every band is split at the pieces: the lines in a piece are folded there, the rest
    # are rule systems in the gaps, the kept bands among them.
    systems, inside = partition(lines, pieces)
    for system in systems:
        adaptive_rule(system, f, e, lines.nu, analysis.decision, settings)
    nodes = np.concatenate([s.energies for s in systems])
    folds = []
    for k, (a, b) in enumerate(pieces):
        parts = inside.get(k, [])
        if not parts:
            continue
        index = np.concatenate([p.index for p in parts])
        below = nodes[nodes < a / HCM].max()
        # Above the highest window no band node is left: the parent's lowest
        # electronic energy, beyond every interval, takes the fold's upper pole.
        higher = nodes[nodes > b / HCM]
        above = higher.min() if len(higher) else float(analysis.parent.energies.min())
        low, high = fold_poles(
            float(f[index].sum()), float(np.sum(f[index] / e[index] ** 2)), below, above
        )
        folds.append(
            System(
                "folded",
                ("*", "folded", ""),
                f"window {k}",
                index,
                np.array([low, high]),
                np.array([below, above]),
                bands=len(parts),
            )
        )
    everything = systems + folds
    poles = (
        np.concatenate([s.strengths for s in everything]),
        np.concatenate([s.energies for s in everything]),
    )
    windows = _qualify(analysis, pieces, poles, analysis.grid, analysis.at)
    return Derivation(
        everything,
        windows,
        analysis.candidates,
        pieces,
        choice,
        kept_all,
        folded_count,
        poles,
    )


def model_error(
    analysis: Analysis, poles: tuple[np.ndarray, np.ndarray], nu: np.ndarray, at: dict
) -> dict[int, np.ndarray]:
    """Return |D_m[model] - D_m[reference]| at points."""
    omega = nu / HCM
    electronic = analysis.parent.electronic(omega)
    infrared = pole_dispersion(*poles, omega)
    return {m: np.abs(electronic[m] + infrared[m] - at[f"d{m}"]) for m in (0, 2, 3)}


def _qualify(
    analysis: Analysis,
    pieces: list[tuple[float, float]],
    poles: tuple[np.ndarray, np.ndarray],
    nu: np.ndarray,
    at: dict,
) -> list[tuple[float, float]]:
    """Return the runs of the pieces where the tolerance holds, wide enough."""
    settings = analysis.settings
    in_piece = np.zeros(len(nu), bool)
    for a, b in pieces:
        in_piece |= (nu >= a) & (nu <= b)
    error = model_error(analysis, poles, nu, at)
    represented = (error[0] <= at["tol0"]) & (error[2] <= at["tol2"])
    ok = (
        in_piece
        & (at["absorption"] < settings.absorption_threshold_m)
        & (represented | ~enforced(analysis, nu))
    )
    return [(a, b) for a, b in runs(ok, nu) if b - a >= settings.min_width_cm]


def enforced(analysis: Analysis, nu: np.ndarray) -> np.ndarray:
    """Return where the representation tolerance qualifies a window.

    Beyond the parent fit's optical range, at longer wavelengths. Inside that range the
    fit's data qualify it; a window there needs only transparency, and the bound is
    stated but not enforced (the user's decision of 2026-09-30, recorded with OEP 3's
    Phase 1: at sigma = 1 cm^-1 the lines of weak near-infrared bands give the GVD a
    structure no pole model follows).
    """
    return nu < analysis.parent.optical_cm


def _near_edges(nu: np.ndarray, pieces: list, width: float) -> np.ndarray:
    """Return where the grid lies within a width of a piece's edge."""
    near = np.zeros(len(nu), bool)
    for a, b in pieces:
        near |= (np.abs(nu - a) <= width) | (np.abs(nu - b) <= width)
    return near


def choose(analysis: Analysis) -> Derivation:
    """Derive with the kept fraction that qualifies the widest span, per window."""
    settings = analysis.settings
    trials = {}
    for fraction in settings.keep_fractions:
        result = derive(
            analysis, dict.fromkeys(range(len(analysis.candidates)), fraction)
        )
        trials[fraction] = result
        width = sum(b - a for a, b in result.windows)
        print(
            f"  keep > {fraction}: {width:.0f} cm^-1 qualified in "
            f"{len(result.windows)} intervals, {len(result.poles[0])} poles"
        )
    choice = {}
    for k, (a, b) in enumerate(analysis.candidates):

        def overlap(result: Derivation, a: float = a, b: float = b) -> float:
            return sum(max(0.0, min(y, b) - max(x, a)) for x, y in result.windows)

        choice[k] = max(settings.keep_fractions, key=lambda fr: overlap(trials[fr]))
    return derive(analysis, choice)


# -------------------------------------------------------------- the final bound


@dataclass
class Final:
    """The qualified intervals, the gaps and the checked bound of one gas."""

    derivation: Derivation
    windows: list[tuple[float, float]]
    gaps: list[dict]
    cells: list[dict]
    checks: dict[str, int]
    peaks: np.ndarray  # each line's self-broadened Lorentz peak at 1 atm, m^-1


def _nearest_pole(poles_e: np.ndarray, nu: np.ndarray) -> np.ndarray:
    """Return the distance from each point to the nearest pole, cm^-1."""
    e_cm = np.sort(poles_e * HCM)
    i = np.searchsorted(e_cm, nu).clip(1, len(e_cm) - 1)
    return np.minimum(np.abs(nu - e_cm[i - 1]), np.abs(e_cm[i] - nu))


def finalize(analysis: Analysis, derivation: Derivation) -> Final:
    """Check the bound between the grid's points, and describe every gap.

    Each window's grid cells (spacing ``grid_cm``) are checked on a grid offset by half
    a step; where a pole is nearer than ten grid steps the cell is sampled at a tenth
    of that distance, down to the fine grid, and a cell whose pole is nearer than ten
    fine steps is dropped. A check sample above its cell's bound, or above the
    tolerance, resamples the cell on the fine grid and raises its bound or drops it.
    The windows are then the runs of qualified cells.
    """
    settings = analysis.settings
    grid, at = analysis.grid, analysis.at
    poles = derivation.poles
    error = model_error(analysis, poles, grid, at)
    ok = np.zeros(len(grid), bool)
    for a, b in derivation.windows:
        ok |= (grid >= a) & (grid <= b)
    cell_bound, checks = _check_cells(analysis, poles, error, ok)
    windows = [(a, b) for a, b in runs(ok, grid) if b - a >= settings.min_width_cm]
    cells = _store_cells(analysis, windows, cell_bound)
    gaps = _gaps(analysis, derivation, windows, error)
    lines = analysis.lines
    peaks = (
        _number_density(settings)
        * lines.strength
        / (math.pi * np.maximum(lines.gamma_self, 1e-6))
        * 100
    )
    return Final(derivation, windows, gaps, cells, checks, peaks)


def _check_cells(
    analysis: Analysis,
    poles: tuple[np.ndarray, np.ndarray],
    error: dict[int, np.ndarray],
    ok: np.ndarray,
) -> tuple[dict[int, np.ndarray], dict[str, int]]:
    """Check every qualified cell between its grid points; update ``ok`` in place.

    Returns the bound per grid point (the larger of a cell's two ends and its checks)
    and the counts of the checks.
    """
    settings, reference, grid = analysis.settings, analysis.reference, analysis.grid
    h, step = settings.fine_cm, settings.grid_cm
    fine_top = reference.fine_nu[-1]
    distance = _nearest_pole(poles[1], grid)
    cell_bound = {m: error[m].copy() for m in (0, 2, 3)}
    cells = np.flatnonzero(ok[:-1] & ok[1:])
    low, high = grid[cells], grid[cells + 1]
    nearest = np.minimum(distance[cells], distance[cells + 1])
    direct = high > fine_top
    dropped = ~direct & (0.1 * nearest < h)
    spacing = np.where(
        direct, high - low, np.minimum(step / 2, np.maximum(h, 0.1 * nearest))
    )
    counts = np.where(direct, 2, np.maximum(2, np.rint((high - low) / spacing))).astype(
        int
    )
    counts[dropped] = 0
    owner = np.repeat(np.arange(len(cells)), np.maximum(counts - 1, 0))
    offset = np.concatenate(
        [np.arange(1, c) / c for c in counts if c > 1] or [np.zeros(0)]
    )
    samples = low[owner] + offset * (high - low)[owner]
    on_fine = ~direct[owner]
    samples[on_fine] = np.round(samples[on_fine] / h) * h
    sampled = {k: np.empty(len(samples)) for k in ("d0", "d2", "d3", "tol0", "tol2")}
    for mask, source in ((on_fine, reference.sample), (~on_fine, None)):
        if not mask.any():
            continue
        values = (
            source(samples[mask]) if source else _direct_sample(analysis, samples[mask])
        )
        for k, target in sampled.items():
            target[mask] = values[k]
    e = model_error(analysis, poles, samples, sampled)
    failed = np.zeros(len(cells), bool)
    outside = (e[0] > sampled["tol0"]) | (e[2] > sampled["tol2"])
    np.logical_or.at(failed, owner, outside & enforced(analysis, samples))
    failed |= dropped
    worst = {m: np.zeros(len(cells)) for m in (0, 2, 3)}
    exceeded = np.zeros(len(cells), bool)
    for m in (0, 2, 3):
        np.maximum.at(worst[m], owner, e[m])
        exceeded |= worst[m] > np.maximum(error[m][cells], error[m][cells + 1])
    redo = np.flatnonzero(exceeded & ~failed)
    _resample(analysis, poles, redo, (low, high, direct), failed, worst)
    for m in (0, 2, 3):
        top = np.maximum(error[m][cells], error[m][cells + 1])
        cell_bound[m][cells] = np.maximum(top, worst[m])
    raised = len(redo)
    ok[cells[failed]] = False
    ok[cells[failed] + 1] = False
    checks = {
        "cells": len(cells),
        "samples": len(samples),
        "refined": int(np.sum(~direct & (spacing < step / 2))),
        "dropped": int(np.sum(failed)),
        "resampled": raised,
    }
    return cell_bound, checks


def _resample(
    analysis: Analysis,
    poles: tuple[np.ndarray, np.ndarray],
    redo: np.ndarray,
    cells: tuple[np.ndarray, np.ndarray, np.ndarray],
    failed: np.ndarray,
    worst: dict[int, np.ndarray],
) -> None:
    """Resample cells whose check exceeded their bound; update failures and bounds."""
    reference, h = analysis.reference, analysis.settings.fine_cm
    low, high, direct = cells
    # A check above its cell's bound: the cell is resampled on the fine grid (on the
    # direct grid, at ten points), and its bound is the largest error found.
    fine_cells = redo[~direct[redo]]
    direct_cells = redo[direct[redo]]
    for chosen, is_fine in ((fine_cells, True), (direct_cells, False)):
        if not len(chosen):
            continue
        per = [
            np.arange(float(low[c]) + h, float(high[c]) - h / 2, h)
            if is_fine
            else low[c] + np.arange(1, 10) / 10 * (high[c] - low[c])
            for c in chosen
        ]
        who = np.repeat(np.arange(len(chosen)), [len(x) for x in per])
        points = np.concatenate(per)
        values = (
            reference.sample(points) if is_fine else _direct_sample(analysis, points)
        )
        found = model_error(analysis, poles, points, values)
        bad = np.zeros(len(chosen), bool)
        outside = (found[0] > values["tol0"]) | (found[2] > values["tol2"])
        np.logical_or.at(bad, who, outside & enforced(analysis, points))
        failed[chosen[bad]] = True
        for m in (0, 2, 3):
            np.maximum.at(worst[m], chosen[who], found[m])


def _direct_sample(analysis: Analysis, nu: np.ndarray) -> dict[str, np.ndarray]:
    """Return the reference at points above the FFT grid, by direct sums."""
    settings, parent = analysis.settings, analysis.parent
    omega = nu / HCM
    lines = dispersion(
        direct_derivatives(
            analysis.lines.nu, analysis.lines.strength, settings.sigma_cm, nu
        ),
        omega,
    )
    electronic = parent.electronic(omega)
    out = {f"d{m}": electronic[m] + lines[m] for m in (0, 2, 3)}
    out["absorption"] = np.zeros(len(nu))
    out["wing"] = np.zeros(len(nu))
    out["tol0"] = settings.alpha_fraction * parent.sigma(omega)
    out["tol2"] = settings.gvd_fraction * np.maximum(
        np.abs(out["d2"]), np.abs(electronic[2])
    )
    return out


def _store_cells(
    analysis: Analysis, windows: list[tuple[float, float]], bound: dict[int, np.ndarray]
) -> list[dict]:
    """Return the bound's stored cells: geometric, the largest bound in each."""
    grid, settings = analysis.grid, analysis.settings
    cells = []
    for a, b in windows:
        edges = [a]
        while edges[-1] < b:
            edges.append(min(b, edges[-1] * settings.cell_ratio))
        for low, high in itertools.pairwise(edges):
            inside = (grid >= low) & (grid <= high)
            if not inside.any():  # a cell narrower than the grid: its two neighbours
                j = np.searchsorted(grid, low)
                inside[max(j - 1, 0) : j + 1] = True
            cells.append(
                {"low": low, "high": high}
                | {f"b{m}": float(bound[m][inside].max()) for m in (0, 2, 3)}
            )
    return cells


def _gaps(
    analysis: Analysis,
    derivation: Derivation,
    windows: list[tuple[float, float]],
    error: dict[int, np.ndarray],
) -> list[dict]:
    """Return the gaps between the windows, with the conditions that fail in each.

    A point fails by absorption (a_sigma at or above the threshold), by representation
    in a kept band's span, by representation elsewhere (the bound above the tolerance,
    where it is enforced), or, where it passes, by width: its run is narrower than the
    minimum. A gap's reason is the condition failing over most of it; ``causes`` gives
    the share of its points each condition fails at.
    """
    settings, grid, at = analysis.settings, analysis.grid, analysis.at
    in_window = np.zeros(len(grid), bool)
    for a, b in windows:
        in_window |= (grid >= a) & (grid <= b)
    kept_span = np.zeros(len(grid), bool)
    for kept in derivation.kept.values():
        for p in kept:
            nu = analysis.lines.nu[p.index]
            kept_span |= (grid >= nu.min() - settings.sigma_cm) & (
                grid <= nu.max() + settings.sigma_cm
            )
    absorbing = at["absorption"] >= settings.absorption_threshold_m
    failing = ((error[0] > at["tol0"]) | (error[2] > at["tol2"])) & enforced(
        analysis, grid
    )
    cause = np.select(
        [absorbing, kept_span, failing],
        ["absorption", "representation: kept band", "representation"],
        default="width",
    )
    edges = [grid[0], *[x for w in windows for x in w], grid[-1]]
    gaps = []
    for low, high in zip(edges[0::2], edges[1::2], strict=True):
        inside = (grid > low) & (grid < high)
        if high <= low or not inside.any():
            continue
        names, counts = np.unique(cause[inside], return_counts=True)
        share = dict(zip(names.tolist(), (counts / counts.sum()).tolist(), strict=True))
        gaps.append(
            {
                "low": float(low),
                "high": float(high),
                "reason": max(share, key=share.__getitem__),
                "causes": "; ".join(
                    f"{name} {share[name]:.2f}"
                    for name in sorted(share, key=share.__getitem__, reverse=True)
                ),
                "max_absorption_m": float(at["absorption"][inside].max()),
                "kept_bands": int(np.any(kept_span[inside])),
            }
        )
    return gaps


# ------------------------------------------------------ the collisional widths


def _voigt_kernel(
    delta: np.ndarray, sigma: float, gamma: float, order: int
) -> np.ndarray:
    """Return the order-th delta derivative of the Voigt-smoothed dispersive 1/delta.

    Im[(pi/2)^1/2 / sigma w(z)], z = (delta + i gamma) / (2^1/2 sigma), w Faddeeva's
    function: the Gaussian average of delta'/(delta'^2 + gamma^2); gamma -> 0 gives
    the Dawson kernel.
    """
    c = 1 / (math.sqrt(2) * sigma)
    z = (delta + 1j * gamma) * c
    w = [wofz(z)]
    w.append(-2 * z * w[0] + 2j / math.sqrt(math.pi))
    w.append(-2 * w[0] - 2 * z * w[1])
    w.append(-4 * w[1] - 2 * z * w[2])
    return np.imag(math.sqrt(math.pi) * c * c**order * w[order])


def collisional(
    analysis: Analysis, final: Final, pressures_atm: tuple[float, ...] = (1, 3, 10)
) -> list[dict]:
    """Return how far collisional widths move the reference, against the tolerance.

    The reference neglects collisional broadening. Here every line within 200 cm^-1 of
    a window point takes the intensity-weighted mean self-broadened HWHM times the
    pressure, and the change of D_0 and D_2 from the Gaussian-smoothed reference is
    compared with the tolerance at each 1 cm^-1 of every window.
    """
    lines, settings = analysis.lines, analysis.settings
    gamma_cm = float(np.sum(lines.strength * lines.gamma_self) / np.sum(lines.strength))
    f, e = oscillator(lines.nu, lines.strength)
    sigma = settings.sigma_cm / HCM
    rows = []
    for a, b in final.windows:
        nu = np.arange(math.ceil(a) + 0.5, min(b, analysis.reference.fine_nu[-1]), 1.0)
        if not len(nu):
            continue
        at = analysis.reference.sample(
            np.round(nu / settings.fine_cm) * settings.fine_cm
        )
        near = (lines.nu > nu.min() - 200) & (lines.nu < nu.max() + 200)
        ee, ff = e[near][None, :], f[near][None, :]
        w = nu[:, None] / HCM
        gauss = [
            np.sum(
                ff / (2 * ee) * (-1) ** j * _resonant_kernel(ee - w, sigma, j), axis=1
            )
            for j in range(4)
        ]
        for p in pressures_atm:
            gamma = gamma_cm * p / HCM
            voigt = [
                np.sum(
                    ff / (2 * ee) * (-1) ** j * _voigt_kernel(ee - w, sigma, gamma, j),
                    axis=1,
                )
                for j in range(4)
            ]
            change = dispersion(np.array(voigt) - np.array(gauss), nu / HCM)
            rows.append(
                {
                    "window": f"{1e4 / b:.3f}-{1e4 / a:.3f} um",
                    "pressure_atm": p,
                    "gamma_cm": gamma_cm * p,
                    "d0_over_tolerance": float(np.max(np.abs(change[0]) / at["tol0"])),
                    "d2_over_tolerance": float(np.max(np.abs(change[2]) / at["tol2"])),
                }
            )
    return rows


# -------------------------------------------------------------------- the files


def _sha256(path: Path) -> str:
    """Return a file's SHA-256."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _commit() -> str:
    """Return the tool's commit, marked when the working tree differs."""
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        [
            "git",
            "status",
            "--porcelain",
            "--",
            "tools/hitran",
            "src/optmat/gas/infrared.py",
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.strip()
    return head + ("+modified" if dirty else "")


def provenance(analysis: Analysis) -> list[str]:
    """Return the ``#`` lines that head every derived file of a gas."""
    gas = analysis.gas
    par = CACHE / f"{gas}_0_11000.par"
    query = CACHE / f"{gas}_0_11000.json"
    parent = analysis.parent
    settings = json.dumps(dataclasses.asdict(analysis.settings), sort_keys=True)
    return [
        f"# Derived by tools/hitran/windows.py (version {VERSION}, commit {_commit()}) "
        f"on {datetime.datetime.now(datetime.UTC).date().isoformat()}; OEP 3.",
        f"# input {par.name}: sha256 {_sha256(par)}",
        f"# input {query.name}: sha256 {_sha256(query)}; "
        f"query {query.read_text().strip()}",
        f"# parent {parent.name}: config_sha256 {parent.config_sha256}",
        f"# settings {settings}",
    ]


def _up(value: float) -> str:
    """Return a bound rounded up to four significant digits, so it stays a bound."""
    if value <= 0 or not math.isfinite(value):
        return repr(float(value))
    exponent = math.floor(math.log10(value)) - 3
    return f"{math.ceil(value / 10.0**exponent) * 10.0**exponent:.4g}"


def _write(path: Path, header: list[str], fields: list[str], rows: list[dict]) -> None:
    """Write a derived table with its provenance lines."""
    with path.open("w", newline="") as handle:
        handle.write("\n".join(header) + "\n")
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: repr(float(v)) if isinstance(v, float | np.floating) else v
                    for k, v in row.items()
                }
            )


def rule_rows(analysis: Analysis, final: Final) -> list[dict]:
    """Return one row per node of every band system."""
    lines = analysis.lines
    f, e = oscillator(lines.nu, lines.strength)
    rows = []
    for number, system in enumerate(final.derivation.systems):
        index = system.index
        iso, upper, lower = system.key
        for j, (strength, energy) in enumerate(
            zip(system.strengths, system.energies, strict=True)
        ):
            rows.append(
                {
                    "system": number,
                    "treatment": system.treatment,
                    "global_iso_id": iso,
                    "upper": upper,
                    "lower": lower,
                    "where": system.where,
                    "bands": system.bands,
                    "n_lines": len(index),
                    "min_cm": float(lines.nu[index].min()),
                    "max_cm": float(lines.nu[index].max()),
                    "m0": float(f[index].sum()),
                    "m_minus2_au": float(np.sum(f[index] / e[index] ** 2)),
                    "node": j,
                    "wavenumber_cm": float(energy * HCM),
                    "strength": float(strength),
                }
            )
    return rows


def window_rows(analysis: Analysis, final: Final) -> list[dict]:
    """Return the windows and gaps in ascending wavenumber, with their criteria."""
    grid, at = analysis.grid, analysis.at
    error = model_error(analysis, final.derivation.poles, grid, at)
    lines = analysis.lines
    rows = []
    for a, b in final.windows:
        inside = (grid >= a) & (grid <= b)
        in_lines = (lines.nu >= a) & (lines.nu <= b)
        candidate = next(
            (k for k, (x, y) in enumerate(analysis.candidates) if x <= a and b <= y),
            None,
        )
        rows.append(
            {
                "kind": "window",
                "low_cm": a,
                "high_cm": b,
                "low_um": 1e4 / b,
                "high_um": 1e4 / a,
                "reason": "",
                "max_absorption_m": float(at["absorption"][inside].max()),
                "absorption_at_low_m": float(at["absorption"][inside][0]),
                "absorption_at_high_m": float(at["absorption"][inside][-1]),
                "max_lorentz_peak_m": float(final.peaks[in_lines].max(initial=0.0)),
                "max_far_wing_m": float(at["wing"][inside].max()),
                **_ratios(analysis, error, inside),
                "kept_fraction": ""
                if candidate is None or final.derivation.choice.get(candidate) is None
                else final.derivation.choice[candidate],
                "kept_bands": len(final.derivation.kept.get(candidate, []))
                if candidate is not None
                else 0,
                "folded_bands": final.derivation.folded.get(candidate, 0)
                if candidate is not None
                else 0,
            }
        )
    for gap in final.gaps:
        rows.append(
            {
                "kind": "gap",
                "low_cm": gap["low"],
                "high_cm": gap["high"],
                "low_um": 1e4 / gap["high"],
                "high_um": 1e4 / gap["low"],
                "reason": gap["reason"],
                "causes": gap["causes"],
                "max_absorption_m": gap["max_absorption_m"],
                "kept_bands": gap["kept_bands"],
            }
        )
    return sorted(rows, key=lambda r: r["low_cm"])


def _ratios(
    analysis: Analysis, error: dict, inside: np.ndarray
) -> dict[str, float | str]:
    """Return a window's largest bound over the tolerance, enforced and reported.

    ``max_*`` over the part beyond the optical range, where the tolerance qualifies the
    window; ``optical_*`` over the part inside it, where the bound is reported only.
    Empty where the window has no such part.
    """
    at, grid = analysis.at, analysis.grid
    out: dict[str, float | str] = {}
    strict = enforced(analysis, grid)
    for prefix, part in (("max", inside & strict), ("optical", inside & ~strict)):
        for m in (0, 2):
            ratio = error[m][part] / at[f"tol{m}"][part]
            out[f"{prefix}_b{m}_over_tolerance"] = (
                float(ratio.max()) if part.any() else ""
            )
    return out


WINDOW_FIELDS = [
    "kind",
    "low_cm",
    "high_cm",
    "low_um",
    "high_um",
    "reason",
    "max_absorption_m",
    "causes",
    "absorption_at_low_m",
    "absorption_at_high_m",
    "max_lorentz_peak_m",
    "max_far_wing_m",
    "max_b0_over_tolerance",
    "max_b2_over_tolerance",
    "optical_b0_over_tolerance",
    "optical_b2_over_tolerance",
    "kept_fraction",
    "kept_bands",
    "folded_bands",
]
RULE_FIELDS = [
    "system",
    "treatment",
    "global_iso_id",
    "upper",
    "lower",
    "where",
    "bands",
    "n_lines",
    "min_cm",
    "max_cm",
    "m0",
    "m_minus2_au",
    "node",
    "wavenumber_cm",
    "strength",
]


def absorption_rows(analysis: Analysis) -> list[dict]:
    """Return the largest a_sigma per 1 cm^-1 cell, where above the floor."""
    reference, settings = analysis.reference, analysis.settings
    nu, a = reference.fine_nu, reference.fine_absorption
    cells = np.floor(nu).astype(int)
    top = np.zeros(cells.max() + 1)
    np.maximum.at(top, cells, a)
    return [
        {"low_cm": int(c), "high_cm": int(c + 1), "max_absorption_m": f"{top[c]:.4g}"}
        for c in np.flatnonzero(top >= settings.absorption_floor_m)
    ]


def write(analysis: Analysis, final: Final) -> None:
    """Write the four derived tables of a gas."""
    header = provenance(analysis)
    state = f"{analysis.settings.pressure_pa} Pa, {analysis.settings.temperature_k} K"
    gas = analysis.gas
    _write(
        OUT / f"{gas}_band_rules.csv", header, RULE_FIELDS, rule_rows(analysis, final)
    )
    _write(
        OUT / f"{gas}_windows.csv", header, WINDOW_FIELDS, window_rows(analysis, final)
    )
    _write(
        OUT / f"{gas}_bound.csv",
        [*header, "# cells: geometric, the largest |D_m - D_m ref| in each, in au"],
        ["low_cm", "high_cm", "b0", "b2", "b3"],
        [
            {
                "low_cm": f"{c['low']:.3f}",
                "high_cm": f"{c['high']:.3f}",
                "b0": _up(c["b0"]),
                "b2": _up(c["b2"]),
                "b3": _up(c["b3"]),
            }
            for c in final.cells
        ],
    )
    _write(
        OUT / f"{gas}_absorption.csv",
        [*header, f"# pure {gas} at {state} (the reference state)"],
        ["low_cm", "high_cm", "max_absorption_m"],
        absorption_rows(analysis),
    )


def main(gases: list[str]) -> None:
    """Derive and write the tables of each gas."""
    settings = Settings()
    for gas in gases:
        print(f"{gas}:")
        analysis = analyse(gas, settings)
        final = finalize(analysis, choose(analysis))
        write(analysis, final)
        for a, b in final.windows:
            print(f"  window {1e4 / b:.4f}-{1e4 / a:.4f} um")
        print(f"  {len(final.derivation.poles[0])} poles; checks {final.checks}")


if __name__ == "__main__":
    main(sys.argv[1:] or list(GASES))

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Vibrational (infrared) contributions to the polarizability, from HITRAN.

An infrared-active molecule has, besides its electronic oscillator-strength
distribution, vibrational bands whose strengths add to the polarizability. A band of
integrated intensity ``S`` (cm/molecule, HITRAN's definition, which includes the
lower-state population at 296 K and stimulated emission) is a dipole oscillator of
strength ``f = S / (pi r_e)``, with ``r_e`` the classical electron radius, in the same
units as the electronic oscillator strengths. Each band is represented by one pole whose
energy keeps both the band's total strength and its static contribution:
``E = sqrt(sum S / sum(S / nu^2))``. At optical frequencies the bands add
``-sum f / omega^2`` to the polarizability, small but not negligible in a difference of
polarizabilities (phase matching); at zero frequency they add ``sum f / E^2`` (the
vibrational polarizability).

Pure rotational bands are left out (their static contribution is orientation, not
polarization). The poles describe the index away from the bands. CH4 and N2O also have
weak overtone and combination bands in the near infrared (the tables reach 11000 cm^-1,
0.91 um), where the index has resonant structure that one pole per band does not
represent.

The band tables are ``data/external/hitran/<gas>_bands.csv``, written by
``tools/hitran/band_strengths.py`` (see ``PROVENANCE.md`` there).
"""

import csv
import hashlib
import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache

import numpy as np
from jaxtyping import Float
from scipy import constants as _constants

from optmat.datasets import DATA_DIRECTORY
from optmat.exceptions import DataError
from optmat.gas.gases import Gas, as_gas

__all__ = [
    "DERIVED_TABLES",
    "GAP_COLUMNS",
    "BandSystem",
    "InfraredBand",
    "WindowRecord",
    "band_systems",
    "bound_cells",
    "derived_digest",
    "derived_provenance",
    "fold_poles",
    "gap_cells",
    "gap_digest",
    "gauss_rule",
    "infrared_bands",
    "infrared_polarizability_au",
    "infrared_windows",
    "qualified_windows",
    "static_infrared_polarizability_au",
]

_DIRECTORY = DATA_DIRECTORY / "external" / "hitran"
_ELECTRON_RADIUS_CM: float = (
    _constants.physical_constants["classical electron radius"][0] * 100
)
_HARTREE_WAVENUMBER: float = (
    _constants.physical_constants["hartree-inverse meter relationship"][0] / 100
)


@dataclass(frozen=True)
class InfraredBand:
    """One vibrational band (or cross-section window) of a gas.

    Attributes
    ----------
    label : str
        The isotopologue and the upper and lower vibrational quanta (or the window).
    strength : float
        Integrated intensity ``S``, cm/molecule at 296 K, per molecule of the natural
        isotopic mixture.
    strength_over_nu2 : float
        ``sum S / nu^2``, cm^3/molecule.
    """

    label: str
    strength: float
    strength_over_nu2: float

    @property
    def oscillator_strength(self) -> float:
        """Return the band's dipole oscillator strength ``S / (pi r_e)``.

        Returns
        -------
        float
            Dimensionless.
        """
        return self.strength / (math.pi * _ELECTRON_RADIUS_CM)

    @property
    def wavenumber(self) -> float:
        """Return the pole's wavenumber, ``sqrt(S / sum(S / nu^2))``, cm^-1.

        Returns
        -------
        float
            The wavenumber that keeps the band's static contribution.
        """
        return math.sqrt(self.strength / self.strength_over_nu2)

    @property
    def energy_au(self) -> float:
        """Return the pole's energy, hartree.

        Returns
        -------
        float
            ``wavenumber`` in atomic units.
        """
        return self.wavenumber / _HARTREE_WAVENUMBER


def _label(row: dict[str, str]) -> str:
    """Return a band's label: isotopologue and transition, or the window.

    Parameters
    ----------
    row : dict
        A row of a band table.

    Returns
    -------
    str
        The label.
    """
    if not row["lower"]:
        return f"{row['global_iso_id']} {row['upper']}"
    return f"{row['global_iso_id']} {row['upper']} <- {row['lower']}"


@cache
def infrared_bands(gas: Gas | str) -> tuple[InfraredBand, ...]:
    """Return a gas's curated infrared bands, strongest first.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    tuple of InfraredBand
        Its bands.

    Raises
    ------
    DataError
        If no band table is curated for the gas (it is infrared-inactive, or not yet
        done).
    """
    species = as_gas(gas)
    path = _DIRECTORY / f"{species.value}_bands.csv"
    if not path.exists():
        msg = f"no infrared bands are curated for {species.value}"
        raise DataError(msg)
    with path.open(newline="") as handle:
        return tuple(
            InfraredBand(
                label=_label(row),
                strength=float(row["strength_cm_per_molecule"]),
                strength_over_nu2=float(row["strength_over_nu2_cm3_per_molecule"]),
            )
            for row in csv.DictReader(handle)
            # Pure rotational bands (the same vibrational state above and below) are
            # left out: their static contribution is the orientation term mu^2/3kT of a
            # polar molecule, not a polarizability (N2O: 1.45 a.u. at 296 K, TC-156),
            # and at optical frequencies they are negligible (f ~ 1e-8).
            if not (row["lower"] and row["upper"] == row["lower"])
        )


def infrared_polarizability_au(
    gas: Gas | str, omega: Float[np.ndarray, "*shape"]
) -> Float[np.ndarray, "*shape"]:
    """Return the vibrational contribution to the polarizability.

    Parameters
    ----------
    gas : Gas or str
        The species.
    omega : numpy.ndarray
        Photon energy, hartree (away from the bands themselves).

    Returns
    -------
    numpy.ndarray
        ``sum_b f_b / (E_b^2 - omega^2)``, atomic units.
    """
    bands = infrared_bands(gas)
    strengths = np.array([band.oscillator_strength for band in bands])
    energies = np.array([band.energy_au for band in bands])
    omega = np.asarray(omega, dtype=float)
    terms = strengths / (energies**2 - omega[..., None] ** 2)
    return np.asarray(np.sum(terms, axis=-1))


def static_infrared_polarizability_au(gas: Gas | str) -> float:
    """Return the vibrational polarizability at zero frequency.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    float
        ``sum_b f_b / E_b^2``, atomic units.

    Examples
    --------
    >>> round(static_infrared_polarizability_au("SF6"), 1)
    13.9
    """
    return float(infrared_polarizability_au(gas, np.zeros(1))[0])


def gauss_rule(
    strengths: Float[np.ndarray, " m"], energies: Float[np.ndarray, " m"], n: int
) -> tuple[Float[np.ndarray, " k"], Float[np.ndarray, " k"]]:
    r"""Return the n-node Gaussian quadrature of a band, as poles.

    A band of lines of oscillator strengths :math:`f_\ell` at energies :math:`E_\ell`
    is a positive measure with weights :math:`w_\ell = f_\ell / E_\ell^2` at
    :math:`x_\ell = E_\ell^2`; its contribution to the polarizability,
    :math:`\sum_\ell f_\ell / (E_\ell^2 - \omega^2)`, is the measure's Stieltjes
    transform at :math:`s = \omega^2` (OEP 3, section 4). The rule's nodes and weights
    are the eigenvalues and the squared first eigenvector components, times the total
    weight, of the measure's Jacobi matrix (Golub and Welsch 1969, eq. 2.6), which the
    Lanczos process builds from the lines. The rule reproduces the moments
    :math:`\mu_0, \ldots, \mu_{2n-1}`: the band's static contribution
    :math:`M_{-2} = \mu_0`, its optical-side limit :math:`M_0 = \mu_1`, and
    :math:`M_2, \ldots, M_{4n-4}`. One node is the single pole per band of
    :func:`infrared_bands`.

    Parameters
    ----------
    strengths : numpy.ndarray
        The lines' oscillator strengths, dimensionless, > 0.
    energies : numpy.ndarray
        Their energies, hartree, > 0.
    n : int
        The largest number of nodes, >= 1. A band with fewer distinct lines gets its
        exact rule with fewer nodes.

    Returns
    -------
    tuple of numpy.ndarray
        The poles' strengths ``W_j X_j`` and energies ``X_j**0.5``, ascending in
        energy, every strength positive and every energy inside the lines' span.

    Raises
    ------
    ValueError
        If ``n < 1``, there are no lines, or a strength or energy is not positive.

    Examples
    --------
    >>> f, e = gauss_rule(np.array([1.0, 1.0]), np.array([1.0, 2.0]), 1)
    >>> round(float(f[0]), 12), round(float(e[0]), 12)  # E = (M0 / M_-2)^(1/2)
    (2.0, 1.264911064067)
    """
    f = np.asarray(strengths, dtype=float)
    e = np.asarray(energies, dtype=float)
    if n < 1 or f.size == 0 or np.any(f <= 0) or np.any(e <= 0):
        msg = "a rule needs n >= 1 and lines of positive strength and energy"
        raise ValueError(msg)
    x = e**2
    scale = float(x.max())
    x = x / scale
    w = f / e**2
    total = float(w.sum())
    q = np.sqrt(w / total)
    basis = [q]
    diagonal: list[float] = []
    off: list[float] = []
    previous, beta = np.zeros_like(q), 0.0
    while True:  # at most n steps: the loop ends at n nodes or an exhausted space
        v = x * basis[-1] - beta * previous
        alpha = float(basis[-1] @ v)
        diagonal.append(alpha)
        v = v - alpha * basis[-1]
        for u in basis:  # full reorthogonalisation keeps the recurrence exact
            v = v - (u @ v) * u
        beta = float(np.linalg.norm(v))
        # The Krylov space is exhausted: the measure has no more distinct nodes.
        if len(diagonal) == n or beta <= 1e-13:
            break
        off.append(beta)
        previous = basis[-1]
        basis.append(v / beta)
    jacobi = np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1)
    nodes, vectors = np.linalg.eigh(jacobi)
    weights = total * vectors[0] ** 2
    nodes = nodes * scale
    return np.asarray(weights * nodes), np.asarray(np.sqrt(nodes))


def fold_poles(
    m0: float, m_minus2: float, below: float, above: float
) -> tuple[float, float]:
    r"""Return two pole strengths that keep a band's :math:`M_0` and :math:`M_{-2}`.

    A weak band inside a window is folded onto two poles outside it, at the nearest
    rule nodes below and above the window (OEP 3, section 4). The strengths solve
    :math:`f_1 + f_2 = M_0` and :math:`f_1/E_1^2 + f_2/E_2^2 = M_{-2}`; both are
    positive because the band's :math:`\bar E = (M_0/M_{-2})^{1/2}` lies between
    the two energies.

    Parameters
    ----------
    m0 : float
        The band's summed oscillator strength, > 0.
    m_minus2 : float
        Its summed ``f / E**2``, hartree^-2, > 0.
    below : float
        The lower pole's energy, hartree.
    above : float
        The upper pole's energy, hartree.

    Returns
    -------
    tuple of float
        The strengths at ``below`` and ``above``.

    Raises
    ------
    ValueError
        If :math:`\bar E` does not lie strictly between the two energies.

    Examples
    --------
    >>> fold_poles(2.0, 1.25, 1.0, 2.0)
    (1.0, 1.0)
    """
    mean = math.sqrt(m0 / m_minus2) if m0 > 0 and m_minus2 > 0 else math.nan
    if not below < mean < above:
        msg = f"the band's mean energy {mean} is not between {below} and {above}"
        raise ValueError(msg)
    inverse = 1 / below**2 - 1 / above**2
    lower = (m_minus2 - m0 / above**2) / inverse
    return lower, m0 - lower


@dataclass(frozen=True)
class BandSystem:
    """A band system of a gas's derived rules (OEP 3), and its poles.

    Attributes
    ----------
    system : int
        Its number in the table.
    treatment : str
        ``rule`` (a band, or its part between two windows), ``rotational`` (every pure
        rotational line) or ``folded`` (the weak bands inside one window, on two poles
        outside it).
    label : str
        The isotopologue and the transition, or the treatment.
    where : str
        The gap or window of the derivation its lines lie in.
    bands : int
        How many bands (or band parts) it represents.
    n_lines : int
        Its number of lines.
    span_cm : tuple of float
        Its lines' wavenumber range, cm^-1.
    m0 : float
        Its lines' summed oscillator strength.
    m_minus2_au : float
        Their summed ``f / E**2``, hartree^-2: the static contribution.
    strengths : tuple of float
        Its poles' strengths.
    energies_au : tuple of float
        Its poles' energies, hartree.
    """

    system: int
    treatment: str
    label: str
    where: str
    bands: int
    n_lines: int
    span_cm: tuple[float, float]
    m0: float
    m_minus2_au: float
    strengths: tuple[float, ...]
    energies_au: tuple[float, ...]


@dataclass(frozen=True)
class WindowRecord:
    """A qualified window of a gas, or a gap between windows, from the derivation.

    Attributes
    ----------
    kind : str
        ``window`` or ``gap``.
    low_cm, high_cm : float
        Its wavenumber range, cm^-1.
    reason : str
        For a gap: ``absorption``, ``representation``, ``representation: kept band``
        or ``width``.
    max_absorption_m : float
        The largest smoothed absorption of the pure gas at the reference state in it,
        m^-1.
    details : Mapping of str to str
        The table's other columns: the criteria at the edges, the largest single-line
        Lorentz peak and far-wing absorption, the bound over the tolerance, the kept
        and folded bands.
    """

    kind: str
    low_cm: float
    high_cm: float
    reason: str
    max_absorption_m: float
    details: Mapping[str, str]


def _derived(gas: Gas | str, table: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a derived table of a gas: its ``#`` provenance lines and its rows.

    Parameters
    ----------
    gas : Gas or str
        The species.
    table : str
        ``band_rules``, ``windows``, ``bound``, ``absorption`` or ``gap_bound``.

    Returns
    -------
    tuple
        The provenance lines and the rows.

    Raises
    ------
    DataError
        If the gas has no such table.
    """
    species = as_gas(gas)
    path = _DIRECTORY / f"{species.value}_{table}.csv"
    if not path.exists():
        msg = f"no derived {table} table for {species.value}"
        raise DataError(msg)
    text = path.read_text().splitlines()
    header = [line for line in text if line.startswith("#")]
    rows = list(csv.DictReader(line for line in text if not line.startswith("#")))
    return header, rows


def derived_provenance(gas: Gas | str) -> tuple[str, ...]:
    """Return the provenance lines of a gas's derived tables.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    tuple of str
        The ``#`` lines: the tool and its commit, the raw inputs' SHA-256 and query,
        the parent fit and every setting.
    """
    return tuple(_derived(gas, "band_rules")[0])


@cache
def band_systems(gas: Gas | str) -> tuple[BandSystem, ...]:
    """Return a gas's band systems and their rules (OEP 3).

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    tuple of BandSystem
        In the table's order.
    """
    grouped: dict[int, list[dict[str, str]]] = {}
    for row in _derived(gas, "band_rules")[1]:
        grouped.setdefault(int(row["system"]), []).append(row)
    out = []
    for number, rows in grouped.items():
        first = rows[0]
        label = (
            f"{first['global_iso_id']} {first['upper']} <- {first['lower']}"
            if first["treatment"] == "rule"
            else first["treatment"]
        )
        out.append(
            BandSystem(
                system=number,
                treatment=first["treatment"],
                label=label,
                where=first["where"],
                bands=int(first["bands"]),
                n_lines=int(first["n_lines"]),
                span_cm=(float(first["min_cm"]), float(first["max_cm"])),
                m0=float(first["m0"]),
                m_minus2_au=float(first["m_minus2_au"]),
                strengths=tuple(float(r["strength"]) for r in rows),
                energies_au=tuple(
                    float(r["wavenumber_cm"]) / _HARTREE_WAVENUMBER for r in rows
                ),
            )
        )
    return tuple(out)


@cache
def infrared_windows(gas: Gas | str) -> tuple[WindowRecord, ...]:
    """Return a gas's derived windows and the gaps between them, ascending (OEP 3).

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    tuple of WindowRecord
        Every window and gap from 12.5 um to the parent fit's shortest wavelength.
    """
    keys = {"kind", "low_cm", "high_cm", "reason", "max_absorption_m"}
    return tuple(
        WindowRecord(
            kind=row["kind"],
            low_cm=float(row["low_cm"]),
            high_cm=float(row["high_cm"]),
            reason=row["reason"],
            max_absorption_m=float(row["max_absorption_m"]),
            details={k: v for k, v in row.items() if k not in keys and v},
        )
        for row in _derived(gas, "windows")[1]
    )


@cache
def bound_cells(gas: Gas | str) -> dict[str, Float[np.ndarray, " n"]]:
    """Return a gas's representation bound per cell of its windows (OEP 3).

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    dict of str to numpy.ndarray
        ``low_cm`` and ``high_cm``, the cells' edges (cm^-1), and ``b0``, ``b2`` and
        ``b3``, the bound on the magnitude of D_m[model] - D_m[reference] in each
        (atomic units).
    """
    rows = _derived(gas, "bound")[1]
    return {
        key: np.array([float(r[key]) for r in rows])
        for key in ("low_cm", "high_cm", "b0", "b2", "b3")
    }


GAP_COLUMNS: tuple[str, ...] = ("low_cm", "high_cm", "b0", "b2", "b3", "absorption_m")
"""The numeric columns of a gas's gap table."""


@cache
def gap_cells(gas: Gas | str) -> dict[str, Float[np.ndarray, " n"]] | None:
    """Return a fit's representation error and absorption in its gaps (OEP 3, Q13).

    Written by ``tools/hitran/gap_bound.py`` for the stored fit, against the same
    smoothed line-by-line reference as :func:`bound_cells`, in the fit's gaps from
    12.5 um to the top of the absorption grid, where a mixture holding the gas as a
    trace may still be qualified.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    dict of str to numpy.ndarray or None
        ``low_cm`` and ``high_cm``, the cells' edges (cm^-1); ``b0``, ``b2`` and
        ``b3``, the largest magnitude of D_m[fit] - D_m[reference] in each (atomic
        units, infinite in a cell holding a pole of the fit); ``absorption_m``, the
        pure gas's largest smoothed absorption at 101325 Pa and 296 K (m^-1); and
        ``pole``, True where a pole lies in the cell. None for a gas without the
        table.
    """
    try:
        rows = _derived(gas, "gap_bound")[1]
    except DataError:
        return None
    cells = {key: np.array([float(r[key]) for r in rows]) for key in GAP_COLUMNS}
    cells["pole"] = np.array([r["pole"] == "1" for r in rows])
    return cells


def gap_digest(gas: Gas | str) -> str:
    """Return the SHA-256 of a gas's gap table.

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    str
        The digest of ``<gas>_gap_bound.csv``.
    """
    path = _DIRECTORY / f"{as_gas(gas).value}_gap_bound.csv"
    return hashlib.sha256(path.read_bytes()).hexdigest()


DERIVED_TABLES: tuple[str, ...] = ("band_rules", "windows", "bound", "absorption")
"""The derived tables of a gas, in the order :func:`derived_digest` hashes them."""


def derived_digest(gas: Gas | str) -> str:
    """Return the SHA-256 of a gas's derived tables, as a fit records them.

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    str
        The digest of the four tables' bytes, in :data:`DERIVED_TABLES` order.

    Raises
    ------
    DataError
        If a table is missing.
    """
    species = as_gas(gas)
    digest = hashlib.sha256()
    for table in DERIVED_TABLES:
        path = _DIRECTORY / f"{species.value}_{table}.csv"
        if not path.exists():
            msg = f"no derived {table} table for {species.value}"
            raise DataError(msg)
        digest.update(path.read_bytes())
    return digest.hexdigest()


def qualified_windows(
    gas: Gas | str,
) -> tuple[tuple[tuple[float, float], ...], tuple[str, ...]]:
    """Return a gas's windows as qualified wavelength intervals, and its gaps' reasons.

    Parameters
    ----------
    gas : Gas or str
        The species: CO, CO2, CH4 or N2O.

    Returns
    -------
    tuple
        The intervals, vacuum wavelength in nm, ascending; and the reason of each gap
        between two neighbouring intervals, in the same order.
    """
    records = infrared_windows(gas)
    intervals: list[tuple[float, float]] = []
    reasons: list[str] = []
    pending = ""
    # Descending wavenumber is ascending wavelength; a gap counts only between two
    # windows, not beyond the longest one.
    for record in reversed(records):
        if record.kind == "window":
            if intervals:
                reasons.append(pending)
            intervals.append((1e7 / record.high_cm, 1e7 / record.low_cm))
        else:
            pending = record.reason
    return tuple(intervals), tuple(reasons)

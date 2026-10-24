# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The second hyperpolarizability: a process-aware dispersion model and its fit.

The third-order response of a gas molecule below resonance is its second
hyperpolarizability :math:`\gamma(-\omega_\sigma; \omega_1, \omega_2, \omega_3)`, which
depends on all four frequencies. Here it is expressed in the Taylor-series ("T",
Buckingham-Orr) convention of Shelton & Rice (Chem. Rev. 94, 3 (1994)),
:math:`p = \dots + \gamma E^3/6`, in atomic units (:data:`AU_HYPERPOLARIZABILITY`).

Bishop showed, and Shelton's data confirm, that to fourth order in frequency the
dispersion depends on the frequencies only through two symmetric sums,

.. math::

   \gamma = \gamma_0 \left[1 + A\,W_2 + B\,W_2^2 + B'\,W_4\right], \qquad
   W_2 = \nu_\sigma^2 + \nu_1^2 + \nu_2^2 + \nu_3^2, \quad
   W_4 = \nu_\sigma^4 + \nu_1^4 + \nu_2^4 + \nu_3^4,

with :math:`\nu` the wavenumbers. :math:`A` is the same for every process; :math:`B'`
separates processes at fourth order (Bishop & De Kee 1996). At a fundamental wavenumber
:math:`\nu`:

- the Kerr effect (the intensity-dependent index, SPM and XPM), :math:`\gamma(-\omega;
  \omega, \omega, -\omega)`, has :math:`W_2 = 4\nu^2`;
- third-harmonic generation, :math:`\gamma(-3\omega; \omega, \omega, \omega)`, has
  :math:`W_2 = 12\nu^2`;
- ESHG :math:`6\nu^2`, dc Kerr :math:`2\nu^2`;
- a pump-probe index at pump :math:`\nu_p` and probe :math:`\nu_r` has
  :math:`2\nu_p^2 + 2\nu_r^2`.

The data enter in their sources' own quantities and conventions (``docs/data.md``,
*Third-order quantities*), converted here:

- ``T`` in au, SI (C^4 m^4 J^-3) or esu;
- Ward's per-molecule chi in esu, gamma_T / 6 (``ward_1969``);
- Lehmeier's gamma, gamma_T / 24 (``lehmeier_1985``);
- CARS chi per amagat, gamma_T = 2 chi / N_0 (``chi_cars_per_amagat``);
- n2, gamma_T = 8 eps0^2 c n2 / N for the self coefficient at the dataset's density
  (half the probe coefficient).

The model is electronic. Data that also carry the vibrational and rotational
contributions (Shelton's H2, D2, N2 and O2 ESHG) are corrected by subtracting Shelton's
gamma^vR, as a dataset of the same rows or as his ``G nu^-2 + H nu^-4`` formula.
"""

import dataclasses
import functools
import hashlib
import json
import math
import tomllib
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar, Final

import numpy as np
from jaxtyping import Float

from optmat import _catalogue
from optmat._catalogue import records
from optmat.conditions import State
from optmat.constants import (
    AU_ENERGY,
    BOHR_RADIUS,
    E_CHARGE,
    EPS0,
    K_B,
    C,
)
from optmat.datasets import Dataset
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY, reduction
from optmat.gas.eos import number_density
from optmat.gas.gases import Gas, as_gas, info
from optmat.identity import Identified, source_ids

__all__ = [
    "AU_HYPERPOLARIZABILITY",
    "ESU_PER_SI",
    "EVALUATION_VERSION",
    "NONLINEAR_DIRECTORY",
    "SUPPORT_DERIVATION_VERSION",
    "SUPPORT_LIMIT_KINDS",
    "W2_MARGIN",
    "AbsorptionEvidence",
    "HyperpolarizabilityModel",
    "HyperpolarizabilitySupport",
    "NonlinearConfig",
    "NonlinearFit",
    "NonlinearObservation",
    "NonlinearOptions",
    "SupportCandidate",
    "SupportLimit",
    "SupportPoint",
    "SupportUncertainty",
    "dumps",
    "fit_hyperpolarizability",
    "frequency_sums",
    "gamma_observation",
    "hyperpolarizability_models",
    "load_nonlinear_config",
    "validation",
]

EVALUATION_VERSION: Final = 1
"""The version of this module's evaluation code, in the identity of every model that
evaluates through it; raised whenever the code changes its output (OEP 1, section 4.8).
"""

AU_HYPERPOLARIZABILITY: Final = E_CHARGE**4 * BOHR_RADIUS**4 / AU_ENERGY**3
"""One atomic unit of the second hyperpolarizability, C^4 m^4 J^-3 (6.235e-65)."""

ESU_PER_SI: Final = C**4 * 1e-9
"""esu (statC cm^4 statV^-3) per SI unit (C m^4 V^-3) of gamma: 8.0776e24."""

W2_MARGIN: Final = 1.1
"""How far beyond the largest fitted W2 a process may be evaluated without a warning."""

NONLINEAR_DIRECTORY: Final = GAS_DATA_DIRECTORY / "nonlinear"
"""Fit configurations of the hyperpolarizability, ``<gas>.toml``."""

_MODEL_DIRECTORY: Final = GAS_DATA_DIRECTORY / "models" / "nonlinear"
_LOSCHMIDT_CM3: Final = 101325.0 / (K_B * 273.15) * 1e-6
# Factor from each (convention, unit) to gamma_T in atomic units.
_TO_AU: Final[dict[tuple[str, str], float]] = {
    ("T", "au"): 1.0,
    ("T", "SI"): 1.0 / AU_HYPERPOLARIZABILITY,
    ("T", "esu"): 1.0 / (AU_HYPERPOLARIZABILITY * ESU_PER_SI),
    ("ward_1969", "esu"): 6.0 / (AU_HYPERPOLARIZABILITY * ESU_PER_SI),
    ("lehmeier_1985", "SI"): 24.0 / AU_HYPERPOLARIZABILITY,
    ("lehmeier_1985", "esu"): 24.0 / (AU_HYPERPOLARIZABILITY * ESU_PER_SI),
}
_MEAN_COMPONENTS: Final = frozenset({"parallel", "mean", "None"})


def frequency_sums(
    process: str,
    wavenumber: Float[np.ndarray, "*shape"],
    other: Float[np.ndarray, "*shape"] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the dispersion variables ``W2`` and ``W4`` of a process.

    Parameters
    ----------
    process : str
        ``static``, ``kerr``, ``eshg``, ``dc_kerr``, ``thg``, ``cars`` or
        ``pump_probe``.
    wavenumber : numpy.ndarray
        The fundamental (or pump) vacuum wavenumber, cm^-1.
    other : numpy.ndarray, optional
        The Stokes wavenumber (``cars``) or the probe wavenumber (``pump_probe``).

    Returns
    -------
    tuple of numpy.ndarray
        ``W2`` (cm^-2) and ``W4`` (cm^-4).

    Raises
    ------
    DataError
        For an unknown process, or ``cars`` and ``pump_probe`` without ``other``.

    Examples
    --------
    >>> w2, w4 = frequency_sums("thg", np.array([12500.0]))
    >>> float(w2[0] / 12500.0**2), float(w4[0] / 12500.0**4)
    (12.0, 84.0)
    """
    nu = np.asarray(wavenumber, dtype=float)
    zero = np.zeros_like(nu)
    if process in {"cars", "pump_probe"} and other is None:
        msg = f"{process} needs a second frequency"
        raise DataError(msg)
    table = {
        "static": (zero, zero, zero, zero),
        "kerr": (nu, nu, nu, nu),
        "eshg": (2 * nu, nu, nu, zero),
        "dc_kerr": (nu, nu, zero, zero),
        "thg": (3 * nu, nu, nu, nu),
    }
    if process == "cars":
        assert other is not None
        stokes = np.asarray(other, dtype=float)
        terms = (2 * nu - stokes, nu, nu, stokes)
    elif process == "pump_probe":
        assert other is not None
        probe = np.asarray(other, dtype=float)
        terms = (probe, probe, nu, nu)
    elif process in table:
        terms = table[process]
    else:
        msg = f"unknown third-order process {process!r}"
        raise DataError(msg)
    return (
        np.asarray(sum(t**2 for t in terms), dtype=float),
        np.asarray(sum(t**4 for t in terms), dtype=float),
    )


@dataclass(frozen=True)
class NonlinearOptions:
    """How one dataset enters a hyperpolarizability fit.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    free_scale : bool
        Fit the dataset's scale without a prior (its shape only).
    systematic : float
        A relative scale uncertainty (one sigma): the scale is fitted with this prior.
    extra_relative : float
        Relative uncertainty added in quadrature to every value.
    range_nm : tuple of float
        Fundamental-wavelength window of the values used, nm.
    subtract_formula : str
        ``source/formula`` of Shelton's ``gamma^vR = G nu^-2 + H nu^-4`` (C^4 m^4 J^-3)
        to subtract, or "".
    subtract_dataset : str
        ``source/dataset`` of gamma^vR on the same wavelengths to subtract, or "".
    scale_group : str
        Datasets naming the same group share one scale (a calculation's processes, whose
        relative dispersion informs B' whatever its level); "" for none.
    rows : tuple of str
        Only these rows (one method of a table of methods, say); empty for all.
    factor : float
        A factor applied to the converted values and uncertainties, for a unit the
        dataset's own conversion does not know (Bishop & Pipin's reduced-mass atomic
        unit of helium, 1.00096 standard units); its reason is stated in the
        configuration.
    """

    key: str
    free_scale: bool = False
    systematic: float = 0.0
    extra_relative: float = 0.0
    range_nm: tuple[float, float] = (0.0, math.inf)
    subtract_formula: str = ""
    subtract_dataset: str = ""
    scale_group: str = ""
    rows: tuple[str, ...] = ()
    factor: float = 1.0


@dataclass(frozen=True)
class NonlinearObservation:
    """One dataset's values as electronic gamma_T in atomic units.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    wavelength_nm : numpy.ndarray
        The fundamental (pump) vacuum wavelength, nm (``inf`` for static).
    w2, w4 : numpy.ndarray
        The dispersion variables, cm^-2 and cm^-4.
    gamma : numpy.ndarray
        gamma_T, atomic units.
    sigma : numpy.ndarray
        One-sigma uncertainty, atomic units (NaN where none is stated).
    process : str
        The process (``pump_probe`` for a nondegenerate index).
    options : NonlinearOptions
        How the dataset enters the fit.
    """

    key: str
    wavelength_nm: Float[np.ndarray, " n"]
    w2: Float[np.ndarray, " n"]
    w4: Float[np.ndarray, " n"]
    gamma: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    process: str
    options: NonlinearOptions


def _formula(key: str) -> Mapping[str, object]:
    """Return a source's formula record by ``source/formula``.

    Parameters
    ----------
    key : str
        ``source/formula`` id.

    Returns
    -------
    mapping
        The ``[[formulas]]`` entry.

    Raises
    ------
    DataError
        If there is no such formula.
    """
    source_id, _, formula_id = key.partition("/")
    source = records().get(source_id)
    for entry in source.formulas if source else ():
        if entry.get("id") == formula_id:
            return entry
    msg = f"no formula {key!r}"
    raise DataError(msg)


def _wavenumber(dataset: Dataset, nm: float | None) -> float | None:
    """Return the vacuum wavenumber (cm^-1) of a wavelength parameter in nm.

    Parameters
    ----------
    dataset : Dataset
        Whose medium applies.
    nm : float or None
        The wavelength, nm.

    Returns
    -------
    float or None
        The wavenumber, or None.
    """
    if nm is None:
        return None
    wavelength = float(nm) * 1e-9
    if dataset.wavelength_medium == "standard_air":
        # A visible line: the air-vacuum difference (about 3e-4) is immaterial here.
        wavelength *= 1.000277
    return 1.0 / (wavelength * 100.0)


def _conversion(
    dataset: Dataset, nu: np.ndarray
) -> tuple[float, str, np.ndarray | None]:
    """Return a dataset's factor to gamma_T (au), its process and second frequency.

    Parameters
    ----------
    dataset : Dataset
        A third-order dataset.
    nu : numpy.ndarray
        Its fundamental wavenumbers, cm^-1.

    Returns
    -------
    tuple
        The factor, the process and the Stokes or probe wavenumbers (or None).

    Raises
    ------
    DataError
        If the dataset's quantity, convention or unit has no conversion.
    """
    parameters = dataset.parameters
    process = str(parameters.get("process", "kerr"))
    other = None
    if dataset.quantity == "hyperpolarizability":
        convention = str(parameters.get("convention"))
        factor = _TO_AU.get((convention, dataset.unit))
        if factor is None:
            msg = f"{dataset.key}: no conversion for ({convention!r}, {dataset.unit!r})"
            raise DataError(msg)
    elif dataset.quantity == "nonlinear_susceptibility":
        if parameters.get("convention") != "chi_cars_per_amagat":
            msg = f"{dataset.key}: only per-amagat CARS chi(3) is converted"
            raise DataError(msg)
        factor = 2.0 / _LOSCHMIDT_CM3 / (AU_HYPERPOLARIZABILITY * ESU_PER_SI)
    elif dataset.quantity == "nonlinear_index":
        # Dataset's quantity validation guarantees one of these two units.
        per_m2 = {"cm2/W": 1e-4, "m2/W": 1.0}[dataset.unit]
        half = 0.5 if parameters.get("coefficient") == "probe" else 1.0
        density = reduction.density_at_state(dataset)
        factor = 8 * EPS0**2 * C * per_m2 * half / density / AU_HYPERPOLARIZABILITY
        process = "kerr"
        probe = _wavenumber(dataset, parameters.get("probe_nm"))  # pyright: ignore[reportArgumentType]
        if probe is not None:
            process, other = "pump_probe", np.full_like(nu, probe)
    else:
        msg = f"{dataset.key}: {dataset.quantity} is not converted to gamma"
        raise DataError(msg)
    if process == "cars":
        stokes = _wavenumber(dataset, parameters.get("stokes_nm"))  # pyright: ignore[reportArgumentType]
        if stokes is None:
            msg = f"{dataset.key}: a CARS value needs parameters.stokes_nm"
            raise DataError(msg)
        other = np.full_like(nu, stokes)
    return factor, process, other


def _vibrational_rotational(
    options: NonlinearOptions, nu: np.ndarray, wavelength: np.ndarray
) -> np.ndarray:
    """Return the gamma^vR (au) an option subtracts, at each value's wavelength.

    Parameters
    ----------
    options : NonlinearOptions
        The dataset's options.
    nu : numpy.ndarray
        Fundamental wavenumbers, cm^-1.
    wavelength : numpy.ndarray
        Fundamental wavelengths, nm.

    Returns
    -------
    numpy.ndarray
        The correction, zero without one.
    """
    correction = np.zeros_like(nu)
    if options.subtract_formula:
        coefficients = _formula(options.subtract_formula)["coefficients"]
        g = float(coefficients["G"])  # pyright: ignore[reportIndexIssue]
        h = float(coefficients["H"])  # pyright: ignore[reportIndexIssue]
        with np.errstate(divide="ignore"):
            formula = np.where(nu > 0, g / nu**2 + h / nu**4, 0.0)
        correction = correction + formula / AU_HYPERPOLARIZABILITY
    if options.subtract_dataset:
        vr = _catalogue.dataset(options.subtract_dataset)
        lookup = dict(
            zip(
                np.round(vr.vacuum_wavelength() * 1e9, 2),
                gamma_observation(vr, NonlinearOptions(vr.key)).gamma,
                strict=True,
            )
        )
        correction = correction + np.array([lookup[w] for w in np.round(wavelength, 2)])
    return correction


def _kept(
    dataset: Dataset, options: NonlinearOptions, wavelength: np.ndarray
) -> np.ndarray:
    """Return which of a dataset's rows its options admit: the window and the rows.

    Parameters
    ----------
    dataset : Dataset
        The dataset.
    options : NonlinearOptions
        Its window and row selection.
    wavelength : numpy.ndarray
        Each row's fundamental vacuum wavelength, nm.

    Returns
    -------
    numpy.ndarray
        A boolean mask over the rows.

    Raises
    ------
    DataError
        If the options name a row the dataset does not have.
    """
    keep = (wavelength >= options.range_nm[0]) & (wavelength <= options.range_nm[1])
    if options.rows:
        unknown = set(options.rows) - set(dataset.rows)
        if unknown:
            msg = f"{dataset.key}: unknown rows {sorted(unknown)}"
            raise DataError(msg)
        keep &= np.isin(np.array(dataset.rows), np.array(options.rows))
    return keep


def gamma_observation(
    dataset: Dataset, options: NonlinearOptions
) -> NonlinearObservation:
    """Convert a third-order dataset to gamma_T in atomic units.

    Parameters
    ----------
    dataset : Dataset
        A ``hyperpolarizability``, ``nonlinear_index`` or per-amagat CARS
        ``nonlinear_susceptibility`` dataset.
    options : NonlinearOptions
        Its window and corrections.

    Returns
    -------
    NonlinearObservation
        The values, with their dispersion variables.

    Raises
    ------
    DataError
        If the dataset's quantity, convention or unit has no conversion, or its
        wavelength is not stated.
    """
    parameters = dataset.parameters
    component = str(parameters.get("component", "None"))
    # For an atom the laboratory zzzz component is the mean gamma.
    atomic = info(dataset.material).atomic and component == "zzzz"
    # Buckingham's Kerr gamma_K of a dc Kerr measurement is the mean gamma up to its
    # dispersion-induced anisotropy: 0.16 % for He and about 0.4 % for H2 at 633 nm,
    # below the measurements' 2-3 % (TC-401).
    kerr = component == "K" and parameters.get("process") == "dc_kerr"
    if component not in _MEAN_COMPONENTS and not atomic and not kerr:
        msg = f"{dataset.key}: component {component!r} is not the mean gamma"
        raise DataError(msg)
    wavelength = dataset.vacuum_wavelength() * 1e9
    if np.any(np.isnan(wavelength)):
        msg = f"{dataset.key}: the wavelength is not stated"
        raise DataError(msg)
    nu = np.where(np.isfinite(wavelength), 1e7 / wavelength, 0.0)
    values = dataset.scaled_values()
    sigma = dataset.standard_uncertainty()
    factor, process, other = _conversion(dataset, nu)
    gamma = values * factor * options.factor
    sigma = sigma * factor * options.factor
    gamma = gamma - _vibrational_rotational(options, nu, wavelength)
    keep = _kept(dataset, options, wavelength)
    w2, w4 = frequency_sums(process, nu, other)
    extra = options.extra_relative * np.abs(gamma)
    sigma = np.where(np.isfinite(sigma), np.hypot(sigma, extra), extra)
    return NonlinearObservation(
        key=dataset.key,
        wavelength_nm=wavelength[keep],
        w2=w2[keep],
        w4=w4[keep],
        gamma=gamma[keep],
        sigma=sigma[keep],
        process=process,
        options=options,
    )


@dataclass(frozen=True)
class NonlinearConfig:
    """A hyperpolarizability fit's specification (``data/gas/nonlinear/<gas>.toml``).

    Attributes
    ----------
    gas : Gas
        The species.
    description : str
        What the fit is for.
    range_nm : tuple of float
        The fundamental-wavelength range over which the model is declared valid, nm.
    b_prime_prior : tuple of float or None
        Prior ``(value, sigma)`` on B' (1e-20 cm^4); None fixes B' at 0.
    a_prior, b_prior : tuple of float or None
        Optional priors ``(value, sigma)`` on A (1e-10 cm^2) and B (1e-20 cm^4): a
        dispersion borrowed from a similar molecule where the gas's own data cannot fix
        it.
    fitted : tuple of NonlinearOptions
        The datasets fitted.
    validators : tuple of NonlinearOptions
        Datasets compared with the fit.
    digest : str
        SHA-256 of the configuration file.
    """

    gas: Gas
    description: str
    range_nm: tuple[float, float]
    b_prime_prior: tuple[float, float] | None
    fitted: tuple[NonlinearOptions, ...]
    validators: tuple[NonlinearOptions, ...]
    digest: str
    a_prior: tuple[float, float] | None = None
    b_prior: tuple[float, float] | None = None


def _pair(value: object) -> tuple[float, float] | None:
    """Read an optional ``[value, sigma]`` prior.

    Parameters
    ----------
    value : object
        The TOML value, or None.

    Returns
    -------
    tuple of float or None
        The prior.
    """
    if value is None:
        return None
    first, second = value  # pyright: ignore[reportGeneralTypeIssues]
    return float(first), float(second)


def _options(entry: object, where: str) -> NonlinearOptions:
    """Read one ``[[fitted]]`` entry or validator.

    Parameters
    ----------
    entry : object
        A key string or a table.
    where : str
        For error messages.

    Returns
    -------
    NonlinearOptions
        The options.

    Raises
    ------
    DataError
        If the entry has unknown keys or names an unknown dataset.
    """
    if isinstance(entry, str):
        entry = {"key": entry}
    if not isinstance(entry, Mapping):
        msg = f"{where}: a dataset entry must be a key or a table"
        raise DataError(msg)
    fields = {f.name for f in dataclasses.fields(NonlinearOptions)}
    unknown = set(entry) - fields
    if unknown:
        msg = f"{where}: unknown keys {sorted(unknown)}"
        raise DataError(msg)
    key = str(entry["key"])
    if key not in _catalogue.by_key():
        msg = f"{where}: unknown dataset {key!r}"
        raise DataError(msg)
    window = entry.get("range_nm", [0.0, "inf"])
    low, high = (float(v) for v in window)  # pyright: ignore[reportGeneralTypeIssues]
    return NonlinearOptions(
        key=key,
        free_scale=bool(entry.get("free_scale", False)),
        systematic=float(entry.get("systematic", 0.0)),  # pyright: ignore[reportArgumentType]
        extra_relative=float(entry.get("extra_relative", 0.0)),  # pyright: ignore[reportArgumentType]
        range_nm=(low, high),
        subtract_formula=str(entry.get("subtract_formula", "")),
        subtract_dataset=str(entry.get("subtract_dataset", "")),
        scale_group=str(entry.get("scale_group", "")),
        rows=tuple(str(r) for r in entry.get("rows", [])),  # pyright: ignore[reportGeneralTypeIssues]
        factor=float(entry.get("factor", 1.0)),  # pyright: ignore[reportArgumentType]
    )


def load_nonlinear_config(name: str | Path) -> NonlinearConfig:
    """Load a hyperpolarizability fit configuration.

    Parameters
    ----------
    name : str or pathlib.Path
        The gas (``data/gas/nonlinear/<gas>.toml``) or a path.

    Returns
    -------
    NonlinearConfig
        The configuration.

    Raises
    ------
    DataError
        If the configuration is malformed.
    """
    path = Path(name)
    if path.suffix != ".toml":
        path = NONLINEAR_DIRECTORY / f"{name}.toml"
    raw = path.read_bytes()
    record = tomllib.loads(raw.decode("utf-8"))
    allowed = {
        "gas",
        "description",
        "range_nm",
        "b_prime_prior",
        "a_prior",
        "b_prior",
        "fitted",
        "validators",
    }
    if set(record) - allowed:
        msg = f"{path.name}: unknown keys {sorted(set(record) - allowed)}"
        raise DataError(msg)
    low, high = (float(v) for v in record["range_nm"])
    return NonlinearConfig(
        gas=as_gas(record["gas"]),
        description=str(record.get("description", "")),
        range_nm=(low, high),
        b_prime_prior=_pair(record.get("b_prime_prior")),
        a_prior=_pair(record.get("a_prior")),
        b_prior=_pair(record.get("b_prior")),
        fitted=tuple(_options(e, path.name) for e in record.get("fitted", [])),
        validators=tuple(_options(e, path.name) for e in record.get("validators", [])),
        digest=hashlib.sha256(raw).hexdigest(),
    )


def _dispersion(
    parameters: Sequence[float] | np.ndarray, w2: np.ndarray, w4: np.ndarray
) -> np.ndarray:
    """Evaluate ``gamma0 (1 + A W2 + B W2^2 + B' W4)``.

    Parameters
    ----------
    parameters : sequence of float
        The gamma0 (au), A (1e-10 cm^2), B and B' (1e-20 cm^4).
    w2, w4 : numpy.ndarray
        The frequency sums, cm^-2 and cm^-4.

    Returns
    -------
    numpy.ndarray
        The gamma_T, atomic units.
    """
    gamma0, a, b, b_prime = parameters
    x2 = w2 * 1e-10
    return np.asarray(gamma0 * (1 + a * x2 + b * x2**2 + b_prime * w4 * 1e-20))


def _gradient(
    parameters: Sequence[float] | np.ndarray, w2: np.ndarray, w4: np.ndarray
) -> np.ndarray:
    """Return the derivatives of the dispersion with respect to its four parameters.

    Parameters
    ----------
    parameters : sequence of float
        The gamma0, A, B, B'.
    w2, w4 : numpy.ndarray
        The frequency sums, cm^-2 and cm^-4.

    Returns
    -------
    numpy.ndarray
        Shape ``(4, n)``.
    """
    gamma0, a, b, b_prime = parameters
    x2 = w2 * 1e-10
    x4 = w4 * 1e-20
    return np.vstack(
        [
            np.atleast_1d(1 + a * x2 + b * x2**2 + b_prime * x4),
            np.atleast_1d(gamma0 * x2),
            np.atleast_1d(gamma0 * x2**2),
            np.atleast_1d(gamma0 * x4),
        ]
    )


@dataclass(frozen=True)
class NonlinearFit:
    """A fitted hyperpolarizability model with its statistics.

    Attributes
    ----------
    config : NonlinearConfig
        The configuration.
    parameters : numpy.ndarray
        gamma0 (au), A (1e-10 cm^2), B, B' (1e-20 cm^4).
    covariance : numpy.ndarray
        Their covariance (4 x 4; B' row and column zero when fixed).
    scales : dict of str to float
        The fitted dataset scales.
    chi2 : float
        Chi-squared of the data.
    dof : int
        Data points minus parameters.
    dataset_chi2 : dict of str to tuple
        ``(n, chi2)`` per dataset.
    """

    config: NonlinearConfig
    parameters: Float[np.ndarray, " 4"]
    covariance: Float[np.ndarray, "4 4"] = field(repr=False)
    scales: Mapping[str, float] = field(default_factory=dict)
    chi2: float = 0.0
    dof: int = 0
    dataset_chi2: Mapping[str, tuple[int, float]] = field(default_factory=dict)


def _observations(options: Sequence[NonlinearOptions]) -> list[NonlinearObservation]:
    """Convert a list of datasets.

    Parameters
    ----------
    options : sequence of NonlinearOptions
        The datasets and how they enter.

    Returns
    -------
    list of NonlinearObservation
        The observations.
    """
    by_key = _catalogue.by_key()
    return [gamma_observation(by_key[o.key], o) for o in options]


def _scale_groups(
    data: Sequence[NonlinearObservation],
) -> tuple[list[str], list[int], dict[int, float]]:
    """Assign the datasets with a fitted scale to scale groups.

    Parameters
    ----------
    data : sequence of NonlinearObservation
        The fitted observations.

    Returns
    -------
    tuple
        The group names, each dataset's group index (-1 for no fitted scale) and the
        prior sigma of each group with a systematic prior.
    """
    groups: list[str] = []
    membership = []
    for obs in data:
        options = obs.options
        if not (options.free_scale or options.systematic):
            membership.append(-1)
            continue
        name = options.scale_group or f"dataset:{obs.key}"
        if name not in groups:
            groups.append(name)
        membership.append(groups.index(name))
    group_prior = {}
    for obs, index in zip(data, membership, strict=True):
        if index >= 0 and obs.options.systematic and not obs.options.free_scale:
            group_prior[index] = obs.options.systematic
    return groups, membership, group_prior


def fit_hyperpolarizability(config: NonlinearConfig) -> NonlinearFit:
    """Fit the dispersion model to a configuration's data.

    Parameters
    ----------
    config : NonlinearConfig
        The fit.

    Returns
    -------
    NonlinearFit
        The parameters, covariance and statistics.

    Raises
    ------
    DataError
        If a fitted value has no uncertainty.
    """
    data = _observations(config.fitted)
    for obs in data:
        if np.any(~(obs.sigma > 0)):
            msg = f"{obs.key}: a fitted value needs an uncertainty (extra_relative)"
            raise DataError(msg)
    groups, membership, group_prior = _scale_groups(data)
    fit_b_prime = config.b_prime_prior is not None
    offset = 4 if fit_b_prime else 3

    def unpack(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return the four model parameters and each dataset's scale.

        Parameters
        ----------
        theta : numpy.ndarray
            The fitted parameters.

        Returns
        -------
        tuple of numpy.ndarray
            The gamma0, A, B, B', and the scales.
        """
        physical = np.array(
            [theta[0], theta[1], theta[2], theta[3] if fit_b_prime else 0.0]
        )
        group_scales = 1 + theta[offset:]
        scales = np.array(
            [1.0 if index < 0 else group_scales[index] for index in membership]
        )
        return physical, scales

    def residuals(theta: np.ndarray) -> np.ndarray:
        """Return the weighted data residuals, then the prior residuals.

        Parameters
        ----------
        theta : numpy.ndarray
            The fitted parameters.

        Returns
        -------
        numpy.ndarray
            The residuals.
        """
        physical, scales = unpack(theta)
        parts = [
            (scale * _dispersion(physical, obs.w2, obs.w4) - obs.gamma) / obs.sigma
            for scale, obs in zip(scales, data, strict=True)
        ]
        if fit_b_prime:
            assert config.b_prime_prior is not None
            value, spread = config.b_prime_prior
            parts.append(np.array([(physical[3] - value) / spread]))
        for position, prior in ((1, config.a_prior), (2, config.b_prior)):
            if prior is not None:
                parts.append(np.array([(physical[position] - prior[0]) / prior[1]]))
        for index, systematic in group_prior.items():
            parts.append(np.array([theta[offset + index] / systematic]))
        return np.concatenate(parts)

    start_gamma = float(np.median(np.concatenate([o.gamma for o in data])))
    initial = (
        [start_gamma, 1.0, 2.0] + ([0.0] if fit_b_prime else []) + [0.0] * len(groups)
    )
    from scipy.optimize import least_squares  # noqa: PLC0415 - needed only to fit

    solution = least_squares(residuals, np.array(initial), method="trf", x_scale="jac")
    physical, scales = unpack(solution.x)
    jacobian = np.asarray(solution.jac, dtype=float)
    covariance_theta = np.linalg.pinv(jacobian.T @ jacobian)
    n_physical = 4 if fit_b_prime else 3
    covariance = np.zeros((4, 4))
    covariance[:n_physical, :n_physical] = covariance_theta[:n_physical, :n_physical]
    dataset_chi2 = {}
    chi2 = 0.0
    for scale, obs in zip(scales, data, strict=True):
        value = float(
            np.sum(
                (
                    (scale * _dispersion(physical, obs.w2, obs.w4) - obs.gamma)
                    / obs.sigma
                )
                ** 2
            )
        )
        dataset_chi2[obs.key] = (len(obs.gamma), value)
        chi2 += value
    n_points = sum(len(o.gamma) for o in data)
    return NonlinearFit(
        config=config,
        parameters=physical,
        covariance=covariance,
        scales={obs.key: float(s) for obs, s in zip(data, scales, strict=True)},
        chi2=chi2,
        dof=n_points - n_physical - len(groups),
        dataset_chi2=dataset_chi2,
    )


def validation(result: NonlinearFit) -> dict[str, dict[str, float | int | str]]:
    """Compare a fit with its validator datasets.

    Parameters
    ----------
    result : NonlinearFit
        The fit.

    Returns
    -------
    dict
        Per validator: the number of values, the mean ratio of data to model and the
        chi-squared where uncertainties are stated; or a note.
    """
    report: dict[str, dict[str, float | int | str]] = {}
    by_key = _catalogue.by_key()
    for options in result.config.validators:
        try:
            obs = gamma_observation(by_key[options.key], options)
        except DataError as error:
            report[options.key] = {"note": str(error)}
            continue
        model = _dispersion(result.parameters, obs.w2, obs.w4)
        stated = obs.sigma > 0
        report[options.key] = {
            "n": len(obs.gamma),
            "mean_ratio": float(np.mean(obs.gamma / model)),
            "chi2": float(
                np.sum(((obs.gamma - model)[stated] / obs.sigma[stated]) ** 2)
            ),
        }
    return report


def dumps(result: NonlinearFit) -> str:
    """Serialise a fit as JSON text.

    Parameters
    ----------
    result : NonlinearFit
        The fit.

    Returns
    -------
    str
        Indented JSON.
    """
    data = _observations(result.config.fitted)
    record = {
        "gas": result.config.gas.value,
        "config_sha256": result.config.digest,
        "range_nm": list(result.config.range_nm),
        "parameters": {
            "gamma0_au": float(result.parameters[0]),
            "A_1e-10_cm2": float(result.parameters[1]),
            "B_1e-20_cm4": float(result.parameters[2]),
            "B_prime_1e-20_cm4": float(result.parameters[3]),
        },
        "covariance": result.covariance.tolist(),
        "max_w2_fitted_cm-2": float(max(float(np.max(o.w2)) for o in data)),
        "scales": dict(result.scales),
        "chi2": result.chi2,
        "dof": result.dof,
        "dataset_chi2": {k: list(v) for k, v in result.dataset_chi2.items()},
        "validation": validation(result),
    }
    return json.dumps(record, indent=1) + "\n"


SUPPORT_DERIVATION_VERSION: Final = 1
"""The version of the support records' derivation (OEP 5), stored in each record and
raised whenever the derivation changes what it writes. It is separate from the models'
:data:`EVALUATION_VERSION`: a record describes a model and changes no evaluation.
"""

SUPPORT_LIMIT_KINDS: Final = (
    "fitted_extent",
    "level",
    "band_origin",
    "lower_bound",
    "upper_bound",
    "absorption_evidence",
)
"""What a :class:`SupportLimit` is (OEP 5, section 6).

- ``fitted_extent``: the largest of a quantity over the fitted points;
- ``level``: a spectroscopic level or threshold energy, as its source gives it;
- ``band_origin``: a molecular state's vibrationless origin, a sum of 0-0 bands;
- ``lower_bound``, ``upper_bound``: an experimental bound, in wavenumber;
- ``absorption_evidence``: absorption detected, or not, over a finite sampled window;
  not a level and not a bound on one.
"""

_EDGE_KINDS: Final = frozenset({"level", "band_origin"})
_LIMIT_NAMES: Final = ("e1", "e2", "e2_any")
_SUPPORT_DIRECTORY: Final = GAS_DATA_DIRECTORY / "models" / "nonlinear_support"
_STANDING: Final = frozenset({"status", "status_reason", "modified", "base"})


@dataclass(frozen=True)
class SupportPoint:
    r"""One fitted value as a quartet: where it comes from and its four signed legs.

    Attributes
    ----------
    dataset : str
        ``source/dataset``.
    row : str
        The dataset's row.
    process : str
        The process, as :func:`frequency_sums` names it.
    wavelength : float
        The fundamental (pump) vacuum wavelength, m; ``inf`` for a static value.
    second_wavelength : float or None
        The Stokes (``cars``) or probe (``pump_probe``) vacuum wavelength, m.
    legs : tuple of float
        The signed wavenumbers :math:`(-\nu_\sigma, \nu_1, \nu_2, \nu_3)`, m^-1:
        inverse vacuum wavelengths, the output with a minus sign, summing to zero.
    """

    dataset: str
    row: str
    process: str
    wavelength: float
    second_wavelength: float | None
    legs: tuple[float, float, float, float]

    def __post_init__(self) -> None:
        """Hold the legs as a tuple of floats."""
        object.__setattr__(self, "legs", tuple(float(x) for x in self.legs))


@dataclass(frozen=True)
class SupportUncertainty:
    """A limit's uncertainty as its source states it, or the fact that it states none.

    Printed precision is not an uncertainty, and an unstated one is not zero.

    Attributes
    ----------
    stated : bool
        Whether the source states an uncertainty.
    value : float or None
        The stated uncertainty, m^-1; None where none is stated.
    level : str
        Its confidence level as the source gives it (``"not stated"`` where it does
        not).
    provenance : str
        Where it comes from, or why there is none.
    """

    stated: bool
    value: float | None
    level: str
    provenance: str


@dataclass(frozen=True)
class SupportCandidate:
    """One candidate of a limit taken as the smallest of several.

    Attributes
    ----------
    datasets : tuple of str
        Its datasets: one level, or the chain of band origins whose sum it is.
    wavenumber : float
        Its energy over hc, m^-1.
    """

    datasets: tuple[str, ...]
    wavenumber: float


@dataclass(frozen=True)
class AbsorptionEvidence:
    """What a finite absorption measurement shows: not a level, and not a bound on one.

    Attributes
    ----------
    detected : bool
        Whether absorption was detected.
    sampled : tuple of float
        The lowest and highest wavenumber of the curated values, m^-1.
    qualification : str
        The source's own qualification of the detection.
    """

    detected: bool
    sampled: tuple[float, float]
    qualification: str


@dataclass(frozen=True)
class SupportLimit:
    """One fact of a support record: a wavenumber, what kind of fact, and its sources.

    Attributes
    ----------
    wavenumber : float
        The limit, m^-1 (an inverse vacuum wavelength; an energy over hc).
    kind : str
        One of :data:`SUPPORT_LIMIT_KINDS`. A bound's direction is in wavenumber;
        inverting to a wavelength reverses it.
    datasets : tuple of str
        The datasets it is taken from: for a fitted extent, those of its witnesses; for
        a smallest of candidates, the winner's.
    location : str
        Where in the datasets: the rows, columns or parameters read.
    derivation : str
        How the wavenumber follows from them.
    registers : tuple of str
        The register entries that record the choice (``TC-n``).
    witnesses : tuple of SupportPoint
        For a fitted extent, every fitted point that attains it, in the order of the
        fit's datasets and their rows; empty otherwise.
    candidates : tuple of SupportCandidate
        For a smallest of several, every candidate; empty otherwise.
    uncertainty : SupportUncertainty or None
        The source's uncertainty, or that it states none; None where none applies (a
        fitted extent, absorption evidence).
    evidence : AbsorptionEvidence or None
        For absorption evidence, what was measured; None otherwise.
    qualifications : tuple of str
        The approximations and selection assumptions behind it.

    Raises
    ------
    DataError
        For an unknown kind, or absorption evidence without its measurement.
    """

    wavenumber: float
    kind: str
    datasets: tuple[str, ...]
    location: str
    derivation: str
    registers: tuple[str, ...]
    witnesses: tuple[SupportPoint, ...] = ()
    candidates: tuple[SupportCandidate, ...] = ()
    uncertainty: SupportUncertainty | None = None
    evidence: AbsorptionEvidence | None = None
    qualifications: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Check the kind and hold every sequence as a tuple."""
        if self.kind not in SUPPORT_LIMIT_KINDS:
            msg = f"unknown support limit kind {self.kind!r}"
            raise DataError(msg)
        if (self.kind == "absorption_evidence") != (self.evidence is not None):
            msg = "absorption evidence, and only it, carries its measurement"
            raise DataError(msg)
        for name in ("datasets", "registers", "witnesses", "candidates"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(self, "qualifications", tuple(self.qualifications))


@dataclass(frozen=True)
class HyperpolarizabilitySupport(Identified):
    r"""Where a gas's hyperpolarizability fit is supported: facts, with their sources.

    The record of OEP 5, read from ``data/gas/models/nonlinear_support/<gas>.json``
    through :attr:`HyperpolarizabilityModel.support`. Wavenumbers are signed inverse
    vacuum wavelengths in m^-1 (:math:`\omega = 2 \pi c \nu`); a quartet's legs are
    :math:`(-\nu_\sigma, \nu_1, \nu_2, \nu_3)`, summing to zero, and
    :math:`W_2` is the sum of their squares, m^-2. It states facts and carries no rule:
    how a propagation code uses them is its own decision (TC-868). It changes no
    evaluation, and is not part of the model's digest; it has its own.

    The resonance limits hold in the electric-dipole approximation, from the ground
    state, with the selection assumptions each limit records (TC-865). They say nothing
    of rotational or vibrational Raman resonances, linewidths, thermally populated
    states or weak forbidden transitions, and being below them establishes no accuracy.

    Attributes
    ----------
    gas : Gas
        The species.
    model_digest : str
        The digest of the model whose fitted data it describes.
    derivation_version : int
        :data:`SUPPORT_DERIVATION_VERSION` when it was derived.
    fitted_points : int
        The number of fitted points, after the fit's windows and row selections.
    max_w2 : float
        The largest fitted :math:`W_2`, m^-2.
    max_w2_witnesses : tuple of SupportPoint
        The fitted points that attain it.
    leg_max : SupportLimit
        The largest single leg :math:`|\nu_k|` of a fitted point, m^-1.
    pair_max : SupportLimit
        The largest pair sum :math:`|\nu_a + \nu_b|`, :math:`a \ne b`, of a fitted
        point, m^-1.
    e1 : SupportLimit or None
        The lowest included level dipole-coupled to the ground state, which no leg (the
        output included) is to reach; or absorption evidence where no level is curated.
    e2 : SupportLimit or None
        The lowest included level two dipole steps from the ground state, which no pair
        sum is to reach, for one linear polarisation common to all legs.
    e2_any : SupportLimit or None
        The same for any polarisation: the smallest of all eligible levels.
    ionisation : SupportLimit
        The ionisation energy over hc, m^-1.
    missing_reasons : Mapping of str to str
        Why each of ``e1``, ``e2`` and ``e2_any`` that is None is missing.
    input_digests : Mapping of str to str
        SHA-256 digests of what it was derived from: the selection, the fit's
        configuration and each selected dataset.

    Raises
    ------
    DataError
        If a resonance limit is both given and missing, or neither.

    Examples
    --------
    >>> from optmat.gas.gases import Gas
    >>> support = hyperpolarizability_models()[Gas.Ar].support
    >>> round(1e9 / support.leg_max.wavenumber, 1), support.leg_max.witnesses[0].process
    (229.0, 'eshg')
    >>> support.leg_max.datasets
    ('shelton_1990/table3_Ar_literature',)
    >>> round(1e9 / support.leg_max_in_w2, 1)  # derived, sqrt(3 W2,max / 4), nm
    215.9
    >>> round(1e9 / support.e1.wavenumber, 1), support.e1.datasets
    (106.7, ('nist_asd_2026/level_Ar_odd_J1',))
    >>> round(1e9 / support.band_edge(), 1), round(1e9 / support.band_edge("any"), 1)
    (189.4, 192.1)
    >>> hyperpolarizability_models()[Gas.SF6].support.band_edge() is None
    True
    """

    gas: Gas
    model_digest: str
    derivation_version: int
    fitted_points: int
    max_w2: float
    max_w2_witnesses: tuple[SupportPoint, ...]
    leg_max: SupportLimit
    pair_max: SupportLimit
    e1: SupportLimit | None
    e2: SupportLimit | None
    e2_any: SupportLimit | None
    ionisation: SupportLimit
    missing_reasons: Mapping[str, str] = field(default_factory=dict)
    input_digests: Mapping[str, str] = field(default_factory=dict, repr=False)

    IDENTITY_KIND: ClassVar[str] = "gas_hyperpolarizability_support"
    EVALUATION_VERSION: ClassVar[int] = 1

    def __post_init__(self) -> None:
        """Take ownership, and check that every absent limit, and only it, has a reason.

        Raises
        ------
        DataError
            If a limit is both given and missing, or neither.
        """
        super().__post_init__()
        unknown = set(self.missing_reasons) - set(_LIMIT_NAMES)
        if unknown:
            msg = f"{self.gas.value}: no limit {sorted(unknown)} to be missing"
            raise DataError(msg)
        for name in _LIMIT_NAMES:
            if (getattr(self, name) is None) != (name in self.missing_reasons):
                msg = f"{self.gas.value}: {name} must be given or have a missing reason"
                raise DataError(msg)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the record's name, gas and the sources of every dataset it cites.

        Returns
        -------
        tuple
            ``("hyperpolarizability_support/<gas>", gas, sources)``.
        """
        limits = [self.leg_max, self.pair_max, self.e1, self.e2, self.e2_any]
        keys = [
            k for limit in [*limits, self.ionisation] if limit for k in limit.datasets
        ]
        keys += [p.dataset for p in self.max_w2_witnesses]
        return (
            f"hyperpolarizability_support/{self.gas.value}",
            self.gas.value,
            source_ids(*(f"{key.partition('#')[0]}/" for key in keys)),
        )

    @property
    def digest(self) -> str:
        """The record's own digest, ``identity.digest`` (OEP 1's machinery).

        Returns
        -------
        str
            64 hexadecimal digits.
        """
        return self.identity.digest

    @property
    def leg_max_in_w2(self) -> float:
        r"""Derived: the largest leg of any quartet within the fitted W2, m^-1.

        Each leg is minus the sum of the other three, so by the Cauchy-Schwarz
        inequality :math:`\nu_k^2 \le 3 W_2 / 4`, with equality for the third
        harmonic's output: :math:`\sqrt{3 W_{2,\max} / 4}`.

        Returns
        -------
        float
            The wavenumber, m^-1.
        """
        return math.sqrt(0.75 * self.max_w2)

    def band_edge(self, polarisation: str = "linear") -> float | None:
        r"""Return the derived :math:`\min(E_1, E_2/2)`, a sufficient bound on all legs.

        If every leg of a quartet, generated outputs included, is below it, no leg
        reaches :math:`E_1` and no pair sum reaches :math:`E_2`. Bounding a field's
        input band alone does not bound the outputs it generates. Exceeding it does not
        prove a quartet resonant. For one linear polarisation the self-phase-modulation
        quartet reaches :math:`E_2` at the edge; for any polarisation it need not (an
        atom's even :math:`J = 1` level is not reached by two photons of equal
        frequency). The number carries no accuracy margin; its inputs' uncertainties
        and qualifications are on the limits.

        Parameters
        ----------
        polarisation : str, optional
            ``"linear"``, one linear polarisation common to all legs (``e2``), or
            ``"any"`` (``e2_any``).

        Returns
        -------
        float or None
            The bound, m^-1; None where either limit is missing or is not a level or a
            band origin (absorption evidence, a bound).

        Raises
        ------
        DataError
            For any other polarisation.
        """
        pairs = {"linear": self.e2, "any": self.e2_any}
        if polarisation not in pairs:
            msg = f"polarisation must be 'linear' or 'any', not {polarisation!r}"
            raise DataError(msg)
        e1, e2 = self.e1, pairs[polarisation]
        if e1 is None or e2 is None or {e1.kind, e2.kind} - _EDGE_KINDS:
            return None
        return min(e1.wavenumber, e2.wavenumber / 2)


def _plain(value: object) -> Any:
    """Return a record's value as plain JSON values (an infinite wavelength as null).

    Parameters
    ----------
    value : object
        A field's value.

    Returns
    -------
    object
        Dicts, lists, strings, numbers, booleans or None.
    """
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            f.name: _plain(getattr(value, f.name))
            for f in dataclasses.fields(value)
            if f.name not in _STANDING
        }
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, tuple):
        return [_plain(v) for v in value]
    if isinstance(value, Gas):
        return value.value
    if isinstance(value, float) and math.isinf(value):
        return None
    return value


def _support_dumps(record: HyperpolarizabilitySupport) -> str:
    """Serialise a support record as JSON text.

    Parameters
    ----------
    record : HyperpolarizabilitySupport
        The record.

    Returns
    -------
    str
        Indented JSON, as stored in ``data/gas/models/nonlinear_support/<gas>.json``.
    """
    return json.dumps(_plain(record), indent=1, allow_nan=False) + "\n"


def _fields(cls: type, raw: Any, where: str) -> Mapping[str, Any]:
    """Return a stored object's fields, checking they are exactly the class's.

    Parameters
    ----------
    cls : type
        The dataclass.
    raw : object
        The JSON object.
    where : str
        For error messages.

    Returns
    -------
    Mapping
        The fields.

    Raises
    ------
    DataError
        If it is not an object with exactly the class's fields.
    """
    names = {f.name for f in dataclasses.fields(cls)} - _STANDING
    if not isinstance(raw, Mapping) or set(raw) != names:
        msg = f"{where}: a {cls.__name__} needs exactly the fields {sorted(names)}"
        raise DataError(msg)
    return raw


def _point(raw: Any, where: str) -> SupportPoint:
    """Read a stored :class:`SupportPoint`.

    Parameters
    ----------
    raw : object
        The JSON object.
    where : str
        For error messages.

    Returns
    -------
    SupportPoint
        The point.
    """
    d = _fields(SupportPoint, raw, where)
    second = d["second_wavelength"]
    return SupportPoint(
        dataset=str(d["dataset"]),
        row=str(d["row"]),
        process=str(d["process"]),
        wavelength=math.inf if d["wavelength"] is None else float(d["wavelength"]),
        second_wavelength=None if second is None else float(second),
        legs=tuple(float(x) for x in d["legs"]),  # pyright: ignore[reportArgumentType]
    )


def _limit(raw: Any, where: str) -> SupportLimit | None:
    """Read a stored :class:`SupportLimit`, or None.

    Parameters
    ----------
    raw : object
        The JSON object, or None.
    where : str
        For error messages.

    Returns
    -------
    SupportLimit or None
        The limit.
    """
    if raw is None:
        return None
    d = _fields(SupportLimit, raw, where)
    uncertainty, evidence = d["uncertainty"], d["evidence"]
    if uncertainty is not None:
        u = _fields(SupportUncertainty, uncertainty, where)
        uncertainty = SupportUncertainty(
            stated=bool(u["stated"]),
            value=None if u["value"] is None else float(u["value"]),
            level=str(u["level"]),
            provenance=str(u["provenance"]),
        )
    if evidence is not None:
        e = _fields(AbsorptionEvidence, evidence, where)
        low, high = (float(x) for x in e["sampled"])
        evidence = AbsorptionEvidence(
            detected=bool(e["detected"]),
            sampled=(low, high),
            qualification=str(e["qualification"]),
        )
    return SupportLimit(
        wavenumber=float(d["wavenumber"]),
        kind=str(d["kind"]),
        datasets=tuple(str(k) for k in d["datasets"]),
        location=str(d["location"]),
        derivation=str(d["derivation"]),
        registers=tuple(str(r) for r in d["registers"]),
        witnesses=tuple(_point(p, where) for p in d["witnesses"]),
        candidates=tuple(
            SupportCandidate(
                datasets=tuple(str(k) for k in c["datasets"]),
                wavenumber=float(c["wavenumber"]),
            )
            for c in (_fields(SupportCandidate, c, where) for c in d["candidates"])
        ),
        uncertainty=uncertainty,
        evidence=evidence,
        qualifications=tuple(str(q) for q in d["qualifications"]),
    )


def _support_loads(
    text: str, where: str = "support record"
) -> HyperpolarizabilitySupport:
    """Read a support record from its JSON text.

    Parameters
    ----------
    text : str
        The JSON, as :func:`_support_dumps` writes it.
    where : str, optional
        For error messages.

    Returns
    -------
    HyperpolarizabilitySupport
        The record, status ``validated``.

    Raises
    ------
    DataError
        If the text is not a record of the current derivation version.
    """
    d = _fields(HyperpolarizabilitySupport, json.loads(text), where)
    if d["derivation_version"] != SUPPORT_DERIVATION_VERSION:
        msg = (
            f"{where}: derivation version {d['derivation_version']}, not "
            f"{SUPPORT_DERIVATION_VERSION}; re-derive it"
        )
        raise DataError(msg)
    leg_max, pair_max = _limit(d["leg_max"], where), _limit(d["pair_max"], where)
    ionisation = _limit(d["ionisation"], where)
    if leg_max is None or pair_max is None or ionisation is None:
        msg = f"{where}: leg_max, pair_max and ionisation are always given"
        raise DataError(msg)
    return HyperpolarizabilitySupport(
        gas=as_gas(d["gas"]),
        model_digest=str(d["model_digest"]),
        derivation_version=int(d["derivation_version"]),
        fitted_points=int(d["fitted_points"]),
        max_w2=float(d["max_w2"]),
        max_w2_witnesses=tuple(_point(p, where) for p in d["max_w2_witnesses"]),
        leg_max=leg_max,
        pair_max=pair_max,
        e1=_limit(d["e1"], where),
        e2=_limit(d["e2"], where),
        e2_any=_limit(d["e2_any"], where),
        ionisation=ionisation,
        missing_reasons={str(k): str(v) for k, v in d["missing_reasons"].items()},
        input_digests={str(k): str(v) for k, v in d["input_digests"].items()},
        status="validated",
    )


@functools.cache
def _stored_supports() -> dict[Gas, HyperpolarizabilitySupport]:
    """Return every stored support record, by gas.

    Returns
    -------
    dict of Gas to HyperpolarizabilitySupport
        The records of ``data/gas/models/nonlinear_support/<gas>.json``.

    Raises
    ------
    DataError
        If a record is malformed or stored under another gas's name.
    """
    records = {}
    for path in sorted(_SUPPORT_DIRECTORY.glob("*.json")):
        record = _support_loads(path.read_text(encoding="utf-8"), path.name)
        if record.gas.value != path.stem:
            msg = f"{path.name} holds the support record of {record.gas.value}"
            raise DataError(msg)
        records[record.gas] = record
    return records


@dataclass(frozen=True)
class HyperpolarizabilityModel(Identified):
    """A gas's fitted electronic second hyperpolarizability, for any process.

    Attributes
    ----------
    gas : Gas
        The species.
    parameters : tuple of float
        gamma0 (au), A (1e-10 cm^2), B, B' (1e-20 cm^4).
    covariance : numpy.ndarray
        Their covariance.
    reduced_chi2 : float
        Of the fit.
    wavelength_range : tuple of float
        The fundamental vacuum-wavelength range over which the model is valid, m.
    record : dict
        The stored record.
    """

    gas: Gas
    parameters: tuple[float, float, float, float]
    covariance: Float[np.ndarray, "4 4"] = field(repr=False)
    reduced_chi2: float = 0.0
    wavelength_range: tuple[float, float] = (0.0, math.inf)
    record: Mapping[str, object] = field(default_factory=dict, repr=False)

    IDENTITY_KIND: ClassVar[str] = "gas_hyperpolarizability"
    EVALUATION_VERSION: ClassVar[int] = 1
    EVALUATES_THROUGH: ClassVar[tuple[str, ...]] = ("optmat.gas.nonlinear",)

    def _identity_label(self) -> tuple[str, str, tuple[str, ...]]:
        """Return the model's name, gas and the sources of its fitted datasets.

        Returns
        -------
        tuple
            ``("hyperpolarizability/<gas>", gas, sources)``.
        """
        fitted = self.record.get("dataset_chi2", {})
        assert isinstance(fitted, Mapping)
        sources = tuple(sorted({str(key).split("/", 1)[0] for key in fitted}))
        return f"hyperpolarizability/{self.gas.value}", self.gas.value, sources

    @property
    def support(self) -> HyperpolarizabilitySupport | None:
        """Where this fit is supported, and the gas's resonance limits (OEP 5).

        The stored record whose ``model_digest`` is this model's digest. It is not part
        of the model: it changes neither the digest nor any evaluated gamma.

        Returns
        -------
        HyperpolarizabilitySupport or None
            The record; None for a model built from a caller's parameters (status
            ``user``) that no curated record describes.

        Raises
        ------
        DataError
            If a curated model has no matching record: the stored records are missing
            or stale (re-derive them with ``tools/fit/nonlinear_support_record.py``).
        """
        record = _stored_supports().get(self.gas)
        if record is not None and record.model_digest == self.identity.digest:
            return record
        if self.status == "user":
            return None
        msg = (
            f"hyperpolarizability/{self.gas.value}: no stored support record describes "
            "this curated model; re-derive it with "
            "tools/fit/nonlinear_support_record.py"
        )
        raise DataError(msg)

    def _sums(
        self, process: str, wavelength: np.ndarray, other: np.ndarray | None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return W2, W4 at vacuum wavelengths (m), warning outside the model's range.

        The model warns (:class:`~optmat.exceptions.ModelRangeWarning`) at a
        fundamental wavelength outside its declared range, and at a W2 more than
        :data:`W2_MARGIN` times the largest fitted W2, where the fourth-order form is
        an extrapolation (a third harmonic of a short fundamental, say).

        Parameters
        ----------
        process : str
            The process.
        wavelength : numpy.ndarray
            The fundamental (pump) wavelength, m; ``inf`` is static.
        other : numpy.ndarray or None
            The second wavelength (Stokes or probe), m.

        Returns
        -------
        tuple of numpy.ndarray
            W2, W4.
        """
        low, high = self.wavelength_range
        finite = wavelength[np.isfinite(wavelength)]
        # The declared ends themselves are inside (400.0 * 1e-9 rounds above 400e-9).
        if np.any(finite < low * (1 - 1e-12)) or np.any(finite > high * (1 + 1e-12)):
            msg = (
                f"hyperpolarizability of {self.gas.value} evaluated outside "
                f"{low * 1e9:.4g}-{high * 1e9:.4g} nm"
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)
        nu = np.where(np.isfinite(wavelength), 1e-2 / wavelength, 0.0)
        second = (
            None if other is None else np.asarray(1e-2 / np.asarray(other, dtype=float))
        )
        w2, w4 = frequency_sums(process, nu, second)
        fitted = self.record.get("max_w2_fitted_cm-2")
        if fitted is not None and np.any(w2 > W2_MARGIN * float(fitted)):  # pyright: ignore[reportArgumentType]
            msg = (
                f"{process} hyperpolarizability of {self.gas.value} evaluated at W2 = "
                f"{float(np.max(w2)):.3g} cm^-2, beyond the fitted data "
                f"({float(fitted):.3g} cm^-2): the dispersion form is extrapolated"  # pyright: ignore[reportArgumentType]
            )
            warnings.warn(msg, ModelRangeWarning, stacklevel=3)
        return w2, w4

    def gamma_au(
        self,
        process: str,
        wavelength: Float[np.ndarray, "*shape"] | float,
        other: Float[np.ndarray, "*shape"] | float | None = None,
    ) -> Float[np.ndarray, "*shape"]:
        """Return gamma_T for a process at a fundamental wavelength.

        Parameters
        ----------
        process : str
            ``static``, ``kerr``, ``eshg``, ``dc_kerr``, ``thg``, ``cars`` or
            ``pump_probe``.
        wavelength : array_like or float
            Fundamental (pump) vacuum wavelength, m.
        other : array_like or float, optional
            The Stokes (``cars``) or probe (``pump_probe``) wavelength, m.

        Returns
        -------
        numpy.ndarray
            The gamma_T, atomic units.
        """
        wavelength = np.asarray(wavelength, dtype=float)
        w2, w4 = self._sums(
            process, wavelength, None if other is None else np.asarray(other)
        )
        return _dispersion(self.parameters, w2, w4)

    def gamma_sigma_au(
        self,
        process: str,
        wavelength: Float[np.ndarray, "*shape"] | float,
        other: Float[np.ndarray, "*shape"] | float | None = None,
        *,
        birge: bool = True,
    ) -> Float[np.ndarray, "*shape"]:
        """Return the one-sigma uncertainty of :meth:`gamma_au` from the covariance.

        Parameters
        ----------
        process : str
            The process.
        wavelength : array_like or float
            Fundamental vacuum wavelength, m.
        other : array_like or float, optional
            The second wavelength, m.
        birge : bool, optional
            Scale the covariance by the reduced chi-squared when it exceeds one.

        Returns
        -------
        numpy.ndarray
            Atomic units.
        """
        wavelength = np.atleast_1d(np.asarray(wavelength, dtype=float))
        w2, w4 = self._sums(
            process, wavelength, None if other is None else np.asarray(other)
        )
        gradient = _gradient(self.parameters, w2, w4)
        covariance = self.covariance
        if birge and self.reduced_chi2 > 1:
            covariance = covariance * self.reduced_chi2
        return np.sqrt(np.einsum("in,ij,jn->n", gradient, covariance, gradient))

    def gamma_si(
        self,
        process: str,
        wavelength: Float[np.ndarray, "*shape"] | float,
        other: Float[np.ndarray, "*shape"] | float | None = None,
    ) -> Float[np.ndarray, "*shape"]:
        """Return gamma_T in SI (C^4 m^4 J^-3).

        Parameters
        ----------
        process : str
            The process.
        wavelength : array_like or float
            Fundamental vacuum wavelength, m.
        other : array_like or float, optional
            The second wavelength, m.

        Returns
        -------
        numpy.ndarray
            The gamma_T, SI.
        """
        return np.asarray(
            self.gamma_au(process, wavelength, other) * AU_HYPERPOLARIZABILITY
        )

    def gamma3(
        self, process: str, wavelength: Float[np.ndarray, "*shape"] | float
    ) -> Float[np.ndarray, "*shape"]:
        """Return the per-molecule coefficient ``gamma3 = gamma_T / (6 eps0)``.

        The third-order polarisation of the real field is
        ``P = eps0 gamma3 N E^3`` (m^5/V^2 per molecule).

        Parameters
        ----------
        process : str
            ``kerr`` for SPM and XPM, ``thg`` for the third harmonic.
        wavelength : array_like or float
            Fundamental vacuum wavelength, m.

        Returns
        -------
        numpy.ndarray
            The coefficient, m^5/V^2.
        """
        return np.asarray(self.gamma_si(process, wavelength) / (6 * EPS0))

    def n2(
        self,
        wavelength: Float[np.ndarray, "*shape"] | float,
        state: State,
        *,
        ideal: bool = False,
    ) -> Float[np.ndarray, "*shape"]:
        """Return the self Kerr coefficient n2 (m^2/W) in the dilute limit.

        ``n2 = N gamma_Kerr / (8 eps0^2 c)`` for linearly polarised light.

        Parameters
        ----------
        wavelength : array_like or float
            Vacuum wavelength, m.
        state : State
            Pressure and temperature.
        ideal : bool, optional
            Use the ideal-gas density instead of the real-gas one (always, for a
            species without an equation of state).

        Returns
        -------
        numpy.ndarray
            The n2, m^2/W.
        """
        # A species without an equation of state (ozone) is taken as ideal.
        ideal = ideal or info(self.gas).fluid is None
        density = number_density(self.gas, state, ideal=ideal)
        return np.asarray(
            density * self.gamma_si("kerr", wavelength) / (8 * EPS0**2 * C)
        )


@functools.cache
def hyperpolarizability_models() -> dict[Gas, HyperpolarizabilityModel]:
    """Return every stored hyperpolarizability model, by gas.

    Returns
    -------
    dict of Gas to HyperpolarizabilityModel
        The models of ``data/gas/models/nonlinear/<gas>.json``.
    """
    models = {}
    for path in sorted(_MODEL_DIRECTORY.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        parameters = record["parameters"]
        low, high = (float(x) * 1e-9 for x in record["range_nm"])
        gas = as_gas(record["gas"])
        models[gas] = HyperpolarizabilityModel(
            gas=gas,
            parameters=(
                float(parameters["gamma0_au"]),
                float(parameters["A_1e-10_cm2"]),
                float(parameters["B_1e-20_cm4"]),
                float(parameters["B_prime_1e-20_cm4"]),
            ),
            covariance=np.array(record["covariance"]),
            reduced_chi2=record["chi2"] / max(record["dof"], 1),
            wavelength_range=(low, high),
            record=record,
            status="validated",
        )
    return models

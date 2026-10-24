# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Fitting a solid's universal dispersion model to the curated data.

A fit starts from a published UDM formula of a source record (its newAD2 parameter
file's values) and frees a chosen set of its parameters. It adjusts them, by least
squares, to curated n and k at their stated uncertainties:

- **n:** `refractive_index`, or a `refractive_index_deviation` from a named formula,
  which the fit adds back to that formula; an index relative to room air
  (``parameters.index_medium = "air"``) is taken to vacuum with the index of that air
  (TC-398);
- **k:** `extinction_coefficient`, compared in log k, since k spans many decades;
- **offsets:** each fitted n dataset may carry a free additive offset with a prior
  (its grade and temperature: silica grades differ by ~1e-5 in n, and dn/dT is ~1e-5
  per kelvin);
- **priors:** free parameters may carry priors.

The model stays causal, passive and sum-rule normalised whatever the parameters, so the
fitted dielectric function is analytic like the published one. The configuration is
``data/solid/fits/<name>.toml``; the result is stored as
``data/solid/models/fitted/<name>.json`` with its covariance (TC-322).
"""

import hashlib
import json
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import numpy as np
from jaxtyping import Float
from scipy import optimize

from optmat import _catalogue
from optmat.conditions import State, air_refractivity
from optmat.constants import E_CHARGE, HBAR, STD_ATMOSPHERE, C
from optmat.datasets import Dataset
from optmat.exceptions import DataError
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models.base import complex_index
from optmat.solid.models.formulas import published_models
from optmat.solid.udm import (
    UniversalDispersionModel,
    from_parameters,
    parse_attributes,
)

__all__ = [
    "FIT_DIRECTORY",
    "DatasetOptions",
    "FitConfig",
    "FitResult",
    "dumps",
    "fit",
    "load_config",
    "validation",
]

FIT_DIRECTORY: Final = SOLID_DATA_DIRECTORY / "fits"
"""The solid fit configurations."""


@dataclass(frozen=True)
class DatasetOptions:
    """How one curated dataset enters a fit.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    range_nm : tuple of float
        The vacuum-wavelength window used, nm.
    sigma : float or None
        Uncertainty used where the dataset states none (n: absolute; k: relative).
    extra : float
        Added in quadrature to every point's uncertainty (same units as ``sigma``).
    offset_sigma : float or None
        The prior (one sigma) of a free additive offset in n; None for no offset.
    offset_group : str
        Datasets naming the same group share one offset (a common reference level,
        such as whether a source's n is relative to air); empty for the dataset's own.
    reference : str
        For a `refractive_index_deviation`, the formula model it is relative to.
    thermo_optic : str
        A `thermo_optic_coefficient` dataset (``source/dataset``) that takes the
        values from their sample temperature to ``temperature_k``, interpolated in
        wavelength; empty for none.
    temperature_k : float or None
        The temperature the fit refers to, K, for ``thermo_optic``.
    """

    key: str
    range_nm: tuple[float, float] = (0.0, math.inf)
    sigma: float | None = None
    extra: float = 0.0
    offset_sigma: float | None = None
    offset_group: str = ""
    reference: str = ""
    thermo_optic: str = ""
    temperature_k: float | None = None


@dataclass(frozen=True)
class FitConfig:
    """A solid fit's configuration.

    Attributes
    ----------
    name : str
        The configuration's name (its file stem).
    material : str
        The material's id.
    description : str
        What the fit is and why.
    range_nm : tuple of float
        The range the fitted model supports, nm.
    start : str
        The starting UDM formula, ``<source>/<formula id>`` of a curated source whose
        formula holds a newAD2 parameter file's values as its terms.
    free : tuple of str
        The freed parameters, in newAD2's names.
    scales : dict of str to tuple of str
        Groups of parameters multiplied by one fitted factor each (the OH bands of a
        grade, say), by the group's name.
    added_attributes : tuple of str
        Terms added to the start model (``loc=1`` adds one Gaussian peak).
    added_values : dict of str to float
        The added terms' starting parameters.
    priors : dict of str to tuple of float
        Priors (value, sigma) on free parameters.
    fitted, validators : tuple of DatasetOptions
        The fitted datasets and the datasets only compared.
    digest : str
        The configuration file's SHA-256.
    """

    name: str
    material: str
    description: str
    range_nm: tuple[float, float]
    start: str
    free: tuple[str, ...]
    scales: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    added_attributes: tuple[str, ...] = ()
    added_values: Mapping[str, float] = field(default_factory=dict)
    priors: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    fitted: tuple[DatasetOptions, ...] = ()
    validators: tuple[DatasetOptions, ...] = ()
    digest: str = ""


def _options(entry: Mapping[str, object], where: str) -> DatasetOptions:
    """Read one dataset's options.

    Parameters
    ----------
    entry : mapping
        The ``[[fitted]]`` or ``[[validators]]`` entry.
    where : str
        The configuration, for messages.

    Returns
    -------
    DatasetOptions
        The options.

    Raises
    ------
    DataError
        If the entry has unknown keys.
    """
    allowed = {
        "key",
        "range_nm",
        "sigma",
        "extra",
        "offset_sigma",
        "offset_group",
        "reference",
        "thermo_optic",
        "temperature_k",
        "note",
    }
    unknown = set(entry) - allowed
    if unknown:
        msg = f"{where}: unknown dataset keys {sorted(unknown)}"
        raise DataError(msg)
    window = entry.get("range_nm", [0.0, math.inf])
    return DatasetOptions(
        key=str(entry["key"]),
        range_nm=(float(window[0]), float(window[1])),  # pyright: ignore[reportIndexIssue, reportArgumentType]
        sigma=None if "sigma" not in entry else float(entry["sigma"]),  # pyright: ignore[reportArgumentType]
        extra=float(entry.get("extra", 0.0)),  # pyright: ignore[reportArgumentType]
        offset_sigma=(
            None if "offset_sigma" not in entry else float(entry["offset_sigma"])  # pyright: ignore[reportArgumentType]
        ),
        offset_group=str(entry.get("offset_group", "")),
        reference=str(entry.get("reference", "")),
        thermo_optic=str(entry.get("thermo_optic", "")),
        temperature_k=(
            None if "temperature_k" not in entry else float(entry["temperature_k"])  # pyright: ignore[reportArgumentType]
        ),
    )


def load_config(name: str | Path) -> FitConfig:
    """Load a solid fit configuration.

    Parameters
    ----------
    name : str or pathlib.Path
        The name (``data/solid/fits/<name>.toml``) or a path.

    Returns
    -------
    FitConfig
        The configuration.

    Raises
    ------
    DataError
        If it is malformed or names an unknown dataset.
    """
    path = Path(name)
    if path.suffix != ".toml":
        path = FIT_DIRECTORY / f"{name}.toml"
    raw = path.read_bytes()
    record = tomllib.loads(raw.decode("utf-8"))
    allowed = {"material", "description", "range_nm", "model", "fitted", "validators"}
    if set(record) - allowed:
        msg = f"{path.name}: unknown top-level keys {sorted(set(record) - allowed)}"
        raise DataError(msg)
    model = record["model"]
    if set(model) - {"start", "free", "priors", "scales", "add"}:
        msg = f"{path.name}: unknown [model] keys"
        raise DataError(msg)
    fitted = tuple(_options(e, path.name) for e in record.get("fitted", []))
    validators = tuple(_options(e, path.name) for e in record.get("validators", []))
    known = _catalogue.by_key()
    for options in (*fitted, *validators):
        if options.key not in known:
            msg = f"{path.name}: unknown dataset {options.key}"
            raise DataError(msg)
    low, high = record["range_nm"]
    priors = {
        str(k): (float(v[0]), float(v[1])) for k, v in model.get("priors", {}).items()
    }
    return FitConfig(
        name=path.stem,
        material=str(record["material"]),
        description=str(record.get("description", "")),
        range_nm=(float(low), float(high)),
        start=str(model["start"]),
        free=tuple(str(p) for p in model["free"]),
        scales={
            str(k): tuple(str(p) for p in v) for k, v in model.get("scales", {}).items()
        },
        added_attributes=tuple(
            str(a) for a in model.get("add", {}).get("attributes", [])
        ),
        added_values={
            str(k): float(v) for k, v in model.get("add", {}).get("values", {}).items()
        },
        priors=priors,
        fitted=fitted,
        validators=validators,
        digest=hashlib.sha256(raw).hexdigest(),
    )


@dataclass(frozen=True)
class _Observations:
    """The fitted points of one dataset."""

    key: str
    wavelength: Float[np.ndarray, " n"]
    value: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    component: str
    offset_sigma: float | None
    offset: str


def _air_index(dataset: Dataset, wavelength: np.ndarray) -> np.ndarray:
    """Return the index of the air a dataset's n is relative to, or 1 if absolute.

    Parameters
    ----------
    dataset : Dataset
        A refractive-index dataset: ``parameters.index_medium`` ``air`` (with
        ``air_temperature_k``, else the sample's ``temperature_k``, and
        ``air_pressure_pa``, else 101 325 Pa) or anything else for an absolute index.
    wavelength : numpy.ndarray
        Vacuum wavelengths, m.

    Returns
    -------
    numpy.ndarray
        The index of that air at each wavelength.

    Raises
    ------
    DataError
        If the air's temperature is not recorded.
    """
    parameters = dataset.parameters
    if parameters.get("index_medium") != "air":
        return np.ones_like(wavelength)
    temperature = parameters.get("air_temperature_k", parameters.get("temperature_k"))
    if temperature is None:
        msg = f"{dataset.key}: an index relative to air needs air_temperature_k"
        raise DataError(msg)
    pressure = float(parameters.get("air_pressure_pa", STD_ATMOSPHERE))  # pyright: ignore[reportArgumentType]
    state = State(pressure, float(temperature))  # pyright: ignore[reportArgumentType]
    return np.asarray(1 + air_refractivity(1 / wavelength, state))


def _temperature_shift(
    dataset: Dataset, options: DatasetOptions, wavelength: np.ndarray
) -> np.ndarray:
    """Return n(T_sample) - n(T_fit) from a measured thermo-optic coefficient.

    Parameters
    ----------
    dataset : Dataset
        The index dataset, with ``parameters.temperature_k``.
    options : DatasetOptions
        Its options: ``thermo_optic`` and ``temperature_k``.
    wavelength : numpy.ndarray
        The dataset's vacuum wavelengths, m.

    Returns
    -------
    numpy.ndarray
        The shift dn/dT (T_sample - T_fit) at each wavelength, dn/dT interpolated in
        wavelength (and held constant beyond its ends).

    Raises
    ------
    DataError
        If a temperature is missing or the named dataset is not a thermo-optic one.
    """
    coefficient = _catalogue.by_key().get(options.thermo_optic)
    sample = dataset.parameters.get("temperature_k")
    if (
        coefficient is None
        or coefficient.quantity != "thermo_optic_coefficient"
        or sample is None
        or options.temperature_k is None
    ):
        msg = (
            f"{dataset.key}: thermo_optic needs a thermo_optic_coefficient dataset, "
            "the sample's temperature_k and the fit's temperature_k"
        )
        raise DataError(msg)
    grid = coefficient.vacuum_wavelength()
    order = np.argsort(grid)
    slope = np.interp(wavelength, grid[order], coefficient.scaled_values()[order])
    return slope * (float(sample) - options.temperature_k)  # pyright: ignore[reportArgumentType]


def _observations(options: DatasetOptions) -> _Observations:
    """Build a dataset's observations: n, or log k, with their uncertainties.

    Parameters
    ----------
    options : DatasetOptions
        How the dataset enters.

    Returns
    -------
    _Observations
        The points inside the window.

    Raises
    ------
    DataError
        If the dataset's quantity cannot be fitted or a point has no uncertainty.
    """
    dataset: Dataset = _catalogue.dataset(options.key)
    wavelength = dataset.vacuum_wavelength()
    values = dataset.scaled_values()
    stated = dataset.standard_uncertainty()
    low, high = options.range_nm
    keep = (
        (wavelength >= low * 1e-9) & (wavelength <= high * 1e-9) & np.isfinite(values)
    )
    if dataset.quantity == "absorption_coefficient":
        # k = alpha lambda / (4 pi) for the intensity's natural-log coefficient; a
        # relative uncertainty carries over to k unchanged.
        per_metre = 100.0 if dataset.unit == "cm-1" else 1.0
        natural = math.log(10) if str(dataset.parameters.get("base")) == "10" else 1.0
        factor = per_metre * natural * wavelength / (4 * math.pi)
        values = values * factor
        stated = stated * factor
    if dataset.quantity in {"extinction_coefficient", "absorption_coefficient"}:
        keep &= values > 0
        positive = np.where(values > 0, values, 1.0)
        component = "log_k"
        observed = np.log(positive)
        sigma = np.where(np.isfinite(stated), stated / positive, np.nan)
    elif dataset.quantity in {"refractive_index", "refractive_index_deviation"}:
        component = "n"
        air = _air_index(dataset, wavelength)
        observed = values * air
        sigma = stated
        if dataset.quantity == "refractive_index_deviation":
            # The reference formula is itself absolute (a formula relative to air
            # carries its air, TC-398), so only the deviation is scaled here.
            reference = published_models()[options.reference]
            observed = observed + np.real(
                reference.refractive_index(wavelength, warn=False)
            )
        if options.thermo_optic:
            observed = observed - _temperature_shift(dataset, options, wavelength)
    else:
        msg = f"{dataset.key}: a {dataset.quantity} cannot be fitted"
        raise DataError(msg)
    if options.sigma is not None:
        sigma = np.where(np.isfinite(sigma), sigma, options.sigma)
    sigma = np.sqrt(sigma**2 + options.extra**2)
    if not np.all(np.isfinite(sigma[keep]) & (sigma[keep] > 0)):
        msg = f"{dataset.key}: points without an uncertainty; set sigma in the config"
        raise DataError(msg)
    return _Observations(
        key=dataset.key,
        wavelength=wavelength[keep],
        value=observed[keep],
        sigma=sigma[keep],
        component=component,
        offset_sigma=options.offset_sigma if component == "n" else None,
        offset=(options.offset_group or dataset.key) if component == "n" else "",
    )


@dataclass(frozen=True)
class FitResult:
    """A fitted model with its statistics.

    Attributes
    ----------
    config : FitConfig
        The configuration.
    attributes : tuple of str
        The UDM's attributes.
    values : dict of str to float
        Every UDM parameter, fitted or fixed.
    covariance : numpy.ndarray
        The covariance of the free parameters, then of the group scales, then of the
        offsets (in their own units), Birge-scaled when the reduced chi-squared
        exceeds 1.
    offsets : dict of str to float
        The fitted n offset of each dataset that has one.
    chi2 : float
        The total chi-squared, priors included.
    dof : int
        Degrees of freedom.
    dataset_chi2 : dict of str to tuple
        Points and chi-squared of each fitted dataset.
    """

    config: FitConfig
    attributes: tuple[str, ...]
    values: Mapping[str, float]
    covariance: Float[np.ndarray, "p p"]
    offsets: Mapping[str, float]
    chi2: float
    dof: int
    dataset_chi2: Mapping[str, tuple[int, float]]

    @property
    def model(self) -> UniversalDispersionModel:
        """The fitted UDM.

        Returns
        -------
        UniversalDispersionModel
            The model at the fitted parameters.
        """
        return from_parameters(self.attributes, self.values)


def _index(model: UniversalDispersionModel, wavelength: np.ndarray) -> np.ndarray:
    """Return n + ik of a UDM at vacuum wavelengths.

    Parameters
    ----------
    model : UniversalDispersionModel
        The model.
    wavelength : numpy.ndarray
        Vacuum wavelength, m.

    Returns
    -------
    numpy.ndarray
        The complex index.
    """
    energy = 2 * math.pi * C * HBAR / (E_CHARGE * wavelength)
    return complex_index(model.permittivity(energy))


def _start(config: FitConfig) -> tuple[tuple[str, ...], dict[str, float]]:
    """Return the starting model's attributes and parameters, with any added terms.

    Parameters
    ----------
    config : FitConfig
        The configuration.

    Returns
    -------
    tuple
        The attributes and the parameters.

    Raises
    ------
    DataError
        If the start is not a curated UDM formula, or an added attribute repeats one
        of the start model's.
    """
    source, _, formula = config.start.partition("/")
    record = _catalogue.records().get(source)
    entries = [] if record is None else record.formulas
    terms = next((f.get("terms") for f in entries if f.get("id") == formula), None)
    if not isinstance(terms, Mapping) or "attributes" not in terms:
        msg = f"{config.name}: the start {config.start!r} is not a curated UDM formula"
        raise DataError(msg)
    attributes = parse_attributes(str(terms["attributes"]), config.start)
    values = {str(k): float(v) for k, v in terms.items() if k != "attributes"}
    added = {a.split("=")[0] for a in config.added_attributes}
    if added & {a.split("=")[0] for a in attributes}:
        msg = f"{config.name}: added terms {sorted(added)} repeat the start model's"
        raise DataError(msg)
    return (*attributes, *config.added_attributes), {**values, **config.added_values}


def _unpack(
    config: FitConfig, start: Mapping[str, float], x: np.ndarray
) -> dict[str, float]:
    """Return the full parameter set for the fitted vector.

    Parameters
    ----------
    config : FitConfig
        The configuration.
    start : mapping
        The starting parameters.
    x : numpy.ndarray
        The logarithms of the free parameters (every UDM parameter is positive), then
        the logarithms of the group scales.

    Returns
    -------
    dict of str to float
        Every parameter.
    """
    values = dict(start)
    for name, value in zip(config.free, x[: len(config.free)], strict=True):
        values[name] = math.exp(value)
    logs = x[len(config.free) : len(config.free) + len(config.scales)]
    for members, log_scale in zip(config.scales.values(), logs, strict=True):
        for name in members:
            values[name] = start[name] * math.exp(log_scale)
    return values


def fit(config: FitConfig) -> FitResult:
    """Fit the configuration's free UDM parameters to its datasets.

    Parameters
    ----------
    config : FitConfig
        The configuration.

    Returns
    -------
    FitResult
        The fitted model, covariance and statistics.
    """
    attributes, start = _start(config)
    observations = [_observations(o) for o in config.fitted]
    # One offset per group (a dataset alone is its own group), with the group's prior.
    offset_priors: dict[str, float] = {}
    for obs in observations:
        if obs.offset_sigma is not None:
            offset_priors.setdefault(obs.offset, obs.offset_sigma)
    offset_keys = list(offset_priors)
    x0 = np.array(
        [math.log(start[p]) for p in config.free]
        + [0.0] * len(config.scales)
        + [0.0] * len(offset_keys)
    )
    n_model = len(config.free) + len(config.scales)

    def residuals(x: np.ndarray) -> np.ndarray:
        """Return the weighted residuals of the data and the priors.

        Parameters
        ----------
        x : numpy.ndarray
            The free parameters, then the datasets' offsets.

        Returns
        -------
        numpy.ndarray
            (model - data) / sigma for every point, then the priors' terms.
        """
        values = _unpack(config, start, x)
        model = from_parameters(attributes, values)
        offsets = dict(zip(offset_keys, x[n_model:], strict=True))
        parts = []
        for obs in observations:
            index = _index(model, obs.wavelength)
            if obs.component == "n":
                predicted = np.real(index) + offsets.get(obs.offset, 0.0)
            else:
                predicted = np.log(np.maximum(np.imag(index), 1e-300))
            parts.append((predicted - obs.value) / obs.sigma)
        for key, sigma in offset_priors.items():
            parts.append(np.array([offsets[key] / sigma]))
        for name, (value, sigma) in config.priors.items():
            parts.append(np.array([(values[name] - value) / sigma]))
        return np.concatenate(parts)

    # The natural scale of each unknown: 1 for the logarithms of the model's
    # parameters, the prior width for an offset. (Scaling by the Jacobian instead takes
    # an enormous step in a parameter the data barely touch, such as an absorption
    # band outside the fitted range.)
    x_scale = np.array([1.0] * n_model + [offset_priors[k] for k in offset_keys])
    solution = optimize.least_squares(residuals, x0, method="trf", x_scale=x_scale)
    r = solution.fun
    jac = np.asarray(solution.jac, dtype=float)
    chi2 = float(r @ r)
    dof = max(len(r) - len(x0), 1)
    covariance = np.linalg.pinv(jac.T @ jac)
    if chi2 / dof > 1:
        covariance *= chi2 / dof
    # Back to the parameters' own units (the chain rule of the logarithms).
    values = _unpack(config, start, solution.x)
    scale_values = np.exp(solution.x[len(config.free) : n_model])
    factors = np.array(
        [values[p] for p in config.free] + list(scale_values) + [1.0] * len(offset_keys)
    )
    covariance = covariance * np.outer(factors, factors)
    offsets = dict(zip(offset_keys, solution.x[n_model:], strict=True))
    dataset_chi2 = {}
    position = 0
    for obs in observations:
        count = len(obs.value)
        part = r[position : position + count]
        dataset_chi2[obs.key] = (count, float(part @ part))
        position += count
    return FitResult(
        config=config,
        attributes=attributes,
        values=values,
        covariance=covariance,
        offsets={k: float(v) for k, v in offsets.items()},
        chi2=chi2,
        dof=dof,
        dataset_chi2=dataset_chi2,
    )


def validation(result: FitResult) -> dict[str, dict[str, float]]:
    """Compare the fitted model with each validator dataset.

    Parameters
    ----------
    result : FitResult
        The fit.

    Returns
    -------
    dict of str to dict
        Per validator: the points, chi-squared and mean deviation (n: absolute;
        log k: relative).
    """
    model = result.model
    report: dict[str, dict[str, float]] = {}
    for options in result.config.validators:
        obs = _observations(options)
        index = _index(model, obs.wavelength)
        if obs.component == "n":
            predicted = np.real(index)
        else:
            predicted = np.log(np.maximum(np.imag(index), 1e-300))
        deviation = predicted - obs.value
        report[obs.key] = {
            "n": float(len(deviation)),
            "chi2": float(np.sum((deviation / obs.sigma) ** 2)),
            "mean_deviation": float(np.mean(deviation)),
        }
    return report


def dumps(result: FitResult) -> str:
    """Serialise a fit result as JSON text.

    Parameters
    ----------
    result : FitResult
        The result.

    Returns
    -------
    str
        Indented JSON: the full parameter set, the covariance of the free parameters,
        the statistics and the configuration's digest.
    """
    record = {
        "name": result.config.name,
        "material": result.config.material,
        "config_sha256": result.config.digest,
        "range_nm": list(result.config.range_nm),
        "start": result.config.start,
        "attributes": list(result.attributes),
        "values": dict(result.values),
        "free": list(result.config.free),
        "scales": {k: list(v) for k, v in result.config.scales.items()},
        "covariance": result.covariance.tolist(),
        "offsets": dict(result.offsets),
        "chi2": result.chi2,
        "dof": result.dof,
        "dataset_chi2": {k: list(v) for k, v in result.dataset_chi2.items()},
        "validation": validation(result),
    }
    return json.dumps(record, indent=1) + "\n"

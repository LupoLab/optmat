# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Fit oscillator-strength models to the curated data.

A fit is specified by a configuration, ``data/gas/fits/<name>.toml``, which names

- the model: the resonance lines kept at their spectroscopic energies (NIST), each
  with a prior on its oscillator strength, and the number of *effective* oscillators
  that stand for the rest of the distribution (the higher Rydberg lines and the
  continuum), above a stated lowest energy;
- optional priors on dipole sums (the static polarizability ``S(-2)``, the
  Thomas-Reiche-Kuhn sum ``S(0)``);
- the datasets fitted, each with its options: a free scale (a dataset whose absolute
  level is less certain than its shape, like Mansfield & Peck's helium), an extra
  relative uncertainty added in quadrature, and a wavelength window;
- the datasets used only as validators.

The fit is a generalised least-squares problem in the polarizability (atomic units).
Each fitted value contributes ``(s_d alpha_model(omega) - alpha)/sigma``, where ``s_d``
is its dataset's scale: 1, free, or constrained by the dataset's stated systematic
uncertainty (a prior residual ``(s_d - 1)/sigma_sys``). Line strengths enter as
``(f - f_prior)/sigma_f``. Effective strengths and energies are fitted as logarithms so
that they stay positive and above the lowest allowed energy.

The parameter covariance is ``(J^T J)^{-1}`` at the solution; when the reduced
chi-squared exceeds one it is also given scaled by it (the Birge ratio), and the
uncertainty band of the polarizability follows by linear propagation.
"""

import dataclasses
import hashlib
import json
import math
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import numpy as np
from jaxtyping import Float
from scipy.optimize import least_squares

from optmat.constants import AU_POLARIZABILITY, BOHR_RADIUS
from optmat.datasets import Dataset, datasets
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY, reduction
from optmat.gas.absorption import MeasuredAbsorption, measured_absorption
from optmat.gas.gases import Gas, as_gas
from optmat.gas.infrared import (
    band_systems,
    derived_digest,
    infrared_bands,
    qualified_windows,
)
from optmat.gas.models.oscillator import OscillatorModel, angular_frequency_au

__all__ = [
    "FIT_DIRECTORY",
    "HARTREE_NM",
    "DatasetOptions",
    "FitConfig",
    "FitResult",
    "Line",
    "Observation",
    "dumps",
    "fit",
    "load_config",
    "observations",
    "result_record",
    "validation",
]

FIT_DIRECTORY: Final = GAS_DATA_DIRECTORY / "fits"
"""Fit configurations, ``<name>.toml``."""

HARTREE_NM: Final = 45.56335252767
"""Vacuum wavelength of a one-hartree photon, nm (2 pi c / (E_h / hbar))."""


@dataclass(frozen=True)
class Line:
    """A resonance line kept at its spectroscopic energy.

    Attributes
    ----------
    label : str
        The transition.
    wavelength_nm : float
        Vacuum wavelength, nm.
    f_prior : float
        Prior oscillator strength.
    f_sigma : float
        Its standard uncertainty.
    """

    label: str
    wavelength_nm: float
    f_prior: float
    f_sigma: float

    @property
    def energy(self) -> float:
        """The line's energy.

        Returns
        -------
        float
            The photon energy, hartree.
        """
        return HARTREE_NM / self.wavelength_nm


@dataclass(frozen=True)
class DatasetOptions:
    """How one dataset enters a fit.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    free_scale : bool
        Fit the dataset's scale without a prior.
    extra_relative : float
        Relative uncertainty added in quadrature to every value.
    range_nm : tuple of float
        Vacuum-wavelength window of the values used, nm.
    electric : bool
        The values are electric polarizabilities (ab initio theory, dielectric
        measurements), compared with the model without its ``chi_0``. Cauchy moments are
        always electric.
    infrared : bool
        The values include the vibrational (infrared) contribution, as measured
        refractivities and dielectric constants do; False for electronic theory. Only
        used when the model has infrared terms; Cauchy moments are always electronic.
    correction_au : tuple of float or None
        ``(value, sigma)`` added to every value, atomic units, with ``sigma`` added in
        quadrature to its uncertainty: a term the values lack and the model includes,
        such as the zero-point vibrational average of a polarizability computed at the
        equilibrium geometry.
    shape_relative : float or None
        With ``free_scale``: the relative uncertainty of the values' shape, used in
        place of the stated one when that states the absolute level (which the free
        scale absorbs) and the register gives the shape's own precision.
    temperature : tuple of float or None
        ``(measured_k, coefficient, sigma)``: the values were measured at
        ``measured_k`` (K), and the molecule's zero-density polarizability changes by
        ``coefficient`` (per K, relative; its uncertainty ``sigma``) with temperature.
        The values are brought to the fit's ``reference_temperature`` by the factor
        ``1 + coefficient (T_ref - measured_k)``, and ``sigma |T_ref - measured_k|``
        is added in quadrature to their relative uncertainty (TC-730).
    """

    key: str
    free_scale: bool = False
    extra_relative: float = 0.0
    range_nm: tuple[float, float] = (0.0, math.inf)
    electric: bool = False
    infrared: bool = True
    correction_au: tuple[float, float] | None = None
    shape_relative: float | None = None
    temperature: tuple[float, float, float] | None = None


@dataclass(frozen=True)
class FitConfig:
    """A fit's specification.

    Attributes
    ----------
    name : str
        The configuration's file stem, and the fitted model's name.
    gas : Gas
        The species.
    description : str
        What the fit is for.
    lines : tuple of Line
        The fixed-energy resonance lines.
    effective : int
        The number of effective oscillators.
    effective_min_energy : float
        The lowest energy an effective oscillator may take, hartree.
    sum_priors : dict of int to tuple of float
        Priors ``{k: (value, sigma)}`` on dipole sums ``S(k)``, atomic units.
    chi0 : float
        The atom's magnetic susceptibility, atomic units: the model's refractive
        polarizability is its electric sum plus ``chi0``.
    fitted : tuple of DatasetOptions
        The datasets fitted.
    validators : tuple of DatasetOptions
        Datasets compared with the fit but not fitted (with their ``electric`` flag).
    range_nm : tuple of float
        The wavelength range over which the fitted model is declared valid, nm.
    digest : str
        SHA-256 of the configuration file, recorded with the result.
    starts : int
        The number of starting points :func:`fit` tries (``[model] starts``, default 6).
    infrared : tuple of tuple
        Fixed vibrational terms ``(label, strength, energy in hartree)``: the gas's
        HITRAN bands (:mod:`optmat.gas.infrared`) up to
        ``[model] infrared_max_wavenumber``; empty if that key is absent.
    infrared_max_wavenumber : float or None
        The cutoff, cm^-1.
    fixed : tuple of tuple
        Fixed electronic terms ``(label, strength, energy in hartree)`` from
        ``[[model.fixed]]`` (a published distribution's interval strengths above the
        measured range, say): always included, outside the fitted parameters.
    reference_temperature : float or None
        The temperature (K) the model's polarizability refers to (``[model]
        reference_temperature_k``): datasets with a ``temperature`` option are brought
        to it. None where no reference is declared (the atoms: their polarizability does
        not depend on temperature at the package's accuracy).
    retardation : float
        The coefficient ``c`` of the retardation term ``c omega^2`` of the refractive
        polarizability (``[model] retardation_au``, atomic units, omega in hartree):
        added, like ``chi0``, to refractive observations and not to electric ones or
        to moments. Known for helium only (TC-705); zero otherwise.
    absorption : MeasuredAbsorption or None
        A fixed measured-absorption term (``[model] absorption = true``, from
        ``data/gas/absorption/<gas>.toml``, smoothed to ``[model] absorption_fwhm_cm``
        if given): its Kramers-Kronig polarizability is added to every prediction and
        its dipole sums to the sum priors.
    band_nodes : int or None
        ``[model] band_nodes`` (OEP 3): the infrared terms are the derived rules of
        every band system (``<gas>_band_rules.csv``), at most this many nodes each,
        in place of ``infrared_max_wavenumber``'s band poles.
    rotational : tuple of tuple
        The rotational system's terms, with ``band_nodes``: a term class of its own,
        added to refractive observations and left out of electric ones and moments.
    windows : bool
        ``[model] windows = true``: the model's qualified intervals are the gas's
        derived windows (``<gas>_windows.csv``), with each gap's reason.
    intervals_nm : tuple of tuple of float
        The qualified intervals, nm: the derived windows, or ``range_nm`` alone.
    gap_reasons : tuple of str
        The reason of each gap between two intervals.
    derived_sha256 : str
        The SHA-256 of the derived tables the fit reads (``band_nodes`` or
        ``windows``); empty otherwise.
    """

    name: str
    gas: Gas
    description: str
    lines: tuple[Line, ...]
    effective: int
    effective_min_energy: float
    sum_priors: Mapping[int, tuple[float, float]]
    chi0: float
    fitted: tuple[DatasetOptions, ...]
    validators: tuple[DatasetOptions, ...]
    range_nm: tuple[float, float]
    digest: str
    starts: int = 6
    infrared: tuple[tuple[str, float, float], ...] = ()
    infrared_max_wavenumber: float | None = None
    fixed: tuple[tuple[str, float, float], ...] = ()
    retardation: float = 0.0
    reference_temperature: float | None = None
    absorption: MeasuredAbsorption | None = None
    band_nodes: int | None = None
    rotational: tuple[tuple[str, float, float], ...] = ()
    windows: bool = False
    intervals_nm: tuple[tuple[float, float], ...] = ()
    gap_reasons: tuple[str, ...] = ()
    derived_sha256: str = ""

    @property
    def fixed_terms(self) -> tuple[tuple[str, float, float], ...]:
        """Every term outside the fitted parameters: rotational, infrared, the rest.

        Returns
        -------
        tuple of tuple
            ``(label, strength, energy in hartree)``.
        """
        return self.rotational + self.infrared + self.fixed


def _temperature(
    entry: Mapping[str, object], where: str
) -> tuple[float, float, float] | None:
    """Read a dataset's ``temperature = {measured_k, coefficient, sigma}`` option.

    Parameters
    ----------
    entry : mapping
        The ``[[fitted]]`` or validator entry.
    where : str
        For error messages.

    Returns
    -------
    tuple of float or None
        ``(measured_k, coefficient, sigma)``, or None without the option.

    Raises
    ------
    DataError
        If the option lacks one of its three keys or has others.
    """
    spec = entry.get("temperature")
    if spec is None:
        return None
    keys = {"measured_k", "coefficient", "sigma"}
    if not isinstance(spec, Mapping) or set(spec) != keys:
        msg = f"{where}: temperature needs exactly {sorted(keys)}"
        raise DataError(msg)
    return (
        float(spec["measured_k"]),  # pyright: ignore[reportArgumentType]
        float(spec["coefficient"]),  # pyright: ignore[reportArgumentType]
        float(spec["sigma"]),  # pyright: ignore[reportArgumentType]
    )


def _pair(value: object, where: str) -> tuple[float, float]:
    """Read a ``[low, high]`` pair ("inf" allowed).

    Parameters
    ----------
    value : object
        The TOML value.
    where : str
        For error messages.

    Returns
    -------
    tuple of float
        The pair.

    Raises
    ------
    DataError
        If the value is not a pair.
    """
    if not isinstance(value, Sequence) or len(value) != 2:
        msg = f"{where}: expected [low, high]"
        raise DataError(msg)
    low, high = (float(v) for v in value)  # pyright: ignore[reportArgumentType]
    return low, high


def _infrared_terms(
    gas: Gas, cutoff: float | None
) -> tuple[tuple[str, float, float], ...]:
    """Return a gas's HITRAN bands up to a wavenumber as fixed oscillator terms.

    Parameters
    ----------
    gas : Gas
        The species.
    cutoff : float or None
        The largest pole wavenumber kept, cm^-1; None for no infrared terms.

    Returns
    -------
    tuple of tuple
        ``(label, strength, energy in hartree)`` per band kept.
    """
    if cutoff is None:
        return ()
    return tuple(
        (f"infrared {band.label}", band.oscillator_strength, band.energy_au)
        for band in infrared_bands(gas)
        if band.wavenumber <= cutoff
    )


def _derived_terms(gas: Gas, model: Mapping[str, object], name: str) -> dict[str, Any]:
    """Return what a configuration takes from the gas's derived tables (OEP 3).

    Parameters
    ----------
    gas : Gas
        The species.
    model : Mapping
        The ``[model]`` table.
    name : str
        The configuration file, for messages.

    Returns
    -------
    dict
        ``band_nodes``, ``infrared`` and ``rotational`` terms, ``windows``,
        ``intervals_nm``, ``gap_reasons`` and ``derived_sha256``; empty values where
        the configuration uses neither key.

    Raises
    ------
    DataError
        If ``band_nodes`` is given with ``infrared_max_wavenumber``, or is below the
        tables' node count.
    """
    nodes = model.get("band_nodes")
    windows = bool(model.get("windows", False))
    out: dict[str, Any] = {
        "band_nodes": None,
        "infrared": (),
        "rotational": (),
        "windows": windows,
        "intervals_nm": (),
        "gap_reasons": (),
        "derived_sha256": "",
    }
    if nodes is not None:
        if "infrared_max_wavenumber" in model:
            msg = f"{name}: band_nodes replaces infrared_max_wavenumber; give one"
            raise DataError(msg)
        systems = band_systems(gas)
        needed = max(len(s.strengths) for s in systems)
        if int(str(nodes)) < needed:
            msg = (
                f"{name}: band_nodes = {nodes}, but {gas.value}'s tables have rules of "
                f"{needed} nodes; rederive them with that cap (tools/hitran/windows.py)"
            )
            raise DataError(msg)
        terms = {
            kind: tuple(
                (f"{kind} {s.system}.{k}", f, e)
                for s in systems
                if (s.treatment == "rotational") == (kind == "rotational")
                for k, (f, e) in enumerate(zip(s.strengths, s.energies_au, strict=True))
            )
            for kind in ("infrared", "rotational")
        }
        out |= {"band_nodes": int(str(nodes)), **terms}
    if windows:
        intervals, reasons = qualified_windows(gas)
        out |= {"intervals_nm": intervals, "gap_reasons": reasons}
    if nodes is not None or windows:
        out["derived_sha256"] = derived_digest(gas)
    return out


def _interval_terms(spec: object, path: Path) -> tuple[tuple[str, float, float], ...]:
    """Return a distribution's interval strengths above an energy as fixed terms.

    Each interval of an ``interval_oscillator_strength`` table in eV becomes one pole
    at the interval's geometric-mean energy; an interval open to infinity is taken to
    end at twice its lower edge.

    Parameters
    ----------
    spec : object
        ``{dataset = "source/dataset", from_ev = <lower edge>}``, or None.
    path : pathlib.Path
        The configuration, for error messages.

    Returns
    -------
    tuple of tuple
        ``(label, strength, energy in hartree)`` per interval at or above ``from_ev``.

    Raises
    ------
    DataError
        If the specification or dataset is malformed.
    """
    if spec is None:
        return ()
    if not isinstance(spec, Mapping) or "dataset" not in spec:
        msg = f"{path.name}: fixed_intervals needs a dataset"
        raise DataError(msg)
    by_key = {d.key: d for d in datasets()}
    key = str(spec["dataset"])
    if (
        key not in by_key
        or by_key[key].quantity != "interval_oscillator_strength"
        or by_key[key].abscissa_unit != "eV"
    ):
        msg = f"{path.name}: {key!r} is not an interval table against eV"
        raise DataError(msg)
    table = by_key[key]
    lower = table.abscissa_values()
    strengths = table.scaled_values()
    upper = table.column_values("x_upper")
    upper = np.where(np.isfinite(upper), upper, 2 * lower)
    start = float(spec.get("from_ev", 0.0))  # pyright: ignore[reportArgumentType]
    hartree_per_ev = HARTREE_NM / 1239.8419843320026  # (hc in eV nm)
    return tuple(
        (
            f"{key} {low:g}-{high:g} eV",
            float(f),
            float(np.sqrt(low * high)) * hartree_per_ev,
        )
        for low, high, f in zip(lower, upper, strengths, strict=True)
        if low >= start
    )


def load_config(name: str | Path) -> FitConfig:
    """Load a fit configuration.

    Parameters
    ----------
    name : str or pathlib.Path
        The configuration's name (``data/gas/fits/<name>.toml``) or path.

    Returns
    -------
    FitConfig
        The configuration.

    Raises
    ------
    DataError
        If the configuration is malformed or names an unknown dataset.
    """
    path = Path(name)
    if path.suffix != ".toml":
        path = FIT_DIRECTORY / f"{name}.toml"
    raw = path.read_bytes()
    record = tomllib.loads(raw.decode("utf-8"))
    model = record["model"]
    known = {d.key for d in datasets()}
    allowed = {"gas", "description", "range_nm", "model", "fitted", "validators"}
    unknown_keys = set(record) - allowed
    if unknown_keys:
        msg = f"{path.name}: unknown top-level keys {sorted(unknown_keys)}"
        raise DataError(msg)
    model_keys = {
        "lines",
        "effective",
        "effective_min_energy_hartree",
        "sum_priors",
        "chi0_au",
        "retardation_au",
        "reference_temperature_k",
        "starts",
        "infrared_max_wavenumber",
        "band_nodes",
        "windows",
        "absorption",
        "absorption_fwhm_cm",
        "fixed",
        "fixed_intervals",
    }
    if set(model) - model_keys:
        msg = f"{path.name}: unknown [model] keys {sorted(set(model) - model_keys)}"
        raise DataError(msg)
    options_keys = {
        "key",
        "free_scale",
        "extra_relative",
        "range_nm",
        "electric",
        "infrared",
        "correction_au",
        "shape_relative",
        "temperature",
    }
    fitted = []
    for entry in record.get("fitted", []):
        if set(entry) - options_keys:
            # A top-level key written after the [[fitted]] tables lands in the last
            # one (TOML): refuse it rather than lose it.
            msg = (
                f"{path.name}: unknown keys {sorted(set(entry) - options_keys)} in "
                f"the [[fitted]] entry for {entry.get('key')!r}"
            )
            raise DataError(msg)
        if entry["key"] not in known:
            msg = f"{path.name}: unknown dataset {entry['key']!r}"
            raise DataError(msg)
        fitted.append(
            DatasetOptions(
                key=entry["key"],
                free_scale=bool(entry.get("free_scale", False)),
                extra_relative=float(entry.get("extra_relative", 0.0)),
                range_nm=_pair(entry.get("range_nm", [0.0, "inf"]), path.name),
                electric=bool(entry.get("electric", False)),
                infrared=bool(entry.get("infrared", True)),
                correction_au=(
                    _pair(entry["correction_au"], path.name)
                    if "correction_au" in entry
                    else None
                ),
                shape_relative=(
                    float(entry["shape_relative"])
                    if "shape_relative" in entry
                    else None
                ),
                temperature=_temperature(entry, path.name),
            )
        )
    validators = tuple(
        DatasetOptions(key=entry)
        if isinstance(entry, str)
        else DatasetOptions(
            key=entry["key"],
            electric=bool(entry.get("electric")),
            infrared=bool(entry.get("infrared", True)),
            range_nm=_pair(entry["range_nm"], path.name)
            if "range_nm" in entry
            else (0.0, math.inf),
            temperature=_temperature(entry, path.name),
        )
        for entry in record.get("validators", [])
    )
    unknown = {v.key for v in validators} - known
    if unknown:
        msg = f"{path.name}: unknown validator datasets {sorted(unknown)}"
        raise DataError(msg)
    gas = as_gas(record["gas"])
    cutoff = model.get("infrared_max_wavenumber")
    cutoff = None if cutoff is None else float(cutoff)
    derived = _derived_terms(gas, model, path.name)
    absorption = None
    if model.get("absorption", False):
        fwhm = model.get("absorption_fwhm_cm")
        absorption = measured_absorption(
            gas, fwhm_cm=None if fwhm is None else float(fwhm)
        )
    return FitConfig(
        name=path.stem,
        gas=gas,
        description=record.get("description", ""),
        lines=tuple(
            Line(
                label=line["label"],
                wavelength_nm=float(line["wavelength_nm"]),
                f_prior=float(line["f_prior"]),
                f_sigma=float(line["f_sigma"]),
            )
            for line in model.get("lines", [])
        ),
        effective=int(model["effective"]),
        effective_min_energy=float(model["effective_min_energy_hartree"]),
        sum_priors={
            int(k): _pair(v, path.name) for k, v in model.get("sum_priors", {}).items()
        },
        chi0=float(model.get("chi0_au", 0.0)),
        retardation=float(model.get("retardation_au", 0.0)),
        reference_temperature=(
            float(model["reference_temperature_k"])
            if "reference_temperature_k" in model
            else None
        ),
        fitted=tuple(fitted),
        validators=validators,
        range_nm=_pair(record["range_nm"], path.name),
        digest=hashlib.sha256(raw).hexdigest(),
        starts=int(model.get("starts", 6)),
        infrared=derived["infrared"] or _infrared_terms(gas, cutoff),
        infrared_max_wavenumber=cutoff,
        band_nodes=derived["band_nodes"],
        rotational=derived["rotational"],
        windows=derived["windows"],
        intervals_nm=derived["intervals_nm"] or (_pair(record["range_nm"], path.name),),
        gap_reasons=derived["gap_reasons"],
        derived_sha256=derived["derived_sha256"],
        fixed=tuple(
            (
                str(term["label"]),
                float(term["f"]),
                HARTREE_NM / float(term["wavelength_nm"]),
            )
            for term in model.get("fixed", [])
        )
        + _interval_terms(model.get("fixed_intervals"), path),
        absorption=absorption,
    )


@dataclass(frozen=True)
class Observation:
    """The fitted values of one dataset, as polarizabilities in atomic units.

    Attributes
    ----------
    key : str
        ``source/dataset``.
    omega : numpy.ndarray
        Photon energies, hartree.
    alpha : numpy.ndarray
        Polarizabilities, atomic units.
    sigma : numpy.ndarray
        One-sigma uncertainties, atomic units.
    systematic : float or None
        The dataset's stated relative systematic uncertainty.
    options : DatasetOptions
        How the dataset enters the fit.
    orders : numpy.ndarray or None
        For Cauchy moments, the order ``k`` of each ``S(k)`` (then ``alpha`` holds the
        sums and ``omega`` is unused); None for polarizabilities.
    harmonic : int
        For a phase-matching density ratio, the harmonic (``omega`` is the fundamental,
        ``alpha`` holds the ratios); 0 otherwise.
    partner_delta : numpy.ndarray or None
        For a ratio, the partner gas's ``alpha(h omega) - alpha(omega)`` from its fixed
        model, atomic units.
    model_in_numerator : bool
        For a ratio, True if the fitted gas's difference is the ratio's numerator
        (``rho_partner / rho_fitted``), False if its denominator.
    """

    key: str
    omega: Float[np.ndarray, " n"]
    alpha: Float[np.ndarray, " n"]
    sigma: Float[np.ndarray, " n"]
    systematic: float | None
    options: DatasetOptions
    orders: Float[np.ndarray, " n"] | None = None
    harmonic: int = 0
    partner_delta: Float[np.ndarray, " n"] | None = None
    model_in_numerator: bool = False


def _ratio(dataset: Dataset, options: DatasetOptions, fitted_gas: Gas) -> Observation:
    """Convert a phase-matching density ratio to an observation of the fitted gas.

    The ratio ``rho_gas / rho_ref = [alpha_ref(h w) - alpha_ref(w)] / [alpha_gas(h w) -
    alpha_gas(w)]`` ties the fitted gas (either ``material`` or ``reference_material``)
    to its partner, whose difference comes from the partner's stored fit (or recommended
    model) and is held fixed: fits are run in an order where the partner comes first.

    Parameters
    ----------
    dataset : Dataset
        A ``phase_match_density_ratio`` dataset.
    options : DatasetOptions
        How it enters the fit.
    fitted_gas : Gas
        The gas being fitted.

    Returns
    -------
    Observation
        The ratios, with the partner's differences.

    Raises
    ------
    DataError
        If neither of the dataset's gases is the fitted one, or a ratio has no stated
        uncertainty.
    """
    # Imported here: the registry imports the fitted models, which import this module's
    # neighbours; the partner is needed only when a ratio is fitted.
    from optmat.gas.models.registry import default_model  # noqa: PLC0415

    reference = as_gas(str(dataset.parameters["reference_material"]))
    if dataset.material == fitted_gas:
        partner, model_in_numerator = reference, False
    elif reference is fitted_gas:
        partner, model_in_numerator = as_gas(dataset.material), True
    else:
        msg = f"{dataset.key}: neither gas is {fitted_gas.value}"
        raise DataError(msg)
    harmonic = int(dataset.parameters.get("harmonic", 2))  # pyright: ignore[reportArgumentType]
    wavelength = dataset.vacuum_wavelength()
    partner_model = default_model(partner)
    delta = (
        partner_model.polarizability(wavelength / harmonic, warn=False)
        - partner_model.polarizability(wavelength, warn=False)
    ) / AU_POLARIZABILITY
    values = dataset.scaled_values()
    stated = dataset.standard_uncertainty()
    extra = options.extra_relative * np.abs(values)
    sigma = np.where(np.isfinite(stated), np.hypot(stated, extra), extra)
    if np.any(~(sigma > 0)):
        msg = f"{dataset.key}: a fitted ratio needs a stated uncertainty"
        raise DataError(msg)
    return Observation(
        key=dataset.key,
        omega=angular_frequency_au(wavelength),
        alpha=values,
        sigma=sigma,
        systematic=None,
        options=options,
        harmonic=harmonic,
        partner_delta=np.asarray(delta),
        model_in_numerator=model_in_numerator,
    )


def _observation(
    dataset: Dataset,
    options: DatasetOptions,
    reference_temperature: float | None = None,
) -> Observation:
    """Convert a dataset to an observation, applying the options' window and extras.

    Parameters
    ----------
    dataset : Dataset
        The dataset (a quantity that determines the polarizability).
    options : DatasetOptions
        How it enters the fit.
    reference_temperature : float or None, optional
        The fit's reference temperature, K; needed when the options carry a
        ``temperature``.

    Returns
    -------
    Observation
        Its values in atomic units.

    Raises
    ------
    DataError
        If a value in the window has no stated uncertainty.
    """
    if dataset.quantity == "cauchy_moment":
        return _moments(dataset, options)
    if dataset.quantity == "rayleigh_cross_section":
        return _rayleigh(dataset, options)
    wavelength = dataset.vacuum_wavelength() * 1e9
    low, high = options.range_nm
    keep = (wavelength >= low) & (wavelength <= high)
    alpha = reduction.polarizability(dataset)[keep] / AU_POLARIZABILITY
    stated = reduction.polarizability_uncertainty(dataset)[keep] / AU_POLARIZABILITY
    if options.shape_relative is not None:
        if not options.free_scale:
            msg = f"{dataset.key}: shape_relative needs free_scale"
            raise DataError(msg)
        stated = options.shape_relative * np.abs(alpha)
    extra = options.extra_relative * np.abs(alpha)
    if options.temperature is not None:
        if reference_temperature is None:
            msg = (
                f"{dataset.key}: a temperature option needs [model] "
                "reference_temperature_k"
            )
            raise DataError(msg)
        measured, coefficient, spread = options.temperature
        # The molecule's own temperature dependence (rotational and vibrational
        # averaging), linear over the few tens of K involved (TC-730).
        gap = reference_temperature - measured
        factor = 1.0 + coefficient * gap
        alpha, stated = alpha * factor, stated * factor
        extra = np.hypot(extra * factor, spread * abs(gap) * np.abs(alpha))
    if options.correction_au is not None:
        value, spread = options.correction_au
        alpha = alpha + value
        extra = np.hypot(extra, spread)
    # A dataset without a stated uncertainty can only be fitted with an assigned one.
    sigma = np.where(np.isfinite(stated), np.hypot(stated, extra), extra)
    if np.any(~(sigma > 0)):
        msg = (
            f"{dataset.key}: a fitted value needs a stated uncertainty or an "
            "extra_relative in the fit configuration"
        )
        raise DataError(msg)
    return Observation(
        key=dataset.key,
        omega=angular_frequency_au(wavelength[keep] * 1e-9),
        alpha=alpha,
        sigma=sigma,
        systematic=dataset.uncertainty.systematic_relative,
        options=options,
    )


def _rayleigh(dataset: Dataset, options: DatasetOptions) -> Observation:
    """Convert Rayleigh cross sections to polarizabilities.

    ``sigma = (8 pi / 3) k^4 (alpha / 4 pi epsilon_0)^2 F_K`` per molecule, with ``F_K``
    the King correction for the anisotropy (``parameters.king_factor``, as the source
    used it; 1 for atoms and spherical tops). The magnitude of the polarizability
    follows, with half the cross section's relative uncertainty. Scattering carries no
    magnetic term, so the observation is electric; it includes the vibrational part.

    Parameters
    ----------
    dataset : Dataset
        A ``rayleigh_cross_section`` dataset (unit cm2).
    options : DatasetOptions
        How it enters the fit.

    Returns
    -------
    Observation
        The polarizabilities, atomic units.

    Raises
    ------
    DataError
        If the dataset has no numeric King factor, lies beyond a resonance (where the
        sign of the polarizability is not determined) or has no stated uncertainty.
    """
    king = dataset.parameters.get("king_factor")
    if isinstance(king, bool) or not isinstance(king, int | float):
        msg = f"{dataset.key}: a Rayleigh cross section needs a numeric king_factor"
        raise DataError(msg)
    if dataset.parameters.get("beyond_resonance"):
        msg = f"{dataset.key}: beyond a resonance the sign of alpha is not determined"
        raise DataError(msg)
    if dataset.unit != "cm2":  # pragma: no cover - the vocabulary fixes the unit
        msg = f"{dataset.key}: unit {dataset.unit!r}, expected cm2"
        raise DataError(msg)
    wavelength = dataset.vacuum_wavelength()
    low, high = options.range_nm
    keep = (wavelength * 1e9 >= low) & (wavelength * 1e9 <= high)
    cross_section = dataset.scaled_values()[keep] * 1e-4  # m^2
    stated = dataset.standard_uncertainty()[keep] * 1e-4
    k = 2 * math.pi / wavelength[keep]
    # alpha / (4 pi epsilon_0) in m^3, then atomic units (a0^3).
    alpha = np.sqrt(cross_section / (8 * math.pi / 3 * k**4 * float(king)))
    alpha = alpha / BOHR_RADIUS**3
    sigma = np.hypot(
        0.5 * alpha * stated / cross_section, options.extra_relative * alpha
    )
    if np.any(~(sigma > 0)):
        msg = f"{dataset.key}: a fitted cross section needs a stated uncertainty"
        raise DataError(msg)
    return Observation(
        key=dataset.key,
        omega=angular_frequency_au(wavelength[keep]),
        alpha=alpha,
        sigma=sigma,
        systematic=dataset.uncertainty.systematic_relative,
        options=dataclasses.replace(options, electric=True),
    )


def _moments(dataset: Dataset, options: DatasetOptions) -> Observation:
    """Convert a dataset of Cauchy moments ``S(k)`` to an observation.

    Parameters
    ----------
    dataset : Dataset
        The dataset (quantity ``cauchy_moment``, orders in the wavelength column).
    options : DatasetOptions
        How it enters the fit.

    Returns
    -------
    Observation
        The sums, with their orders.

    Raises
    ------
    DataError
        If a sum has no stated uncertainty.
    """
    orders = dataset.abscissa_values()
    values = dataset.scaled_values()
    stated = dataset.standard_uncertainty()
    extra = options.extra_relative * np.abs(values)
    sigma = np.where(np.isfinite(stated), np.hypot(stated, extra), extra)
    if np.any(~(sigma > 0)):
        msg = f"{dataset.key}: a fitted sum needs a stated uncertainty"
        raise DataError(msg)
    return Observation(
        key=dataset.key,
        omega=np.zeros_like(values),
        alpha=values,
        sigma=sigma,
        systematic=dataset.uncertainty.systematic_relative,
        options=options,
        orders=orders,
    )


def observations(config: FitConfig) -> tuple[Observation, ...]:
    """Return the fitted observations of a configuration.

    Parameters
    ----------
    config : FitConfig
        The fit.

    Returns
    -------
    tuple of Observation
        One per fitted dataset.
    """
    by_key = {d.key: d for d in datasets()}
    return tuple(
        _ratio(by_key[o.key], o, config.gas)
        if by_key[o.key].quantity == "phase_match_density_ratio"
        else _observation(by_key[o.key], o, config.reference_temperature)
        for o in config.fitted
    )


@dataclass(frozen=True)
class _Layout:
    """Where each parameter sits in the parameter vector."""

    n_lines: int
    n_effective: int
    scaled: tuple[int, ...]  # indices of observations with a fitted scale

    @property
    def size(self) -> int:
        """The number of parameters.

        Returns
        -------
        int
            The length of the parameter vector.
        """
        return self.n_lines + 2 * self.n_effective + len(self.scaled)


def _unpack(
    theta: Float[np.ndarray, " p"], config: FitConfig, layout: _Layout
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return strengths, energies and dataset scales from the parameter vector.

    Parameters
    ----------
    theta : numpy.ndarray
        Line strengths, then log effective strengths, then log excess energies of the
        effective oscillators above the minimum, then the fitted scales minus one.
    config : FitConfig
        The fit.
    layout : _Layout
        The vector's layout.

    Returns
    -------
    tuple of numpy.ndarray
        Strengths, energies (hartree) and scales (one per observation).
    """
    nl, ne = layout.n_lines, layout.n_effective
    f_lines = theta[:nl]
    # Clipped so that a wild trial step of the solver cannot overflow: e^40 hartree
    # and e^40 in strength are far outside any physical solution.
    f_eff = np.exp(np.clip(theta[nl : nl + ne], -40.0, 40.0))
    e_eff = config.effective_min_energy + np.exp(
        np.clip(theta[nl + ne : nl + 2 * ne], -40.0, 40.0)
    )
    energies = np.concatenate([[line.energy for line in config.lines], e_eff])
    return np.concatenate([f_lines, f_eff]), energies, theta[nl + 2 * ne :]


def _alpha(
    strengths: np.ndarray, energies: np.ndarray, omega: np.ndarray
) -> np.ndarray:
    """Evaluate ``sum_i f_i / (E_i^2 - omega^2)``.

    Parameters
    ----------
    strengths : numpy.ndarray
        Oscillator strengths.
    energies : numpy.ndarray
        Energies, hartree.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    numpy.ndarray
        The polarizability at each ``omega``, atomic units.
    """
    return np.sum(strengths[:, None] / (energies[:, None] ** 2 - omega**2), axis=0)


def _infrared(config: FitConfig, omega: np.ndarray) -> np.ndarray:
    """Return the fixed vibrational terms' polarizability.

    Parameters
    ----------
    config : FitConfig
        The fit.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    numpy.ndarray
        ``sum f / (E^2 - omega^2)`` over the infrared terms (zero without them).
    """
    return _terms(config.infrared, omega)


def _terms(
    terms: tuple[tuple[str, float, float], ...], omega: np.ndarray
) -> np.ndarray:
    """Return fixed terms' polarizability.

    Parameters
    ----------
    terms : tuple of tuple
        ``(label, strength, energy in hartree)``.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    numpy.ndarray
        ``sum f / (E^2 - omega^2)`` (zero without terms).
    """
    if not terms:
        return np.zeros_like(omega, dtype=float)
    strengths = np.array([term[1] for term in terms])
    energies = np.array([term[2] for term in terms])
    return _alpha(strengths, energies, omega)


def _measured(config: FitConfig, omega: np.ndarray) -> np.ndarray:
    """Return the fixed electronic terms' and measured absorption's polarizability.

    Parameters
    ----------
    config : FitConfig
        The fit.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    numpy.ndarray
        The fixed electronic terms plus ``Re alpha`` of the measured spectrum (zero
        without either), atomic units.
    """
    fixed = np.zeros_like(omega, dtype=float)
    if config.fixed:
        strengths = np.array([term[1] for term in config.fixed])
        energies = np.array([term[2] for term in config.fixed])
        fixed = _alpha(strengths, energies, omega)
    if config.absorption is None:
        return fixed
    return fixed + config.absorption.complex_polarizability_au(omega)[0]


def _fixed(config: FitConfig, omega: np.ndarray) -> np.ndarray:
    """Return the fixed terms: rotational, infrared and the measured absorption.

    Parameters
    ----------
    config : FitConfig
        The fit.
    omega : numpy.ndarray
        Photon energies, hartree.

    Returns
    -------
    numpy.ndarray
        Their polarizability, atomic units.
    """
    return (
        _terms(config.rotational, omega)
        + _infrared(config, omega)
        + _measured(config, omega)
    )


def _model_values(
    strengths: np.ndarray, energies: np.ndarray, obs: Observation, config: FitConfig
) -> np.ndarray:
    """Return the model's prediction of an observation's values.

    Parameters
    ----------
    strengths : numpy.ndarray
        The fitted (electronic) oscillator strengths.
    energies : numpy.ndarray
        Their energies, hartree.
    obs : Observation
        The observation.
    config : FitConfig
        The fit: its magnetic term and fixed infrared terms.

    Returns
    -------
    numpy.ndarray
        The electronic sums ``S(k)`` for moments; otherwise the polarizability, with
        the magnetic term unless the observation is electric and the infrared terms
        unless it excludes them.
    """
    if obs.orders is not None:
        sums = np.sum(strengths[:, None] * energies[:, None] ** obs.orders, axis=0)
        for _, strength, energy in config.fixed:
            sums = sums + strength * energy**obs.orders
        if config.absorption is not None:
            sums = sums + np.array(
                [config.absorption.sum_rule(float(k)) for k in obs.orders]
            )
        return sums
    if obs.harmonic:
        # chi_0 cancels in the difference; the fixed and retardation terms do not.
        high, low = obs.harmonic * obs.omega, obs.omega
        delta = (
            _alpha(strengths, energies, high)
            + _fixed(config, high)
            - _alpha(strengths, energies, low)
            - _fixed(config, low)
            + config.retardation * (high**2 - low**2)
        )
        assert obs.partner_delta is not None
        if obs.model_in_numerator:
            return delta / obs.partner_delta
        return obs.partner_delta / delta
    # Refraction measures the refractive polarizability: chi_0 and, where known, the
    # retardation term; electric data and moments carry neither (TC-71, TC-705).
    offset = (
        0.0 if obs.options.electric else config.chi0 + config.retardation * obs.omega**2
    )
    if obs.options.infrared:
        offset = offset + _infrared(config, obs.omega)
        # The rotational system is refractive only: electric (static dielectric) data
        # hold the orientation polarizability, or have it removed (OEP 3, TC-751).
        if not obs.options.electric:
            offset = offset + _terms(config.rotational, obs.omega)
    offset = offset + _measured(config, obs.omega)
    return _alpha(strengths, energies, obs.omega) + offset


def _residuals(
    theta: np.ndarray,
    config: FitConfig,
    layout: _Layout,
    data: Sequence[Observation],
) -> np.ndarray:
    """Return the weighted residual vector of the fit.

    Parameters
    ----------
    theta : numpy.ndarray
        Parameters (see :func:`_unpack`).
    config : FitConfig
        The fit.
    layout : _Layout
        The parameter layout.
    data : sequence of Observation
        The fitted observations.

    Returns
    -------
    numpy.ndarray
        Data residuals, then prior residuals.
    """
    strengths, energies, scale_offsets = _unpack(theta, config, layout)
    parts = []
    scales = np.ones(len(data))
    scales[list(layout.scaled)] += scale_offsets
    for scale, obs in zip(scales, data, strict=True):
        model = _model_values(strengths, energies, obs, config)
        parts.append((scale * model - obs.alpha) / obs.sigma)
    for i, line in enumerate(config.lines):
        parts.append(np.array([(strengths[i] - line.f_prior) / line.f_sigma]))
    for index, obs_index in enumerate(layout.scaled):
        obs = data[obs_index]
        if not obs.options.free_scale and obs.systematic:
            parts.append(np.array([scale_offsets[index] / obs.systematic]))
    for k, (value, sigma) in config.sum_priors.items():
        moment = float(np.sum(strengths * energies**k))
        moment += sum(strength * energy**k for _, strength, energy in config.fixed)
        if config.absorption is not None:
            moment += config.absorption.sum_rule(k)
        parts.append(np.array([(moment - value) / sigma]))
    return np.concatenate(parts)


@dataclass(frozen=True)
class FitResult:
    """A fitted model with its uncertainty and statistics.

    Attributes
    ----------
    config : FitConfig
        The fit's configuration.
    model : OscillatorModel
        The fitted model.
    parameters : numpy.ndarray
        The fitted strengths (lines, then effective) and effective energies, hartree.
    covariance : numpy.ndarray
        Their covariance, from ``(J^T J)^{-1}`` (not scaled by the Birge ratio).
    scales : dict of str to float
        The fitted scale of each dataset (1 where not fitted).
    chi2 : float
        Chi-squared of the data residuals at the solution.
    dof : int
        Degrees of freedom: data points minus parameters.
    dataset_chi2 : dict of str to tuple
        ``(number of points, chi-squared)`` per dataset.
    """

    config: FitConfig
    model: OscillatorModel
    parameters: Float[np.ndarray, " q"]
    covariance: Float[np.ndarray, "q q"] = field(repr=False)
    scales: Mapping[str, float] = field(default_factory=dict)
    chi2: float = 0.0
    dof: int = 0
    dataset_chi2: Mapping[str, tuple[int, float]] = field(default_factory=dict)

    @property
    def reduced_chi2(self) -> float:
        """Chi-squared per degree of freedom.

        Returns
        -------
        float
            ``chi2 / dof`` (``dof`` at least 1).
        """
        return self.chi2 / max(self.dof, 1)

    def polarizability_sigma_au(
        self, omega: Float[np.ndarray, " n"], *, birge: bool = True
    ) -> Float[np.ndarray, " n"]:
        """Return the one-sigma uncertainty of the fitted polarizability.

        Parameters
        ----------
        omega : numpy.ndarray
            Photon energies, hartree.
        birge : bool, optional
            Scale the covariance by the reduced chi-squared when it exceeds one.

        Returns
        -------
        numpy.ndarray
            One-sigma uncertainty, atomic units.
        """
        n_terms = len(self.model.strengths) - len(self.config.fixed_terms)
        strengths = np.asarray(self.model.strengths)[:n_terms]
        energies = np.asarray(self.model.energies)[:n_terms]
        denominator = energies[:, None] ** 2 - omega**2
        d_strength = 1.0 / denominator
        d_energy = -2.0 * strengths[:, None] * energies[:, None] / denominator**2
        n_lines = len(self.config.lines)
        gradient = np.vstack([d_strength, d_energy[n_lines:n_terms]])
        covariance = self.covariance
        if birge and self.reduced_chi2 > 1:
            covariance = covariance * self.reduced_chi2
        return np.sqrt(np.einsum("in,ij,jn->n", gradient, covariance, gradient))


def _initial(config: FitConfig, layout: _Layout, start: int) -> np.ndarray:
    """Return a starting vector: line priors, effective terms spread in energy.

    Parameters
    ----------
    config : FitConfig
        The fit.
    layout : _Layout
        The parameter layout.
    start : int
        Which of the spread starting points to use.

    Returns
    -------
    numpy.ndarray
        The starting parameters.
    """
    f_lines = [line.f_prior for line in config.lines]
    excess = np.geomspace(0.05, 3.0, layout.n_effective) * (1.0 + 0.7 * start)
    f_eff = np.full(layout.n_effective, 1.0)
    return np.concatenate(
        [f_lines, np.log(f_eff), np.log(excess), np.zeros(len(layout.scaled))]
    )


def _warm_start(
    previous: FitResult,
    config: FitConfig,
    layout: _Layout,
    data: Sequence[Observation],
) -> np.ndarray:
    """Return the parameter vector of a previous solution in this fit's layout.

    Parameters
    ----------
    previous : FitResult
        The previous fit.
    config : FitConfig
        This fit.
    layout : _Layout
        This fit's parameter layout.
    data : sequence of Observation
        This fit's observations (their scales are taken from the previous fit where it
        had them).

    Returns
    -------
    numpy.ndarray
        The starting parameters.

    Raises
    ------
    ValueError
        If the previous model has a different number of lines or effective terms.
    """
    # The fixed (infrared) terms come last and are not parameters.
    n_fitted = len(previous.model.strengths) - len(previous.config.fixed_terms)
    strengths = np.asarray(previous.model.strengths)[:n_fitted]
    energies = np.asarray(previous.model.energies)[:n_fitted]
    if len(strengths) != layout.n_lines + layout.n_effective:
        msg = "the previous fit has a different model layout"
        raise ValueError(msg)
    n = layout.n_lines
    # A starting excess energy must be positive for its logarithm.
    excess = np.maximum(energies[n:] - config.effective_min_energy, 1e-12)
    offsets = [previous.scales.get(data[i].key, 1.0) - 1.0 for i in layout.scaled]
    return np.concatenate(
        [strengths[:n], np.log(strengths[n:]), np.log(excess), offsets]
    )


def fit(
    config: FitConfig,
    data: Sequence[Observation] | None = None,
    *,
    starts: int | None = None,
    start: FitResult | None = None,
) -> FitResult:
    """Fit a configuration's model to its data.

    Several starting points spread the effective oscillators over different energies;
    the solution with the lowest cost is kept.

    Parameters
    ----------
    config : FitConfig
        The fit.
    data : sequence of Observation, optional
        The observations; by default those the configuration names
        (:func:`observations`).
    starts : int, optional
        The number of spread starting points; by default the configuration's.
    start : FitResult, optional
        A previous fit of the same model (for instance on more data), whose solution is
        tried first; with ``starts=0`` the fit follows it to the nearest minimum.

    Returns
    -------
    FitResult
        The fitted model, its covariance and statistics.

    Raises
    ------
    ValueError
        If there is no starting point, or ``start`` has a different model layout.
    """
    data = tuple(observations(config) if data is None else data)
    scaled = tuple(
        i for i, obs in enumerate(data) if obs.options.free_scale or obs.systematic
    )
    layout = _Layout(len(config.lines), config.effective, scaled)
    initials = [
        _initial(config, layout, i)
        for i in range(config.starts if starts is None else starts)
    ]
    if start is not None:
        initials.insert(0, _warm_start(start, config, layout, data))
    if not initials:
        msg = "fit needs at least one starting point"
        raise ValueError(msg)
    best = None
    for initial in initials:
        solution = least_squares(
            _residuals,
            initial,
            args=(config, layout, data),
            x_scale="jac",
            method="trf",
        )
        if best is None or solution.cost < best.cost:
            best = solution
    assert best is not None
    strengths, energies, offsets = _unpack(best.x, config, layout)
    # Covariance in the physical parameters (strengths, then effective energies),
    # from the Jacobian in the fitted ones through the chain rule of the logarithms.
    # The dense solver (method "trf" on a dense problem) returns a dense Jacobian.
    jacobian = np.asarray(best.jac, dtype=float)
    covariance_theta = np.linalg.pinv(jacobian.T @ jacobian)
    n_lines, n_eff = layout.n_lines, layout.n_effective
    chain = np.concatenate(
        [
            np.ones(n_lines),
            strengths[n_lines:],
            energies[n_lines:] - config.effective_min_energy,
        ]
    )
    n_physical = n_lines + 2 * n_eff
    covariance = covariance_theta[:n_physical, :n_physical] * np.outer(chain, chain)
    scales = {obs.key: 1.0 for obs in data}
    for index, obs_index in enumerate(layout.scaled):
        scales[data[obs_index].key] = 1.0 + float(offsets[index])
    dataset_chi2 = {}
    chi2 = 0.0
    for obs in data:
        model = _model_values(strengths, energies, obs, config) * scales[obs.key]
        value = float(np.sum(((model - obs.alpha) / obs.sigma) ** 2))
        dataset_chi2[obs.key] = (len(obs.alpha), value)
        chi2 += value
    n_points = sum(len(obs.alpha) for obs in data)
    labels = [line.label for line in config.lines] + [
        f"effective {i + 1}" for i in range(n_eff)
    ]
    model = OscillatorModel(
        gas=config.gas,
        name=f"fit/{config.name}",
        citation=f"fitted here ({config.name}.toml, sha256 {config.digest[:12]})",
        wavelength_range=(config.range_nm[0] * 1e-9, config.range_nm[1] * 1e-9),
        kind="fitted",
        strengths=(
            *(float(f) for f in strengths),
            *(t[1] for t in config.fixed_terms),
        ),
        energies=(*(float(e) for e in energies), *(t[2] for t in config.fixed_terms)),
        labels=(*labels, *(t[0] for t in config.fixed_terms)),
        offset_au=config.chi0,
        retardation_au=config.retardation,
        absorption=config.absorption,
    )
    return FitResult(
        config=config,
        model=model,
        parameters=np.concatenate([strengths, energies[n_lines:]]),
        covariance=covariance,
        scales=scales,
        chi2=chi2,
        dof=n_points - layout.size,
        dataset_chi2=dataset_chi2,
    )


def _validator_observations(config: FitConfig) -> dict[str, Observation | str]:
    """Return the observations of a configuration's validators.

    Parameters
    ----------
    config : FitConfig
        The fit.

    Returns
    -------
    dict of str to Observation or str
        Per validator, its observation at its stated uncertainty alone, or the reason
        it has none (a quantity with no comparison yet).
    """
    by_key = {d.key: d for d in datasets()}
    found: dict[str, Observation | str] = {}
    for validator in config.validators:
        key = validator.key
        dataset = by_key[key]
        options = DatasetOptions(
            key=key,
            extra_relative=1e-300,
            range_nm=validator.range_nm,
            electric=validator.electric,
            infrared=validator.infrared,
            temperature=validator.temperature,
        )
        try:
            if dataset.quantity == "phase_match_density_ratio":
                found[key] = _ratio(dataset, options, config.gas)
            else:
                found[key] = _observation(
                    dataset, options, config.reference_temperature
                )
        except DataError as error:
            found[key] = str(error)
    return found


def validation(result: FitResult) -> dict[str, dict[str, float | int | str]]:
    """Compare the fit with its validator datasets.

    Parameters
    ----------
    result : FitResult
        The fit.

    Returns
    -------
    dict
        Per validator: the number of values, their chi-squared against the fit (where
        uncertainties are stated) and their mean relative deviation from it; or a note
        when the dataset's quantity has no comparison yet. A validator may take a
        wavelength window (``range_nm``); one that is also fitted, in another window,
        is compared at its fitted scale.
    """
    report: dict[str, dict[str, float | int | str]] = {}
    for key, obs in _validator_observations(result.config).items():
        if isinstance(obs, str):
            report[key] = {"note": obs}
            continue
        n_fitted = len(result.model.strengths) - len(result.config.fixed_terms)
        # A validator that is also a fitted dataset (its values beyond the fitted
        # window) is compared at the dataset's fitted scale.
        model = result.scales.get(key, 1.0) * _model_values(
            np.asarray(result.model.strengths)[:n_fitted],
            np.asarray(result.model.energies)[:n_fitted],
            obs,
            result.config,
        )
        keep = np.ones(len(obs.alpha), dtype=bool)
        if obs.orders is not None:
            # Only the sums refraction determines: S(k) for k <= -2. The model's
            # effective oscillators are not a distribution for S(-1) and above.
            keep = obs.orders <= -2
        model, obs_alpha, obs_sigma = model[keep], obs.alpha[keep], obs.sigma[keep]
        deviation = obs_alpha / model - 1
        stated = obs_sigma > 1e-200
        residual = (model - obs_alpha)[stated] / obs_sigma[stated]
        report[key] = {
            "n": int(np.count_nonzero(keep)),
            "chi2": float(np.sum(residual**2)),
            "mean_relative_deviation": float(np.mean(deviation)),
        }
    return report


def _absorption_record(absorption: MeasuredAbsorption | None) -> object:
    """Return what a stored fit records of its measured-absorption term.

    Parameters
    ----------
    absorption : MeasuredAbsorption or None
        The term.

    Returns
    -------
    dict or None
        Its smoothing, grid step, temperature and pieces, or None.
    """
    if absorption is None:
        return None
    return {
        "fwhm_cm": absorption.fwhm_cm,
        "step_cm": absorption.step_cm,
        "temperature_k": absorption.temperature_k,
        "pieces": [[p.source, p.from_nm, p.to_nm] for p in absorption.pieces],
    }


def result_record(result: FitResult) -> dict[str, object]:
    """Return a fit result as a JSON-serialisable record.

    Parameters
    ----------
    result : FitResult
        The result.

    Returns
    -------
    dict
        The model's terms, covariance, statistics and provenance.
    """
    return {
        "name": result.config.name,
        "gas": result.config.gas.value,
        "config_sha256": result.config.digest,
        "range_nm": list(result.config.range_nm),
        "labels": list(result.model.labels),
        "strengths": list(result.model.strengths),
        "energies_hartree": list(result.model.energies),
        "chi0_au": result.config.chi0,
        "retardation_au": result.config.retardation,
        "reference_temperature_k": result.config.reference_temperature,
        "n_fixed": len(result.config.fixed_terms),
        "infrared_max_wavenumber": result.config.infrared_max_wavenumber,
        "absorption": _absorption_record(result.config.absorption),
        "covariance": result.covariance.tolist(),
        "scales": dict(result.scales),
        "chi2": result.chi2,
        "dof": result.dof,
        "dataset_chi2": {k: list(v) for k, v in result.dataset_chi2.items()},
        "validation": validation(result),
    } | _derived_record(result.config)


def _derived_record(config: FitConfig) -> dict[str, object]:
    """Return the record fields of a fit that reads the derived tables (OEP 3).

    Parameters
    ----------
    config : FitConfig
        The fit.

    Returns
    -------
    dict
        ``band_nodes``, ``windows``, ``intervals_nm``, ``gap_reasons``,
        ``n_rotational`` and ``derived_sha256``; empty for a fit that reads neither,
        whose record is unchanged.
    """
    if not config.derived_sha256:
        return {}
    return {
        "band_nodes": config.band_nodes,
        "windows": config.windows,
        "intervals_nm": [list(i) for i in config.intervals_nm],
        "gap_reasons": list(config.gap_reasons),
        "n_rotational": len(config.rotational),
        "derived_sha256": config.derived_sha256,
    }


def dumps(result: FitResult) -> str:
    """Serialise a fit result as JSON text.

    Parameters
    ----------
    result : FitResult
        The result.

    Returns
    -------
    str
        Indented JSON.
    """
    return json.dumps(result_record(result), indent=1) + "\n"

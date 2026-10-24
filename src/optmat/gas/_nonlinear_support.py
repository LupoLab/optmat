# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""The derivation of the hyperpolarizability support records (OEP 5).

A record states where a gas's fit is supported and the gas's resonance limits
(:class:`~optmat.gas.nonlinear.HyperpolarizabilitySupport`). It is derived from:

- the stored fit and its configuration: every fitted point, after the fit's windows and
  row selections, written as four signed legs :math:`(-\nu_\sigma, \nu_1, \nu_2,
  \nu_3)` that sum to zero; the largest :math:`W_2`, single leg and pair sum, with the
  points that attain them;
- ``data/gas/nonlinear_support.toml``: the datasets of each resonance limit, why each
  was chosen, and why any is missing.

``tools/fit/nonlinear_support_record.py`` writes the records, and a test checks that
they are fresh. A fitted dataset that cannot be read, or a recomputed extent or point
count that disagrees with the stored fit, is an error, never a smaller support.
"""

import functools
import math
import tomllib
from collections.abc import Mapping
from typing import Any, Final

import numpy as np

from optmat import _catalogue
from optmat.constants import HBAR, C
from optmat.datasets import Dataset
from optmat.exceptions import DataError
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.gases import Gas, as_gas
from optmat.gas.ionisation import _energy, _row
from optmat.gas.nonlinear import (
    SUPPORT_DERIVATION_VERSION,
    AbsorptionEvidence,
    HyperpolarizabilitySupport,
    SupportCandidate,
    SupportLimit,
    SupportPoint,
    SupportUncertainty,
    _conversion,
    _kept,
    gamma_observation,
    hyperpolarizability_models,
    load_nonlinear_config,
)
from optmat.identity import sha256_digest

__all__ = ["SELECTION_FILE", "derive_support", "fitted_points", "signed_legs"]

SELECTION_FILE: Final = GAS_DATA_DIRECTORY / "nonlinear_support.toml"
"""The per-gas choices of the resonance limits, each with its reasons (TC-870)."""

# Points whose extent agrees with the largest to this relative tolerance are its
# witnesses too: the same wavelength reached by different arithmetic (a standard-air
# conversion, 2 nu against nu + nu) differs only by rounding, about 1e-16.
_TIE: Final = 1e-12
_HC: Final = 2 * math.pi * HBAR * C  # J m: an energy over hc, m^-1
_LIMIT_KEYS: Final = frozenset({"kind", "candidates", "registers", "qualifications"})
_EVIDENCE_KEYS: Final = frozenset(
    {"kind", "dataset", "at", "detected", "registers", "qualifications"}
)
_GAS_KEYS: Final = frozenset(
    {
        "ionisation",
        "ionisation_registers",
        "ionisation_qualifications",
        "e1",
        "e2",
        "e2_any",
        "missing",
    }
)
_ORIGIN: Final = (
    "The vibrationless origin, from v = 0 of the ground state: it does not bound "
    "hot-band lines, and the band system extends above it."
)


def signed_legs(
    process: str, nu: np.ndarray, other: np.ndarray | None = None
) -> np.ndarray:
    r"""Return the signed legs :math:`(-\nu_\sigma, \nu_1, \nu_2, \nu_3)` of a process.

    The four sum to zero; the sum of their squares and of their fourth powers are
    :func:`~optmat.gas.nonlinear.frequency_sums`'s ``W2`` and ``W4``.

    Parameters
    ----------
    process : str
        ``static``, ``kerr``, ``eshg``, ``dc_kerr``, ``thg``, ``cars`` or
        ``pump_probe``.
    nu : numpy.ndarray
        The fundamental (pump) wavenumbers, in any unit.
    other : numpy.ndarray, optional
        The Stokes (``cars``) or probe (``pump_probe``) wavenumbers, same unit.

    Returns
    -------
    numpy.ndarray
        Shape ``(4, n)``, the output first, same unit.

    Raises
    ------
    DataError
        For an unknown process, or ``cars`` and ``pump_probe`` without ``other``.

    Examples
    --------
    >>> signed_legs("thg", np.array([1.0]))[:, 0].tolist()
    [-3.0, 1.0, 1.0, 1.0]
    """
    nu = np.asarray(nu, dtype=float)
    zero = np.zeros_like(nu)
    if process in {"cars", "pump_probe"} and other is None:
        msg = f"{process} needs a second frequency"
        raise DataError(msg)
    table = {
        "static": (zero, zero, zero, zero),
        "kerr": (-nu, nu, nu, -nu),
        "thg": (-3 * nu, nu, nu, nu),
        "eshg": (-2 * nu, nu, nu, zero),
        "dc_kerr": (-nu, nu, zero, zero),
    }
    if process == "cars":
        stokes = np.asarray(other, dtype=float)
        legs = (-(2 * nu - stokes), nu, nu, -stokes)
    elif process == "pump_probe":
        probe = np.asarray(other, dtype=float)
        legs = (-probe, probe, nu, -nu)
    elif process in table:
        legs = table[process]
    else:
        msg = f"unknown third-order process {process!r}"
        raise DataError(msg)
    return np.vstack([np.broadcast_to(leg, nu.shape) for leg in legs])


def _pair_max(legs: tuple[float, ...]) -> float:
    r"""Return the largest :math:`|\nu_a + \nu_b|` over a quartet's six pairs.

    Parameters
    ----------
    legs : tuple of float
        The four signed legs.

    Returns
    -------
    float
        The largest pair sum, same unit.
    """
    return max(abs(legs[a] + legs[b]) for a in range(4) for b in range(a + 1, 4))


def fitted_points(gas: Gas | str) -> list[tuple[SupportPoint, float]]:
    """Return every fitted point of a gas's stored fit, with its W2.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    list of tuple
        Each point (SI legs, m^-1) and its :math:`W_2` (m^-2), in the order of the
        fit's datasets and their rows.

    Raises
    ------
    DataError
        If the stored fit is not of the current configuration, a fitted dataset cannot
        be converted, or the points disagree with the stored fit's counts or its
        largest W2.
    """
    gas = as_gas(gas)
    model = hyperpolarizability_models()[gas]
    config = load_nonlinear_config(gas.value)
    if config.digest != model.record["config_sha256"]:
        msg = f"{gas.value}: the stored fit is not of the current configuration"
        raise DataError(msg)
    counts = model.record["dataset_chi2"]
    assert isinstance(counts, Mapping)
    by_key = _catalogue.by_key()
    out: list[tuple[SupportPoint, float]] = []
    largest = 0.0
    for options in config.fitted:
        dataset = by_key[options.key]
        observation = gamma_observation(dataset, options)
        metres = dataset.vacuum_wavelength()
        wavelength = metres * 1e9
        nu = np.where(np.isfinite(wavelength), 1e7 / wavelength, 0.0)
        _, process, other = _conversion(dataset, nu)
        keep = _kept(dataset, options, wavelength)
        legs = signed_legs(process, nu, other)[:, keep]
        if int(np.count_nonzero(keep)) != int(counts[options.key][0]):  # pyright: ignore[reportIndexIssue]
            msg = f"{options.key}: the fitted points disagree with the stored fit"
            raise DataError(msg)
        # The same squares as frequency_sums', summed in another order: they agree to
        # rounding, about 1e-16 relative.
        if not np.allclose(np.sum(legs**2, axis=0), observation.w2, rtol=1e-12, atol=0):
            msg = f"{options.key}: the legs do not give the fit's W2"
            raise DataError(msg)
        rows = np.array(dataset.rows)[keep]
        seconds = None if other is None else other[keep]
        for i in range(legs.shape[1]):
            point = SupportPoint(
                dataset=options.key,
                row=str(rows[i]),
                process=process,
                wavelength=float(metres[keep][i]),
                second_wavelength=None if seconds is None else 1e-2 / seconds[i],
                legs=tuple(100.0 * legs[:, i]),  # pyright: ignore[reportArgumentType]
            )
            out.append((point, float(observation.w2[i]) * 1e4))
        largest = max(largest, float(np.max(observation.w2)))
    if largest != model.record["max_w2_fitted_cm-2"]:
        msg = f"{gas.value}: the fitted points' largest W2 is not the stored fit's"
        raise DataError(msg)
    return out


def _extent(
    points: list[SupportPoint], values: list[float]
) -> tuple[float, tuple[SupportPoint, ...]]:
    """Return the largest of a quantity over the points, and the points attaining it.

    Parameters
    ----------
    points : list of SupportPoint
        The fitted points.
    values : list of float
        The quantity at each.

    Returns
    -------
    tuple
        The largest value and its witnesses, in order.
    """
    top = max(values)
    witnesses = tuple(
        p for p, v in zip(points, values, strict=True) if v >= top * (1 - _TIE)
    )
    return top, witnesses


def _fitted_extent(
    points: list[SupportPoint], values: list[float], what: str
) -> SupportLimit:
    """Return a fitted extent as a limit.

    Parameters
    ----------
    points : list of SupportPoint
        The fitted points.
    values : list of float
        The quantity at each, m^-1.
    what : str
        What the quantity is.

    Returns
    -------
    SupportLimit
        Kind ``fitted_extent``, with its witnesses.
    """
    top, witnesses = _extent(points, values)
    return SupportLimit(
        wavenumber=top,
        kind="fitted_extent",
        datasets=tuple(dict.fromkeys(p.dataset for p in witnesses)),
        location="the witnesses' rows",
        derivation=(
            f"the largest {what} over the {len(points)} fitted points, after the fit's "
            "windows and row selections; the legs are inverse vacuum wavelengths"
        ),
        registers=("TC-864",),
        witnesses=witnesses,
    )


@functools.cache
def _selection() -> Mapping[str, Any]:
    """Return the parsed selection file.

    Returns
    -------
    Mapping
        Its tables, by gas.
    """
    return tomllib.loads(SELECTION_FILE.read_text(encoding="utf-8"))


def _dataset(by_key: Mapping[str, Dataset], key: str, quantity: str) -> Dataset:
    """Return a selected dataset, checking its quantity.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    key : str
        ``source/dataset``.
    quantity : str
        The quantity it must be.

    Returns
    -------
    Dataset
        The dataset.

    Raises
    ------
    DataError
        If it does not exist or is of another quantity.
    """
    dataset = by_key.get(key)
    if dataset is None or dataset.quantity != quantity:
        msg = f"nonlinear_support.toml: {key!r} is not a {quantity} dataset"
        raise DataError(msg)
    return dataset


def _level(
    by_key: Mapping[str, Dataset], key: str
) -> tuple[float, str, SupportUncertainty]:
    """Return a ``level_energy`` dataset's level, where it is and its uncertainty.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    key : str
        ``source/dataset``, one level in cm^-1.

    Returns
    -------
    tuple
        The level (m^-1), its location and its uncertainty.

    Raises
    ------
    DataError
        If the dataset is not one level in cm^-1.
    """
    dataset = _dataset(by_key, key, "level_energy")
    if dataset.unit != "cm-1" or len(dataset.values) != 1:
        msg = f"nonlinear_support.toml: {key} is not one level in cm-1"
        raise DataError(msg)
    p = dataset.parameters
    location = (
        f"{key}: row {dataset.rows[0]}, {p.get('configuration')} {p.get('term')} "
        f"J = {p.get('j')}, the level energy"
    )
    sigma = float(dataset.standard_uncertainty()[0])
    stated = math.isfinite(sigma)
    uncertainty = SupportUncertainty(
        stated=stated,
        value=100.0 * sigma if stated else None,
        level=dataset.uncertainty.level,
        provenance=(
            f"{key}: the source's stated uncertainty"
            if stated
            else f"{key}: the source states none (its printed precision is not one)"
        ),
    )
    return 100.0 * float(dataset.scaled_values()[0]), location, uncertainty


def _origin(
    by_key: Mapping[str, Dataset], chain: str
) -> tuple[float, str, SupportUncertainty]:
    """Return a molecular state's origin: the sum of a chain of 0-0 bands.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    chain : str
        ``a + b + ...``, ``electronic_state_energy`` datasets whose ``nu00_cm1`` sum to
        the state's origin above the ground state's v = 0.

    Returns
    -------
    tuple
        The origin (m^-1), its location and its (unstated) uncertainty.

    Raises
    ------
    DataError
        If a dataset of the chain gives no 0-0 band.
    """
    keys = [k.strip() for k in chain.split("+")]
    total, parts = 0.0, []
    for key in keys:
        p = _dataset(by_key, key, "electronic_state_energy").parameters
        if "nu00_cm1" not in p:
            msg = f"nonlinear_support.toml: {key} gives no 0-0 band (nu00_cm1)"
            raise DataError(msg)
        total += float(p["nu00_cm1"])  # pyright: ignore[reportArgumentType]
        parts.append(f"{key}: nu00_cm1 of {p.get('state')} ({p.get('transition')})")
    uncertainty = SupportUncertainty(
        stated=False,
        value=None,
        level="not stated",
        provenance=(
            "no uncertainty is stated for the 0-0 bands; the origin is not the "
            "potential minimum Te, and their difference is not an uncertainty"
        ),
    )
    return 100.0 * total, "; ".join(parts), uncertainty


def _strings(entry: Mapping[str, Any], name: str) -> tuple[str, ...]:
    """Return a list of strings from a table.

    Parameters
    ----------
    entry : Mapping
        The table.
    name : str
        The key.

    Returns
    -------
    tuple of str
        The strings, empty where the key is absent.
    """
    return tuple(str(x) for x in entry.get(name, ()))


def _resonance(
    by_key: Mapping[str, Dataset], entry: Mapping[str, Any], where: str
) -> SupportLimit:
    """Return a level or band-origin limit: the smallest of its candidates.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    entry : Mapping
        The limit's table in the selection.
    where : str
        For error messages.

    Returns
    -------
    SupportLimit
        The limit.

    Raises
    ------
    DataError
        For unknown keys or no candidate.
    """
    if set(entry) - _LIMIT_KEYS or not entry.get("candidates"):
        msg = (
            f"nonlinear_support.toml: {where} needs candidates and only "
            f"{sorted(_LIMIT_KEYS)}"
        )
        raise DataError(msg)
    kind = entry["kind"]
    read = _level if kind == "level" else _origin
    found = [(text, *read(by_key, str(text))) for text in entry["candidates"]]
    text, wavenumber, location, uncertainty = min(found, key=lambda c: c[1])
    keys = tuple(k.strip() for k in str(text).split("+"))
    several = len(found) > 1
    qualifications = _strings(entry, "qualifications")
    return SupportLimit(
        wavenumber=wavenumber,
        kind=kind,
        datasets=keys,
        location=location,
        derivation=(
            ("the smallest of the candidates; " if several else "")
            + (
                "the level energy in cm^-1, times 100"
                if kind == "level"
                else "the sum of the chain's 0-0 bands in cm^-1, times 100"
            )
        ),
        registers=_strings(entry, "registers"),
        candidates=tuple(
            SupportCandidate(
                datasets=tuple(k.strip() for k in str(c[0]).split("+")),
                wavenumber=c[1],
            )
            for c in found
        )
        if several
        else (),
        uncertainty=uncertainty,
        qualifications=qualifications + ((_ORIGIN,) if kind == "band_origin" else ()),
    )


def _evidence(
    by_key: Mapping[str, Dataset], entry: Mapping[str, Any], where: str
) -> SupportLimit:
    """Return absorption evidence: one end of a measurement's sampled window.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    entry : Mapping
        The limit's table in the selection.
    where : str
        For error messages.

    Returns
    -------
    SupportLimit
        Kind ``absorption_evidence``.

    Raises
    ------
    DataError
        For unknown keys, or an end that is neither ``lowest`` nor ``highest``.
    """
    if set(entry) != _EVIDENCE_KEYS or entry["at"] not in {"lowest", "highest"}:
        msg = (
            f"nonlinear_support.toml: {where} needs exactly {sorted(_EVIDENCE_KEYS)}, "
            "with at = 'lowest' or 'highest'"
        )
        raise DataError(msg)
    key = str(entry["dataset"])
    dataset = by_key.get(key)
    if dataset is None:
        msg = f"nonlinear_support.toml: no dataset {key!r}"
        raise DataError(msg)
    wavenumbers = 1.0 / dataset.vacuum_wavelength()
    low, high = float(np.min(wavenumbers)), float(np.max(wavenumbers))
    highest = entry["at"] == "highest"
    detected = bool(entry["detected"])
    end = "shortest" if highest else "longest"
    return SupportLimit(
        wavenumber=high if highest else low,
        kind="absorption_evidence",
        datasets=(key,),
        location=f"{key}: the {end} wavelength of its {len(wavenumbers)} values",
        derivation="the inverse vacuum wavelength of that value, m^-1",
        registers=_strings(entry, "registers"),
        evidence=AbsorptionEvidence(
            detected=detected,
            sampled=(low, high),
            qualification=(
                f"absorption {'detected' if detected else 'not detected'} at the "
                f"curated values, {1e9 / high:.1f}-{1e9 / low:.1f} nm; a finite "
                "window, not a level or a bound on one"
            ),
        ),
        qualifications=_strings(entry, "qualifications"),
    )


def _ionisation(
    by_key: Mapping[str, Dataset], entry: Mapping[str, Any]
) -> SupportLimit:
    """Return the ionisation energy over hc as a limit.

    Parameters
    ----------
    by_key : Mapping
        The curated datasets.
    entry : Mapping
        The gas's table in the selection.

    Returns
    -------
    SupportLimit
        Kind ``level``.
    """
    reference = str(entry["ionisation"])
    dataset, index = _row(by_key, reference)
    energy, sigma = _energy(dataset, index)
    stated = math.isfinite(sigma)
    kind = dataset.parameters.get("kind")
    return SupportLimit(
        wavenumber=energy / _HC,
        kind="level",
        datasets=(reference,),
        location=f"{dataset.key}: row {dataset.rows[index]}",
        derivation="the ionisation energy over hc, m^-1",
        registers=_strings(entry, "ionisation_registers"),
        uncertainty=SupportUncertainty(
            stated=stated,
            value=sigma / _HC if stated else None,
            level=dataset.uncertainty.level,
            provenance=(
                f"{dataset.key}: the source's stated uncertainty"
                if stated
                else f"{dataset.key}: the source states none"
            ),
        ),
        qualifications=_strings(entry, "ionisation_qualifications")
        + ((f"The ionisation energy's kind: {kind}.",) if kind else ()),
    )


def derive_support(gas: Gas | str) -> HyperpolarizabilitySupport:
    """Derive a gas's support record from its stored fit, the selection and the data.

    Parameters
    ----------
    gas : Gas or str
        The species.

    Returns
    -------
    HyperpolarizabilitySupport
        The record, status ``validated``.

    Raises
    ------
    DataError
        If the gas has no selection, the selection is malformed, or the fitted points
        cannot be derived (:func:`fitted_points`).
    """
    gas = as_gas(gas)
    entry = _selection().get(gas.value)
    if not isinstance(entry, Mapping) or set(entry) - _GAS_KEYS:
        msg = f"nonlinear_support.toml: no valid entry for {gas.value}"
        raise DataError(msg)
    model = hyperpolarizability_models()[gas]
    pairs = fitted_points(gas)
    points = [p for p, _ in pairs]
    w2 = [v for _, v in pairs]
    max_w2, max_w2_witnesses = _extent(points, w2)
    by_key = _catalogue.by_key()
    limits: dict[str, SupportLimit | None] = {}
    for name in ("e1", "e2", "e2_any"):
        spec = entry.get(name)
        if spec is None:
            limits[name] = None
        elif spec.get("kind") == "absorption_evidence":
            limits[name] = _evidence(by_key, spec, f"{gas.value}.{name}")
        elif spec.get("kind") in {"level", "band_origin"}:
            limits[name] = _resonance(by_key, spec, f"{gas.value}.{name}")
        else:
            msg = f"nonlinear_support.toml: {gas.value}.{name} has no known kind"
            raise DataError(msg)
    ionisation = _ionisation(by_key, entry)
    cited = [k for limit in limits.values() if limit for k in limit.datasets]
    cited += [str(entry["ionisation"]).partition("#")[0]]
    return HyperpolarizabilitySupport(
        gas=gas,
        model_digest=model.identity.digest,
        derivation_version=SUPPORT_DERIVATION_VERSION,
        fitted_points=len(points),
        max_w2=max_w2,
        max_w2_witnesses=max_w2_witnesses,
        leg_max=_fitted_extent(
            points, [max(abs(x) for x in p.legs) for p in points], "single leg |nu_k|"
        ),
        pair_max=_fitted_extent(
            points, [_pair_max(p.legs) for p in points], "pair sum |nu_a + nu_b|"
        ),
        e1=limits["e1"],
        e2=limits["e2"],
        e2_any=limits["e2_any"],
        ionisation=ionisation,
        missing_reasons={str(k): str(v) for k, v in entry.get("missing", {}).items()},
        input_digests={
            "selection": sha256_digest(entry),
            "configuration": str(model.record["config_sha256"]),
            **{key: sha256_digest(by_key[key]) for key in sorted(set(cited))},
        },
        status="validated",
    )

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Where each hyperpolarizability fit is supported, and how a per-leg cap behaves.

    uv run python tools/fit/nonlinear_support.py [--json PATH] [--markdown]
        [--update-docs]

solprop evaluates the fitted gamma = gamma0 [1 + A W2 + B W2^2 + B' W4] as a full
frequency-dependent chi(3), through per-leg spectral filters, so any rule for the
frequencies beyond the fitted data must act on each leg separately. This tool:

1. lists the four signed leg wavenumbers (-nu_o; nu_1, nu_2, nu_3) of every fitted
   value and validator, and the largest single leg among the fitted values (its
   dataset, process and fundamental);
2. reads each gas's resonance limits from its support record (OEP 5, TC-871): the lowest
   dipole-allowed level E1, which a single leg (the output included) reaches, and the
   lowest two-photon level E2, which a pair sum |nu_a + nu_b| reaches, with the
   ionisation energy;
3. checks that the fit increases with each nu_j^2;
4. evaluates leg caps nu_j -> min(|nu_j|, nu_c) at (a) the Kerr diagonal of the largest
   fitted W2, (b) the largest fitted leg, (b') the largest leg any quartet within the
   fitted W2 can have, and (c) the resonance limit min(E1, E2 / 2): their effect on the
   fitted values and on quartets within the fitted W2, whether they are no worse than a
   constant Kerr gamma at a 400, 800 or 1030 nm carrier, the windowed-out calculated
   values beyond the support, and solprop's demonstration quartets.

The manual pages *Third-order nonlinearity* and *Third-order accuracy* tabulate the
results (``--markdown`` prints their tables; ``tests/test_nonlinear_support_page.py``
checks them).
"""

import argparse
import dataclasses
import json
import math
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from optmat.datasets import datasets
from optmat.gas.gases import Gas
from optmat.gas.nonlinear import (
    HyperpolarizabilityModel,
    NonlinearOptions,
    SupportLimit,
    _conversion,
    _dispersion,
    gamma_observation,
    hyperpolarizability_models,
    load_nonlinear_config,
)

GASES = (
    "He", "Ne", "Ar", "Kr", "Xe", "H2", "D2", "N2", "O2", "CO2", "CH4", "SF6", "CO",
    "N2O", "O3",
)  # fmt: skip
NEAR_EDGE = 0.9  # a value is "near the edge" with its largest leg within 10 %
CARRIERS_NM = (400.0, 800.0, 1030.0)
SAMPLES = 20000  # random quartets per gas and test; seeded, so the tables reproduce
SEED = 20261007
OUT = Path("build") / "nonlinear_support"
DOCS = Path(__file__).parents[2] / "docs"
PAGES = (
    DOCS / "user-guide/nonlinear-support.md",
    DOCS / "user-guide/nonlinear-extrapolation.md",
)


# ----------------------------------------------------------------------------------
# Legs


def signed_legs(
    process: str, nu: np.ndarray, other: np.ndarray | None = None
) -> np.ndarray:
    """Return the signed legs (-nu_o, nu_1, nu_2, nu_3) of a process, cm^-1.

    The four sum to zero; W2 is the sum of their squares, as ``frequency_sums``.
    """
    nu = np.asarray(nu, dtype=float)
    zero = np.zeros_like(nu)
    if process == "static":
        legs = (zero, zero, zero, zero)
    elif process == "kerr":
        legs = (-nu, nu, nu, -nu)
    elif process == "thg":
        legs = (-3 * nu, nu, nu, nu)
    elif process == "eshg":
        legs = (-2 * nu, nu, nu, zero)
    elif process == "dc_kerr":
        legs = (-nu, nu, zero, zero)
    elif process == "cars":
        assert other is not None
        legs = (-(2 * nu - other), nu, nu, -other)
    elif process == "pump_probe":
        assert other is not None
        legs = (-other, other, nu, -nu)
    else:
        msg = f"unknown process {process!r}"
        raise ValueError(msg)
    return np.vstack([np.broadcast_to(leg, nu.shape) for leg in legs])


def pair_sums(legs: np.ndarray) -> np.ndarray:
    """Return the largest |nu_a + nu_b| over the six pairs of each quartet, cm^-1.

    A pair and its complement have the same magnitude, the legs summing to zero.
    """
    pairs = [legs[a] + legs[b] for a in range(4) for b in range(a + 1, 4)]
    return np.max(np.abs(np.vstack(pairs)), axis=0)


def gamma(
    model: HyperpolarizabilityModel, legs: np.ndarray, cap: float = math.inf
) -> np.ndarray:
    """Evaluate the polynomial on (capped) legs, au, without the range warnings."""
    magnitude = np.minimum(np.abs(legs), cap)
    return _dispersion(
        model.parameters, np.sum(magnitude**2, axis=0), np.sum(magnitude**4, axis=0)
    )


# ----------------------------------------------------------------------------------
# Step 1: the legs of the data


@dataclass(frozen=True)
class Value:
    """One fitted value or validator with its legs."""

    key: str
    process: str
    wavelength_nm: float
    legs: tuple[float, float, float, float]
    w2: float
    w4: float
    role: str  # fitted, validator, or windowed (a fitted dataset's excluded rows)
    gamma: float = math.nan
    scale: float = 1.0

    @property
    def largest(self) -> float:
        """The largest single-leg wavenumber, cm^-1."""
        return max(abs(leg) for leg in self.legs)


def _values(
    options: Iterable[NonlinearOptions], role: str, scales: dict[str, float]
) -> list[Value]:
    """Convert a configuration's datasets to values with their legs."""
    by_key = {d.key: d for d in datasets()}
    out = []
    for option in options:
        dataset = by_key[option.key]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                obs = gamma_observation(dataset, option)
            except Exception:
                # A validator without a conversion is skipped; a fitted dataset that
                # cannot be read is an error, never a smaller support.
                if role != "validator":
                    raise
                continue
        nu = np.where(np.isfinite(obs.wavelength_nm), 1e7 / obs.wavelength_nm, 0.0)
        _, process, other = _conversion(dataset, nu)
        legs = signed_legs(process, nu, other)
        for i in range(len(nu)):
            out.append(
                Value(
                    key=option.key,
                    process=process,
                    wavelength_nm=float(obs.wavelength_nm[i]),
                    legs=tuple(float(x) for x in legs[:, i]),  # pyright: ignore[reportArgumentType]
                    w2=float(obs.w2[i]),
                    w4=float(obs.w4[i]),
                    role=role,
                    gamma=float(obs.gamma[i]),
                    scale=scales.get(option.key, 1.0),
                )
            )
    return out


def _windowed(options: Iterable[NonlinearOptions]) -> list[NonlinearOptions]:
    """Return each windowed fitted dataset with its excluded fundamentals instead."""
    out = []
    for option in options:
        low, high = option.range_nm
        if low > 0:
            out.append(dataclasses.replace(option, range_nm=(0.0, low * (1 - 1e-9))))
        if math.isfinite(high):
            out.append(
                dataclasses.replace(option, range_nm=(high * (1 + 1e-9), math.inf))
            )
    return out


@dataclass(frozen=True)
class Support:
    """Step 1 for one gas."""

    gas: str
    max_w2: float
    fitted: list[Value]
    validators: list[Value]
    windowed: list[Value]

    @staticmethod
    def _edge(values: Sequence[Value]) -> Value:
        return max(values, key=lambda v: (v.largest, v.w2))

    @property
    def leg_max(self) -> float:
        """nu_leg,max: the largest single leg among the fitted values, cm^-1."""
        return self._edge(self.fitted).largest

    @property
    def pair_max(self) -> float:
        """The largest pair sum |nu_a + nu_b| among the fitted values, cm^-1."""
        return max(float(pair_sums(np.array(v.legs)[:, None])[0]) for v in self.fitted)

    @property
    def leg_max_value(self) -> Value:
        """The fitted value that sets nu_leg,max."""
        return self._edge(self.fitted)

    def near_edge(self, values: Sequence[Value]) -> tuple[int, int, int]:
        """Count the values with their largest leg, or W2, within 10 % of the edge."""
        legs = sum(v.largest >= NEAR_EDGE * self.leg_max for v in values)
        w2 = sum(v.w2 >= NEAR_EDGE * self.max_w2 for v in values)
        return int(legs), int(w2), len(values)


def support(gas: str) -> Support:
    """Collect the legs of a gas's fitted values, validators and windowed rows."""
    model = hyperpolarizability_models()[_gas(gas)]
    config = load_nonlinear_config(gas)
    stored = model.record.get("scales", {})
    assert isinstance(stored, Mapping)
    scales = {str(k): float(v) for k, v in stored.items()}
    return Support(
        gas=gas,
        max_w2=float(model.record["max_w2_fitted_cm-2"]),  # pyright: ignore[reportArgumentType]
        fitted=_values(config.fitted, "fitted", scales),
        validators=_values(config.validators, "validator", scales),
        windowed=_values(_windowed(config.fitted), "windowed", scales),
    )


def _gas(name: str) -> Gas:
    return next(g for g in hyperpolarizability_models() if g.value == name)


# ----------------------------------------------------------------------------------
# Step 2: resonance limits


def _cm1(limit: SupportLimit | None) -> float | None:
    """Return a record's limit in cm^-1 (the tool's unit), or None."""
    return None if limit is None else limit.wavenumber / 100


def _keys(limit: SupportLimit) -> str:
    """Return a limit's datasets (all candidates', where it is the smallest)."""
    if limit.candidates:
        return " + ".join(k for c in limit.candidates for k in c.datasets)
    return " + ".join(limit.datasets)


@dataclass(frozen=True)
class Resonances:
    """Step 2 for one gas, cm^-1, read from its support record (OEP 5, TC-870).

    The selections, their sources and their qualifications are curated in
    ``data/gas/nonlinear_support.toml``; the record is their one copy.
    """

    gas: str
    e1: float | None
    e2: float | None
    e2_any: float | None
    ionisation: float
    e1_kind: str  # level, no_detection, absorption_detected, or missing
    leg_cap: float | None  # min(E1, E2_any / 2): the record's band_edge("any")
    band_edge_parallel: float | None  # min(E1, E2 / 2): its band_edge("linear")
    sources: dict[str, str] = field(default_factory=dict)
    notes: dict[str, str] = field(default_factory=dict)


def resonances(gas: str) -> Resonances:
    """Collect a gas's resonance limits from its support record."""
    record = hyperpolarizability_models()[_gas(gas)].support
    assert record is not None
    limits = {"E1": record.e1, "E2": record.e2, "E2_any": record.e2_any}
    e1 = record.e1
    if e1 is None:
        kind = "missing"
    elif e1.evidence is None:
        kind = "level"
    else:
        kind = "absorption_detected" if e1.evidence.detected else "no_detection"
    sources = {k: _keys(v) for k, v in limits.items() if v is not None}
    sources["IE"] = _keys(record.ionisation)
    names = {"e1": "E1", "e2": "E2", "e2_any": "E2_any"}
    notes = {names[k]: v for k, v in record.missing_reasons.items()}
    if e1 is not None and e1.evidence is not None:
        notes["E1"] = e1.evidence.qualification
    edges = [record.band_edge(p) for p in ("any", "linear")]
    return Resonances(
        gas=gas,
        e1=_cm1(record.e1),
        e2=_cm1(record.e2),
        e2_any=_cm1(record.e2_any),
        ionisation=record.ionisation.wavenumber / 100,
        e1_kind=kind,
        leg_cap=None if edges[0] is None else edges[0] / 100,
        band_edge_parallel=None if edges[1] is None else edges[1] / 100,
        sources=sources,
        notes=notes,
    )


# ----------------------------------------------------------------------------------
# Step 3: monotonicity


def slope_minimum(model: HyperpolarizabilityModel, w2_max: float) -> float:
    """Return the smallest d(gamma/gamma0)/d(nu_j^2), 1e-10 cm^2.

    Over 0 <= nu_j^2 <= W2 <= w2_max.

    The slope is A + 2 B W2 + 2 B' nu_j^2 (in the module's units); scanned on a grid.
    """
    _, a, b, b_prime = model.parameters
    w2 = np.linspace(0.0, w2_max * 1e-10, 201)[:, None]
    fraction = np.linspace(0.0, 1.0, 51)[None, :]
    return float(np.min(a + 2 * b * w2 + 2 * b_prime * fraction * w2))


# ----------------------------------------------------------------------------------
# Step 4: caps


def caps(s: Support, r: Resonances) -> dict[str, float | None]:
    """Return the candidate per-leg caps, cm^-1."""
    return {
        "a": math.sqrt(s.max_w2 / 4),
        "b": s.leg_max,
        "b'": math.sqrt(3 * s.max_w2 / 4),
        "c": r.leg_cap,
        "c1": r.e1,
    }


def _random_quartets(
    rng: np.random.Generator, n: int, low: float, high: float, *, output_free: bool
) -> np.ndarray:
    """Frequency-conserving quartets with inputs of magnitude in [low, high], cm^-1.

    Each input is positive (absorbed) or negative (emitted); the output is their sum,
    kept non-negative, and, unless ``output_free``, also within [low, high].
    """
    chunks: list[np.ndarray] = []
    count = 0
    while count < n:
        size = 4 * n
        magnitude = rng.uniform(low, high, (3, size))
        sign = rng.choice((-1.0, 1.0), (3, size), p=(0.3, 0.7))
        inputs = magnitude * sign
        output = inputs.sum(axis=0)
        keep = output > 0
        if not output_free:
            keep &= (output >= low) & (output <= high)
        chunk = np.vstack([-output[keep], inputs[:, keep]])
        chunks.append(chunk)
        count += chunk.shape[1]
    return np.hstack(chunks)[:, :n]


def _in_w2(rng: np.random.Generator, n: int, w2_max: float) -> np.ndarray:
    """Frequency-conserving quartets with W2 <= w2_max, cm^-1."""
    chunks: list[np.ndarray] = []
    count = 0
    edge = math.sqrt(w2_max)
    while count < n:
        inputs = rng.uniform(-edge, edge, (3, 4 * n))
        legs = np.vstack([-inputs.sum(axis=0), inputs])
        chunk = legs[:, np.sum(legs**2, axis=0) <= w2_max]
        chunks.append(chunk)
        count += chunk.shape[1]
    # The analytic extreme: the third harmonic whose W2 is w2_max has the largest leg
    # any quartet within w2_max can have, sqrt(3 w2_max / 4).
    nu = math.sqrt(w2_max / 12)
    extreme = signed_legs("thg", np.array([nu]))
    return np.hstack([*chunks, extreme])[:, : n + 1]


@dataclass(frozen=True)
class CapResult:
    """Step 4 for one gas and one cap."""

    cap: float | None
    fitted_changed: int
    fitted_total: int
    fitted_worst: float  # relative change
    in_w2_worst: float  # relative change, the most negative
    in_range_violations: dict[float, int]
    in_range_worst_excess: dict[float, float]  # (|cap error| - |carrier error|) / gamma
    beyond: dict[float, dict[str, int]]


def _relative(
    model: HyperpolarizabilityModel, legs: np.ndarray, cap: float
) -> np.ndarray:
    return gamma(model, legs, cap) / gamma(model, legs) - 1


def evaluate_cap(
    model: HyperpolarizabilityModel,
    s: Support,
    r: Resonances,
    cap: float | None,
    rng_seed: int,
) -> CapResult:
    """Evaluate one cap on the data, within the fitted W2 and against the carriers."""
    if cap is None:
        return CapResult(None, 0, len(s.fitted), math.nan, math.nan, {}, {}, {})
    rng = np.random.default_rng(rng_seed)
    data = np.array([v.legs for v in s.fitted]).T
    change = _relative(model, data, cap)
    changed = np.abs(change) > 1e-12
    inside = _in_w2(rng, SAMPLES, s.max_w2)
    in_w2 = _relative(model, inside, cap)
    low_nm, high_nm = (x * 1e9 for x in model.wavelength_range)
    low, high = 1e7 / high_nm, 1e7 / low_nm
    in_range = _random_quartets(rng, SAMPLES, low, high, output_free=False)
    free = _random_quartets(rng, SAMPLES, low, high, output_free=True)
    truth = gamma(model, in_range)
    capped = gamma(model, in_range, cap)
    capped_free = gamma(model, free, cap)
    true_free = gamma(model, free)
    w2_free = np.sum(free**2, axis=0)
    w2_capped = np.sum(np.minimum(np.abs(free), cap) ** 2, axis=0)
    largest = np.max(np.abs(free), axis=0)
    pairs = pair_sums(free)
    resonant = np.zeros(free.shape[1], dtype=bool)
    if r.e1 is not None:
        resonant |= largest >= r.e1
    e2 = r.e2_any if r.e2_any is not None else r.e2
    if e2 is not None:
        resonant |= pairs >= e2
    violations: dict[float, int] = {}
    excess: dict[float, float] = {}
    beyond: dict[float, dict[str, int]] = {}
    for carrier_nm in CARRIERS_NM:
        carrier = float(gamma(model, np.full((4, 1), 1e7 / carrier_nm))[0])
        error_cap = np.abs(capped - truth)
        error_carrier = np.abs(carrier - truth)
        # A rounding allowance: 1e-12 of gamma.
        worse = error_cap > error_carrier + 1e-12 * truth
        violations[carrier_nm] = int(np.sum(worse))
        excess[carrier_nm] = float(np.max((error_cap - error_carrier) / truth))
        supported = w2_free <= s.max_w2
        ok_inside = supported & ~(
            np.abs(capped_free - true_free)
            > np.abs(carrier - true_free) + 1e-12 * true_free
        )
        demonstrated = (~supported) & ~resonant & (carrier <= capped_free)
        beyond[carrier_nm] = {
            "in_support_no_worse": int(np.sum(ok_inside)),
            "in_support_worse": int(np.sum(supported & ~ok_inside)),
            "beyond_demonstrated_cap_in_support": int(
                np.sum(demonstrated & (w2_capped <= s.max_w2))
            ),
            "beyond_demonstrated_cap_beyond_support": int(
                np.sum(demonstrated & (w2_capped > s.max_w2))
            ),
            "beyond_not_demonstrated": int(
                np.sum((~supported) & ~resonant & (carrier > capped_free))
            ),
            "past_resonance": int(np.sum((~supported) & resonant)),
            "total": int(free.shape[1]),
        }
    return CapResult(
        cap=cap,
        fitted_changed=int(np.sum(changed)),
        fitted_total=len(s.fitted),
        fitted_worst=float(change[np.argmax(np.abs(change))]) if len(change) else 0.0,
        in_w2_worst=float(np.min(in_w2)),
        in_range_violations=violations,
        in_range_worst_excess=excess,
        beyond=beyond,
    )


# The demonstration quartets: (name, gas, wavelengths of the inputs in nm, signs).
CASES = (
    ("Ar 400 + 400 - 800 -> 266.7 nm", "Ar", (400.0, 400.0, 800.0), (1, 1, -1)),
    ("Ar SPM at 266.7 nm", "Ar", (800.0 / 3, 800.0 / 3, 800.0 / 3), (1, 1, -1)),
    ("He 250 + 250 - 1500 -> 136.4 nm", "He", (250.0, 250.0, 1500.0), (1, 1, -1)),
    ("He 200 + 200 - 600 -> 120 nm", "He", (200.0, 200.0, 600.0), (1, 1, -1)),
    ("He SPM at 120 nm", "He", (120.0, 120.0, 120.0), (1, 1, -1)),
)


def case_legs(wavelengths: Sequence[float], signs: Sequence[int]) -> np.ndarray:
    """Return the signed legs of a quartet from its inputs, cm^-1."""
    inputs = np.array([s * 1e7 / w for w, s in zip(wavelengths, signs, strict=True)])
    return np.concatenate([[-inputs.sum()], inputs])[:, None]


def case(
    name: str, gas: str, wavelengths: Sequence[float], signs: Sequence[int]
) -> dict[str, object]:
    """Evaluate one demonstration quartet under every cap."""
    model = hyperpolarizability_models()[_gas(gas)]
    s, r = support(gas), resonances(gas)
    legs = case_legs(wavelengths, signs)
    output_nm = 1e7 / abs(float(legs[0, 0]))
    largest = float(np.max(np.abs(legs)))
    pair = float(pair_sums(legs)[0])
    result: dict[str, object] = {
        "name": name,
        "gas": gas,
        "output_nm": output_nm,
        "w2_over_max": float(np.sum(legs**2) / s.max_w2),
        "leg_over_leg_max": largest / s.leg_max,
        "leg_over_e1": largest / r.e1 if r.e1 else math.nan,
        "pair_over_e2": pair / r.e2 if r.e2 else math.nan,
        "pair_over_e2_any": pair / r.e2_any if r.e2_any else math.nan,
        "uncapped_au": float(gamma(model, legs)[0]),
        "carrier_au": {
            nm: float(gamma(model, np.full((4, 1), 1e7 / nm))[0]) for nm in CARRIERS_NM
        },
        "capped_au": {
            k: (float(gamma(model, legs, c)[0]) if c is not None else math.nan)
            for k, c in caps(s, r).items()
        },
    }
    result["past_resonance"] = bool(
        (r.e1 is not None and largest >= r.e1)
        or (r.e2_any is not None and pair >= r.e2_any)
    )
    return result


def beyond_support_evidence(s: Support) -> list[dict[str, float | str]]:
    """Compare the polynomial with the windowed-out and validator values beyond it.

    The model carries the dataset's fitted scale.

    Only values beyond the fitted W2; the ratio is model over datum.
    """
    model = hyperpolarizability_models()[_gas(s.gas)]
    out = []
    for value in [*s.windowed, *s.validators]:
        if value.w2 <= s.max_w2:
            continue
        if not math.isfinite(value.gamma) or value.gamma == 0:
            continue
        legs = np.array(value.legs)[:, None]
        scale = value.scale if value.role == "windowed" else 1.0
        out.append(
            {
                "key": value.key,
                "role": value.role,
                "process": value.process,
                "wavelength_nm": value.wavelength_nm,
                "w2_over_max": value.w2 / s.max_w2,
                "leg_over_leg_max": value.largest / s.leg_max,
                "ratio": float(scale * gamma(model, legs)[0] / value.gamma),
            }
        )
    return out


def derived(s: Support, r: Resonances, range_low_nm: float) -> dict[str, object]:
    """Return the quantities that follow from a gas's facts by definition.

    sqrt(3 W2,max / 4) is the largest leg any quartet within the fitted W2 can have.
    min(E1, E2 / 2) is the lowest frequency whose own SPM quartet reaches a resonance
    (a leg at E1 or a pair sum at E2) for one polarisation, and a sufficient bound for
    any; None where E1 or E2 is not curated or E1 is only absorption evidence. How to
    use them is the propagation code's decision, not optmat's.
    """
    return {
        "leg_max_in_w2": math.sqrt(3 * s.max_w2 / 4),
        "band_edge": r.band_edge_parallel,
        "band_edge_any": r.leg_cap,
        "e1_in_range": r.e1 is not None and r.e1 < 1e7 / range_low_nm,
    }


# ----------------------------------------------------------------------------------
# The whole analysis


def analyse(gas: str) -> dict[str, object]:
    """Run steps 1 to 4 for one gas and return them as plain data."""
    model = hyperpolarizability_models()[_gas(gas)]
    s, r = support(gas), resonances(gas)
    edge = s.leg_max_value
    validators = s._edge(s.validators) if s.validators else None
    if validators is not None and validators.largest == 0:
        validators = None  # static validators only
    candidate = caps(s, r)
    seeds = {name: SEED + i for i, name in enumerate(candidate)}
    return {
        "gas": gas,
        "range_nm": [x * 1e9 for x in model.wavelength_range],
        "max_w2": s.max_w2,
        "kerr_diagonal_nm": 1e7 / math.sqrt(s.max_w2 / 4),
        "leg_max": s.leg_max,
        "leg_max_nm": 1e7 / s.leg_max,
        "leg_max_source": {
            "key": edge.key,
            "process": edge.process,
            "fundamental_nm": edge.wavelength_nm,
            "legs": edge.legs,
            "w2": edge.w2,
            "w4": edge.w4,
        },
        "pair_max": s.pair_max,
        "near_edge_fitted": s.near_edge(s.fitted),
        "validators": None
        if validators is None
        else {
            "leg_max": validators.largest,
            "leg_max_nm": 1e7 / validators.largest,
            "key": validators.key,
            "process": validators.process,
            "fundamental_nm": validators.wavelength_nm,
            "w2": validators.w2,
            "w4": validators.w4,
            "near_edge": s.near_edge(s.validators),
        },
        "fitted_values": [dataclasses.asdict(v) for v in s.fitted],
        "validator_values": [dataclasses.asdict(v) for v in s.validators],
        "resonances": {
            "E1": r.e1,
            "E1_kind": r.e1_kind,
            "E2": r.e2,
            "E2_any": r.e2_any,
            "IE": r.ionisation,
            "leg_cap": r.leg_cap,
            "band_edge_parallel": r.band_edge_parallel,
            "sources": r.sources,
            "notes": r.notes,
        },
        "slope_minimum_in_support": slope_minimum(model, s.max_w2),
        "slope_minimum_to_4x": slope_minimum(model, 4 * s.max_w2),
        "parameters": list(model.parameters),
        "caps": {
            name: dataclasses.asdict(evaluate_cap(model, s, r, cap, seeds[name]))
            for name, cap in candidate.items()
        },
        "beyond_support": beyond_support_evidence(s),
        "derived": derived(s, r, model.wavelength_range[0] * 1e9),
    }


# ----------------------------------------------------------------------------------
# The manual's tables


def _nm(wavenumber: float | None, factor: float = 1.0) -> str:
    """Format a wavenumber (cm^-1) as the wavelength of `factor` photons sharing it."""
    if wavenumber is None:
        return "n/a"
    return f"{factor * 1e7 / wavenumber:.1f}"


def _pct(x: float) -> str:
    return f"{100 * x:+.1f} %" if abs(x) >= 5e-4 else "0"


PROCESS = {
    "kerr": "Kerr",
    "thg": "THG",
    "eshg": "ESHG",
    "dc_kerr": "dc Kerr",
    "cars": "CARS",
    "pump_probe": "pump-probe",
    "static": "static",
}


B_PRIME = "b'"


def _table(header: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _e1_cell(x: dict) -> str:
    """Format E1, or the absorption evidence that stands in for it."""
    e1, kind = x["E1"], x["E1_kind"]
    if e1 is None:
        return "not curated"
    text = {
        "level": f"{e1:.0f} ({_nm(e1)})",
        "no_detection": f"none detected to {_nm(e1)} nm",
        "absorption_detected": f"detected in a curated window ending at {e1:.0f} "
        f"({_nm(e1)})",
    }
    return text[kind]


def _notes(gas_result: dict) -> str:
    """Return the facts that qualify a gas's derived quantities, for the table."""
    x, d = gas_result["resonances"], gas_result["derived"]
    notes = []
    if x["E1"] is None:
        notes.append("E1 not curated")
    elif x["E1_kind"] != "level":
        notes.append("E1: absorption evidence only")
    if d["e1_in_range"]:
        notes.append("absorbs inside the declared range")
    if x["E2"] is None:
        notes.append("E2 not curated")
    return "; ".join(notes) or "-"


def render(result: dict) -> dict[str, str]:
    """Return the manual's marker tables, by name, from the analysis result."""
    gases = [g for g in GASES if g in result]
    tables = {}
    rows = []
    for g in gases:
        r = result[g]
        s = r["leg_max_source"]
        v = r["validators"]
        near = r["near_edge_fitted"]
        rows.append(
            [
                g,
                f"{r['leg_max_nm']:.1f} ({r['leg_max']:.0f})",
                f"`{s['key']}`, {PROCESS[s['process']]} at "
                f"{s['fundamental_nm']:.1f} nm",
                f"{s['w2'] / 1e9:.3f}",
                f"{s['w4'] / 1e18:.3f}",
                f"{near[0]} / {near[1]} / {near[2]}",
                f"{r['max_w2'] / 1e9:.3f}",
                f"{r['kerr_diagonal_nm']:.1f}",
                "none"
                if v is None
                else f"{v['leg_max_nm']:.1f}, `{v['key']}` ({PROCESS[v['process']]})",
            ]
        )
    tables["nonlinear-support"] = _table(
        [
            "Gas",
            "Largest fitted leg, nm (cm⁻¹)",
            "Set by",
            "W2 there (1e9 cm⁻²)",
            "W4 there (1e18 cm⁻⁴)",
            "Near the edge: leg / W2 / values",
            "Largest fitted W2 (1e9 cm⁻²)",
            "Its Kerr diagonal (nm)",
            "Largest validator leg (nm)",
        ],
        rows,
    )
    rows = []
    for g in gases:
        x = result[g]["resonances"]
        sources = "; ".join(
            f"{k}: `{v}`" for k, v in x["sources"].items() if v != "not curated"
        )
        rows.append(
            [
                g,
                _e1_cell(x),
                "not curated"
                if x["E2"] is None
                else f"{x['E2']:.0f} ({_nm(x['E2'], 2)})",
                "not curated"
                if x["E2_any"] is None
                else f"{x['E2_any']:.0f} ({_nm(x['E2_any'], 2)})",
                f"{x['IE']:.0f} ({_nm(x['IE'])}, {_nm(x['IE'], 3)})",
                f"{_nm(x['band_edge_parallel'])} / {_nm(x['leg_cap'])}",
                sources,
            ]
        )
    tables["nonlinear-resonances"] = _table(
        [
            "Gas",
            "E1, cm⁻¹ (one leg, nm)",
            "E2, cm⁻¹ (a degenerate pair, nm)",
            "E2, any polarisation",
            "Ionisation, cm⁻¹ (one leg, three photons, nm)",
            "Band edge min(E1, E2/2), nm: one polarisation / any",
            "Datasets",
        ],
        rows,
    )
    rows = []
    for g in gases:
        c = result[g]["caps"]

        def cell(name: str, c: dict = c) -> str:
            k = c[name]
            if k["cap"] is None:
                return "n/a"
            changed = f"{k['fitted_changed']}/{k['fitted_total']}"
            worst = _pct(k["fitted_worst"]) if k["fitted_changed"] else "exact"
            return f"{_nm(k['cap'])}: {changed}, {worst}; {_pct(k['in_w2_worst'])}"

        rows.append([g, cell("a"), cell("b"), cell("b'"), cell("c")])
    tables["nonlinear-caps"] = _table(
        [
            "Gas",
            "(a) Kerr diagonal",
            "(b) largest fitted leg",
            "(b') √(3 W2,max / 4)",
            "(c) min(E1, E2/2)",
        ],
        rows,
    )
    rows = []
    for g in gases:
        c = result[g]["caps"]
        row = [g]
        for name in ("a", "b", "b'"):
            b = c[name]["beyond"]["400.0"]
            row.append(f"{b['in_support_worse']} / {b['beyond_not_demonstrated']}")
        b = c["b'"]["beyond"]["400.0"]
        row += [
            str(
                b["beyond_demonstrated_cap_in_support"]
                + b["beyond_demonstrated_cap_beyond_support"]
            ),
            str(b["past_resonance"]),
            str(
                sum(
                    c["b'"]["beyond"][k]["in_support_worse"]
                    + c["b'"]["beyond"][k]["beyond_not_demonstrated"]
                    for k in ("800.0", "1030.0")
                )
            ),
        ]
        rows.append(row)
    tables["nonlinear-carrier"] = _table(
        [
            "Gas",
            "(a): worse in support / not shown beyond",
            "(b)",
            "(b')",
            "(b'): shown no worse beyond",
            "Past a resonance",
            "(b'), 800 and 1030 nm carriers: worse or not shown",
        ],
        rows,
    )
    rows = []
    for g in gases:
        by_key: dict[str, list[dict]] = {}
        for e in result[g]["beyond_support"]:
            by_key.setdefault(e["key"], []).append(e)
        for key, entries in by_key.items():
            far = max(entries, key=lambda e: e["w2_over_max"])
            near = min(entries, key=lambda e: e["w2_over_max"])
            rows.append(
                [
                    g,
                    f"`{key}`",
                    entries[0]["role"],
                    f"{PROCESS[far['process']]}, {near['wavelength_nm']:.1f}-"
                    f"{far['wavelength_nm']:.1f} nm",
                    f"{near['w2_over_max']:.2f}-{far['w2_over_max']:.2f}",
                    f"{near['ratio']:.3f} → {far['ratio']:.3f}",
                ]
            )
    tables["nonlinear-beyond"] = _table(
        [
            "Gas",
            "Dataset",
            "Role",
            "Process, fundamentals",
            "W2 / W2,max",
            "Model / value, nearest → farthest",
        ],
        rows,
    )
    rows = []
    for case_result in result["cases"]:
        caps_au = case_result["capped_au"]
        carriers = case_result["carrier_au"]
        rows.append(
            [
                case_result["name"],
                f"{case_result['w2_over_max']:.2f}",
                f"{case_result['leg_over_leg_max']:.2f}",
                f"{case_result['leg_over_e1']:.2f}",
                f"{case_result['pair_over_e2']:.3f}",
                f"{case_result['uncapped_au']:.1f}",
                f"{caps_au['a']:.1f}",
                f"{caps_au['b']:.1f}",
                f"{caps_au[B_PRIME]:.1f}",
                " / ".join(f"{carriers[k]:.1f}" for k in ("400.0", "800.0", "1030.0")),
                "yes" if case_result["past_resonance"] else "no",
            ]
        )
    tables["nonlinear-cases"] = _table(
        [
            "Quartet",
            "W2 / W2,max",
            "Largest leg / fitted",
            "Largest leg / E1",
            "Largest pair / E2",
            "Polynomial (au)",
            "(a)",
            "(b)",
            "(b')",
            "Kerr carrier 400 / 800 / 1030 nm",
            "Past a resonance",
        ],
        rows,
    )
    rows = []
    for g in gases:
        d = result[g]["derived"]
        edge = (
            "n/a"
            if d["band_edge"] is None
            else f"{_nm(d['band_edge'])} / {_nm(d['band_edge_any'])}"
        )
        rows.append(
            [
                g,
                f"{_nm(d['leg_max_in_w2'])} ({d['leg_max_in_w2']:.0f})",
                edge,
                f"{_nm(result[g]['pair_max'], 2)}",
                _notes(result[g]),
            ]
        )
    tables["nonlinear-derived"] = _table(
        [
            "Gas",
            "Largest leg within the fitted W2, √(3 W2,max / 4), nm (cm⁻¹)",
            "min(E1, E2/2), nm: one polarisation / any",
            "Largest fitted pair, as a degenerate fundamental (nm)",
            "Note",
        ],
        rows,
    )
    return tables


def main(argv: Sequence[str] | None = None) -> None:
    """Run the analysis for every gas and write the JSON (and print the tables)."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument("--json", type=Path, default=OUT / "support.json")
    parser.add_argument("--markdown", action="store_true")
    parser.add_argument(
        "--update-docs",
        action="store_true",
        help="rewrite the marker tables of the two manual pages",
    )
    parser.add_argument("gases", nargs="*", default=list(GASES))
    args = parser.parse_args(argv)
    result = {gas: analyse(gas) for gas in args.gases}
    result["cases"] = [case(*c) for c in CASES]  # pyright: ignore[reportArgumentType]
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(f"wrote {args.json}")
    # JSON round trip: the tables read the result as the page tests do.
    tables = render(json.loads(json.dumps(result, default=float)))
    if args.markdown:
        for name, table in tables.items():
            print(f"<!-- {name}:start -->\n{table}\n<!-- {name}:end -->\n")
    if args.update_docs:
        for page in PAGES:
            text = page.read_text(encoding="utf-8")
            for name, table in tables.items():
                start, end = f"<!-- {name}:start -->", f"<!-- {name}:end -->"
                if start in text:
                    head, rest = text.split(start, 1)
                    _, tail = rest.split(end, 1)
                    text = f"{head}{start}\n{table}\n{end}{tail}"
            page.write_text(text, encoding="utf-8")
            print(f"updated {page}")


if __name__ == "__main__":
    main()

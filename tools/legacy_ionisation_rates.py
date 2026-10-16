# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Add solprop's and optmat's ionisation rates to the legacy-rates fixture.

The fixture ``tests/data/legacy_ionisation_rates.json`` compares the instantaneous
ionisation rates of the noble gases at 800 nm, and the yields of Gaussian pulses, in
the two codes earlier propagation used, Luna.jl and solprop, and in optmat. It is
written in three steps, the first of which writes it afresh:

    julia --project=<Luna checkout> tools/legacy_ionisation_rates.jl
    cd <solprop checkout>
    uv run python <optmat checkout>/tools/legacy_ionisation_rates.py solprop
    cd <optmat checkout>
    uv run python tools/legacy_ionisation_rates.py optmat

The first (``tools/legacy_ionisation_rates.jl``) sets the fields and the pulses and
records Luna's rate. This script adds, on the same fields and pulses:

- ``solprop`` (run in *solprop's* environment from its checkout; this package never
  imports solprop): the rate solprop's capillary driver uses with ``plasma=True``,
  ``solprop.driver.simulation._ionisation_rate``, which is ``PPT(gas, wavelength)`` with
  its defaults for the noble gases, evaluated directly. The driver evaluates it through
  a spline table converged to 1e-6 (``build_rate_table`` over
  ``default_field_range``); the table's largest deviation from the direct rate on these
  fields is recorded. ``rate_bsi_none`` is the same rate without its barrier-suppression
  correction, the form Luna has;
- ``optmat`` (run in this package's environment): ``ionisation_rate_model(gas, 800e-9,
  model).instantaneous_rate`` for the models ``default`` and ``one_electron``, with the
  fields at which it warns (``ModelRangeWarning``, above ``maximum_field``) and the
  Keldysh parameter of each field. The pulses are ``gaussian_pulse``'s, checked against
  the ones the Julia step used; ``pulse_yield(model, peak, fwhm, carrier=True)`` is
  recorded as a cross-check of the integration here.

Every rate is per atom, in 1/s, at the instantaneous field; every yield is
``1 - exp(-int W dt)`` by the trapezoid rule on the pulse's own grid.
"""

import json
import math
import subprocess
import sys
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "data"
    / "legacy_ionisation_rates.json"
)
GASES = ("He", "Ne", "Ar", "Kr", "Xe")
C = 299792458.0
SAMPLES_PER_CYCLE = 96  # optmat's ionisation_validation._SAMPLES_PER_CYCLE
HALF_WIDTH = 3.0  # FWHM either side of the peak (_PULSE_HALF_WIDTH)

type Rate = Callable[[np.ndarray], np.ndarray]


def _git(directory: Path) -> dict[str, object]:
    def run(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=directory, capture_output=True, text=True, check=True
        ).stdout.strip()

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
        "clean": run("status", "--porcelain") == "",
    }


def _pulse(peak_amplitude: float, fwhm: float, wavelength: float) -> tuple[Any, Any]:
    """Return the times and the field of optmat's ``gaussian_pulse``, CEP 0."""
    omega = 2 * math.pi * C / wavelength
    period = wavelength / C
    half = HALF_WIDTH * fwhm
    count = 2 * math.ceil(half / period * SAMPLES_PER_CYCLE) + 1
    times = np.linspace(-half, half, count)
    shape = np.exp(-2 * math.log(2) * times**2 / fwhm**2)
    slope = -4 * math.log(2) * times / fwhm**2 * shape
    phase = omega * times
    # -dA/dt of A = -(F0 / omega) e(t) sin(phase).
    return times, peak_amplitude * (
        shape * np.cos(phase) + slope / omega * np.sin(phase)
    )


def _yield(rate: Rate, times: np.ndarray, field: np.ndarray) -> float:
    return float(-np.expm1(-np.trapezoid(rate(field), times)))


def _pulses(data: dict[str, Any]) -> list[tuple[str, float, list[float]]]:
    """Return the pulses as (label, fwhm, peak amplitudes), checked against Julia's."""
    pulses = data["pulses"]
    out = []
    for fwhm in pulses["fwhm_s"]:
        label = f"{round(fwhm * 1e15)}fs"
        times, field = _pulse(1.0, fwhm, data["wavelength"])
        check = pulses["checks"][label]
        assert len(times) == check["n"]
        assert math.isclose(times[0], check["t_first"], rel_tol=1e-12)
        assert math.isclose(
            float(np.sum(field**2)), check["sum_field_squared_unit_peak"], rel_tol=1e-10
        )
        out.append((label, fwhm, pulses["peak_amplitude_v_per_m"]))
    return out


def _solprop(data: dict[str, Any]) -> dict[str, object]:
    # solprop is not a dependency of this package: this mode runs in solprop's
    # environment.
    from solprop.driver.simulation import (  # noqa: PLC0415  # pyright: ignore[reportMissingImports]
        _ionisation_rate,  # pyright: ignore[reportPrivateUsage]
        default_field_range,
    )
    from solprop.nonlinear.ionisation import (  # noqa: PLC0415  # pyright: ignore[reportMissingImports]
        PPT,
    )
    from solprop.translate.ionisation import (  # noqa: PLC0415  # pyright: ignore[reportMissingImports]
        build_rate_table,
    )

    wavelength = data["wavelength"]
    fields = np.asarray(data["field_v_per_m"])
    gases: dict[str, object] = {}
    for gas in GASES:
        rate = _ionisation_rate(gas, True, wavelength)
        assert rate == PPT(gas, wavelength), "the driver's rate is not PPT's default"
        table = build_rate_table(rate, default_field_range(float(fields.max())))
        direct = rate(fields)
        gases[gas] = {
            "class": type(rate).__name__,
            "citation": rate.citation,
            "parameters": dict(rate.parameters),
            "table_field_range_v_per_m": [table.field_min, table.field_max],
            "table_max_relative_difference_from_direct": float(
                np.max(np.abs(table(fields) / direct - 1))
            ),
            "rate": direct.tolist(),
            "rate_bsi_none": PPT(gas, wavelength, bsi="none")(fields).tolist(),
            "pulses": {
                label: {
                    "yield": [
                        _yield(rate, *_pulse(peak, fwhm, wavelength)) for peak in peaks
                    ]
                }
                for label, fwhm, peaks in _pulses(data)
            },
        }
    return {
        **_git(Path.cwd()),
        "function": "solprop.driver.simulation._ionisation_rate(gas, True, 800e-9) = "
        "solprop.nonlinear.ionisation.PPT(gas, 800e-9), the capillary driver's "
        "plasma=True rate for the noble gases",
        "options": "defaults: sum_integral=True, bsi='auto' (Tong-Lin, alpha 7 for He "
        "and 9 for the others), cycle_averaged=False, stark_shift=True, "
        "dipole_correction=True (both inactive for Kr and Xe: no polarisability in "
        "solprop), m_sum=True, occupancy=2, coulomb_factor=None (the asymptotic "
        "C_n*l* with l* = n* - 1), maximum 1/au_time",
        "gases": gases,
    }


def _warned(rate: Rate, field: np.ndarray) -> tuple[np.ndarray, bool]:
    from optmat.exceptions import ModelRangeWarning  # noqa: PLC0415

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ModelRangeWarning)
        values = rate(field)
    return values, any(issubclass(w.category, ModelRangeWarning) for w in caught)


def _optmat(data: dict[str, Any]) -> dict[str, object]:
    import optmat  # noqa: PLC0415
    from optmat.gas.ionisation_rates import ionisation_rate_model  # noqa: PLC0415
    from optmat.gas.ionisation_validation import (  # noqa: PLC0415
        gaussian_pulse,
        pulse_yield,
    )

    wavelength = data["wavelength"]
    fields = np.asarray(data["field_v_per_m"])
    pulses = _pulses(data)
    peak_intensity = np.asarray(data["pulses"]["peak_intensity_w_cm2"]) * 1e4
    # optmat's pulse against the one the Julia step used (their eps0 and c differ in
    # the last digits only).
    pulse_difference = 0.0
    for _, fwhm, peaks in pulses:
        for peak, intensity in zip(peaks, peak_intensity, strict=True):
            _, _, field = gaussian_pulse(float(intensity), fwhm, wavelength)
            _, ours = _pulse(peak, fwhm, wavelength)
            pulse_difference = max(
                pulse_difference, float(np.max(np.abs(field - ours)) / peak)
            )
    assert pulse_difference < 1e-9, pulse_difference
    gases: dict[str, object] = {}
    for gas in GASES:
        entry: dict[str, object] = {}
        for name in ("default", "one_electron"):
            model = ionisation_rate_model(gas, wavelength, name)
            rate = model.instantaneous_rate
            warned = [_warned(rate, np.asarray(f))[1] for f in fields]
            values, _ = _warned(rate, fields)
            yields: dict[str, object] = {}
            for label, fwhm, _ in pulses:
                ours, check, pulse_warned = [], [], []
                for intensity in peak_intensity:
                    times, _, field = gaussian_pulse(float(intensity), fwhm, wavelength)
                    rates, flag = _warned(rate, field)
                    ours.append(float(-np.expm1(-np.trapezoid(rates, times))))
                    pulse_warned.append(flag)
                    check.append(
                        pulse_yield(model, float(intensity), fwhm, carrier=True)
                    )
                yields[label] = {
                    "yield": ours,
                    "pulse_yield_carrier": check,
                    "warned": pulse_warned,
                }
            entry[name] = {
                "identity": model.identity.name,
                "digest": model.identity.digest,
                "citation": model.citation,
                "settings": {
                    "coulomb": model.coulomb,
                    "static_correction": model.static_correction,
                    "static_constants": list(model.static_constants),
                    "structure": model.structure,
                    "structure_coefficient": model.structure_coefficient,
                    "electrons": model.electrons,
                    "channels": [
                        [c.label, c.ionization_energy, c.weight] for c in model.channels
                    ],
                    "channel_sum": model.channel_sum,
                    "averaging": model.averaging,
                },
                "maximum_field_v_per_m": model.maximum_field,
                "barrier_suppression_field_v_per_m": model.barrier_suppression_field,
                "keldysh": model.keldysh(fields).tolist(),
                "rate": values.tolist(),
                "warned": warned,
                "pulses": yields,
            }
        gases[gas] = entry
    return {
        **_git(FIXTURE.parents[2]),
        "version": optmat.__version__,
        "function": "optmat.gas.ionisation_rates.ionisation_rate_model(gas, 800e-9, "
        "model).instantaneous_rate(field)",
        "pulse_max_difference_over_peak": pulse_difference,
        "gases": gases,
    }


def main() -> None:
    """Add the section named on the command line (solprop or optmat) to the fixture."""
    mode = sys.argv[1] if len(sys.argv) == 2 else ""
    sections = {"solprop": _solprop, "optmat": _optmat}
    if mode not in sections:
        sys.exit("usage: legacy_ionisation_rates.py solprop|optmat")
    data = json.loads(FIXTURE.read_text())
    if "luna" not in data:
        sys.exit(f"{FIXTURE} has no Luna section: run legacy_ionisation_rates.jl first")
    data[mode] = sections[mode](data)
    FIXTURE.write_text(json.dumps(data, indent=1) + "\n")


if __name__ == "__main__":
    main()

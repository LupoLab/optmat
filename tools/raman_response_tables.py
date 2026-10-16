# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The tables of the page *The Raman response* (``docs/raman_response.md``).

    uv run python tools/raman_response_tables.py

prints the tables; ``tests/test_raman_response_page.py`` checks that the page holds
exactly what this script prints:

- ``response-terms``: each gas's response at 295 K: its lines, damping groups, the
  density and temperature ranges, the anisotropic share of the Q branch;
- ``response-gain``: the gains from the response's transform against the measured;
- ``response-dephasing``: the decay of the rotational revivals against the measured;
- ``response-solprop``: optmat's responses against solprop's
  (``tests/data/solprop_raman_response.json``);
- ``response-rank``: how many separable terms represent a response over a range of
  densities, for a propagation code (the rank of the kernel's family);
- ``silica-models``, ``silica-validation``: fused silica's models and their checks.
"""

import json
import math
import warnings
from pathlib import Path

import numpy as np

from optmat.constants import C
from optmat.exceptions import ModelRangeWarning
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman_response import raman_response
from optmat.gas.raman_response_validation import gas_gains, rotational_dephasing
from optmat.solid.raman_response import silica_raman_responses
from optmat.solid.raman_response_validation import (
    silica_gains,
    silica_low_frequency_gain,
    silica_raman_fractions,
)

ROOT = Path(__file__).resolve().parents[1]
GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O", "dry air")
WAVENUMBER = 2 * math.pi * C * 100
# The rank study: gas, density range (amagat), time step (s), number of samples
# (2^15: a window of 131 ps at 4 fs, 65 ps at 2 fs, as a propagation grid might be).
RANKS = (
    ("N2", (0.5, 2.0), 4e-15),
    ("N2", (0.1, 10.0), 4e-15),
    ("O2", (0.1, 10.0), 4e-15),
    ("dry air", (0.1, 10.0), 4e-15),
    ("H2", (1.0, 30.0), 2e-15),
    ("D2", (1.0, 30.0), 2e-15),
)
RANK_SAMPLES = 2**15
RANK_DENSITIES = 64
TOLERANCES = (1e-3, 1e-4, 1e-6)


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _amagat(value: float) -> str:
    if value == 0:
        return "0"
    if math.isinf(value):
        return "-"
    return f"{value / AMAGAT:.3g}"


def _density_range(low: float, high: float) -> str:
    if low == 0 and math.isinf(high):
        return "no limit"
    if math.isinf(high):
        return f"above {_amagat(low)}"
    return f"{_amagat(low)}-{_amagat(high)}"


def terms_table() -> str:
    """Return each gas's response: its size, groups and ranges.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas in GASES:
        response = raman_response(gas, 295.0)
        branch = np.array(response.branch)
        q = branch == "Q"
        share = (
            np.sum(response.anisotropic[q]) / np.sum(response.isotropic[q])
            if q.any()
            else math.nan
        )
        low, high = response.density_range
        t_low, t_high = response.temperature_range
        rows.append(
            [
                gas,
                str(int(np.sum(branch == "S"))),
                str(int(np.sum(q))),
                str(len(set(response.damping_group.tolist()))),
                _density_range(low, high),
                f"{t_low:g}-{t_high:g}",
                "" if math.isnan(share) else f"{100 * share:.1f}%",
            ]
        )
    return _table(
        [
            "Gas",
            "S lines",
            "Q lines",
            "Damping groups",
            "Densities (amagat)",
            "Temperatures (K)",
            "Anisotropic Q share",
        ],
        rows,
    )


def _comparisons(comparisons, digits: str = ".3g") -> list[list[str]]:  # noqa: ANN001
    rows = []
    for c in comparisons:
        measured = f"{c.expected:{digits}}"
        if not math.isnan(c.sigma):
            measured += f" +- {c.sigma:{digits}}"
        rows.append(
            [
                c.material,
                c.label,
                c.reference,
                measured,
                f"{c.predicted:{digits}}",
                f"{100 * (c.ratio - 1):+.0f}%",
            ]
        )
    return rows


def gain_table() -> str:
    """Return the gases' gains from the response against the measured.

    Returns
    -------
    str
        A Markdown table.
    """
    return _table(
        [
            "Gas",
            "Line",
            "Measurement",
            "Measured (cm/GW)",
            "Response (cm/GW)",
            "Response / measured - 1",
        ],
        _comparisons(gas_gains()),
    )


def dephasing_table() -> str:
    """Return the decay of the rotational revivals against the measured.

    Returns
    -------
    str
        A Markdown table.
    """
    return _table(
        [
            "Gas",
            "Condition",
            "Measurement",
            "Measured (ps)",
            "Response (ps)",
            "Response / measured - 1",
        ],
        _comparisons(rotational_dephasing()),
    )


def _optmat_part(response, t: np.ndarray, density: float, part: str) -> np.ndarray:  # noqa: ANN001
    """Return the sum of a response's rotational or vibrational lines."""
    rotation = np.array(response.band) == "rotation"
    select = rotation if part == "rotation" else ~rotation
    omega = response.frequencies(density)[select]
    rate = response.decay_rates(density)[select]
    amplitude = response.amplitude[select]
    return np.array(
        [np.sum(amplitude * np.sin(omega * x) * np.exp(-rate * x)) for x in t]
    )


def solprop_table() -> str:
    """Return optmat's responses against solprop's, part by part, at 1 amagat.

    The early amplitude is the least-squares scale of optmat's response on solprop's
    over the first 100 fs, where no line has dephased; the late one the ratio of the
    two responses' largest magnitudes between 19.5 and 20.5 ps.

    Returns
    -------
    str
        A Markdown table.
    """
    data = json.loads(
        (ROOT / "tests" / "data" / "solprop_raman_response.json").read_text(
            encoding="utf-8"
        )
    )
    t = np.array(data["delays_s"])
    early = t <= 100e-15
    late = (t > 19.5e-12) & (t < 20.5e-12)
    rows = []
    for gas, parts in data["gases"].items():
        response = raman_response(gas, 295.0)
        for part, by_density in parts.items():
            reference = np.array(by_density["1.0"])
            ours = _optmat_part(response, t, AMAGAT, part)
            scale = float(
                np.dot(ours[early], reference[early])
                / np.dot(reference[early], reference[early])
            )
            ratio = float(np.max(np.abs(ours[late])) / np.max(np.abs(reference[late])))
            rows.append([gas, part, f"{scale:.3f}", f"{ratio:.3g}"])
    return _table(
        [
            "Gas",
            "Part",
            "Early amplitude: optmat / solprop",
            "At 20 ps: optmat / solprop",
        ],
        rows,
    )


def _rank(gas: str, densities: tuple[float, float], dt: float) -> list[int]:
    """Return the smallest number of separable terms reaching each tolerance."""
    response = raman_response(gas, 295.0)
    t = np.arange(RANK_SAMPLES) * dt
    grid = np.geomspace(*densities, RANK_DENSITIES) * AMAGAT
    spectra = []
    for density in grid:
        omega = response.frequencies(density)
        rate = response.decay_rates(density)
        h = np.zeros(RANK_SAMPLES)
        for start in range(0, len(response), 16):
            block = slice(start, start + 16)
            h += np.sum(
                response.amplitude[block]
                * np.sin(omega[block] * t[:, None])
                * np.exp(-rate[block] * t[:, None]),
                axis=1,
            )
        spectra.append(np.fft.rfft(h) * dt)
    matrix = np.array(spectra)
    u, s, vh = np.linalg.svd(matrix, full_matrices=False)
    scale = np.max(np.abs(matrix))
    ranks = []
    for tolerance in TOLERANCES:
        for r in range(1, len(s) + 1):
            approx = (u[:, :r] * s[:r]) @ vh[:r]
            if np.max(np.abs(matrix - approx)) <= tolerance * scale:
                ranks.append(r)
                break
    return ranks


def rank_table() -> str:
    """Return the rank of each response's family of kernels over a density range.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for gas, densities, dt in RANKS:
        response = raman_response(gas, 295.0)
        ranks = _rank(gas, densities, dt)
        rows.append(
            [
                gas,
                f"{densities[0]:g}-{densities[1]:g}",
                f"{RANK_SAMPLES * dt * 1e12:.0f} ps at {dt * 1e15:g} fs",
                str(len(response)),
                str(len(set(response.damping_group.tolist()))),
                *(str(r) for r in ranks),
            ]
        )
    return _table(
        [
            "Gas",
            "Densities (amagat)",
            "Window",
            "Lines",
            "Damping groups",
            *(f"Rank at {tol:.0e}" for tol in TOLERANCES),
        ],
        rows,
    )


def silica_models_table() -> str:
    """Return fused silica's models: f_R, the pairing, the gain peak.

    Returns
    -------
    str
        A Markdown table.
    """
    grid = np.linspace(1.0, 1300.0, 12991) * WAVENUMBER
    rows = []
    for name, model in silica_raman_responses().items():
        peak = int(np.argmax(np.imag(model.spectrum(grid))))
        gains = []
        for pump in (526e-9, 795.5e-9, 1064e-9):
            gains.append(f"{np.max(model.gain(pump, grid)) * 100 * 1e11:.2f}")
        rows.append(
            [
                name,
                f"{model.raman_fraction:.3f}",
                "yes" if model.paired else "no",
                "yes" if model.has_tensor else "no",
                f"{grid[peak] / WAVENUMBER:.0f}",
                *gains,
            ]
        )
    return _table(
        [
            "Model",
            "f_R",
            "Paired with the Kerr fit",
            "Isotropic and anisotropic parts",
            "Gain peak (cm-1)",
            "Peak gain at 526 nm (1e-11 cm/W)",
            "at 795.5 nm",
            "at 1064 nm",
        ],
        rows,
    )


def silica_validation_table() -> str:
    """Return fused silica's models against the measurements.

    Returns
    -------
    str
        A Markdown table.
    """
    rows = []
    for c in (*silica_gains(), *silica_raman_fractions(), *silica_low_frequency_gain()):
        measured = f"{c.expected:.3g}"
        if not math.isnan(c.sigma):
            measured += f" +- {c.sigma:.2g}"
        rows.append(
            [
                c.model,
                c.quantity,
                c.label,
                c.reference,
                measured,
                f"{c.predicted:.3g}",
                f"{100 * (c.ratio - 1):+.0f}%",
            ]
        )
    return _table(
        [
            "Model",
            "Quantity",
            "Condition",
            "Measurement",
            "Measured",
            "Model",
            "Model / measured - 1",
        ],
        rows,
    )


TABLES = {
    "response-terms": terms_table,
    "response-gain": gain_table,
    "response-dephasing": dephasing_table,
    "response-solprop": solprop_table,
    "response-rank": rank_table,
    "silica-models": silica_models_table,
    "silica-validation": silica_validation_table,
}


def main() -> None:
    """Print every table between its markers."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        for name, function in TABLES.items():
            print(f"<!-- {name}:start -->")
            print(function())
            print(f"<!-- {name}:end -->\n")


if __name__ == "__main__":
    main()

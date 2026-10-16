# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The numerics of ``tools/hitran/windows.py``, on synthetic lines (no HITRAN access).

The smoothed reference of OEP 3 is computed as an FFT convolution of the far lines plus
a direct sum of the near ones; these tests compare it with the direct sum over every
line, and check the absorption, the Voigt kernel and the helpers that split the bands
and the windows.
"""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parents[1]


def _tool():
    path = ROOT / "tools" / "hitran" / "windows.py"
    spec = importlib.util.spec_from_file_location("hitran_windows", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("hitran_windows", module)
    spec.loader.exec_module(module)
    return module


W = _tool()
SETTINGS = W.Settings(top_cm=4000.0)


@pytest.fixture(scope="module")
def lines():
    rng = np.random.default_rng(0)
    nu = np.sort(rng.uniform(1000, 3000, 3000))
    return nu, rng.uniform(0.1, 1, 3000) * 1e-20


@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_the_fft_reference_equals_the_direct_sum(lines, order):
    nu, strength = lines
    grid, far = W.far_derivatives(nu, strength, SETTINGS)
    points = np.array([800.0, 1500.0, 2000.3, 2999.5, 3100.0, 3500.0])
    index = np.rint(points / SETTINGS.fine_cm).astype(int)
    assert np.allclose(grid[index], points)
    total = far[:, index] + W.near_derivatives(nu, strength, SETTINGS, points)
    direct = W.direct_derivatives(nu, strength, SETTINGS.sigma_cm, points)
    # The far lines' deposition on the 0.005 cm^-1 grid errs by (h / 30 cm^-1)^2 of
    # their contribution, and the FFT's roundoff stays below the far values: 1e-8.
    np.testing.assert_allclose(total[order], direct[order], rtol=1e-8)


def test_the_smoothed_kernel_is_continuous_where_it_turns_asymptotic():
    sigma = 1.0 / W.HCM
    c = 1 / (math.sqrt(2) * sigma)
    for order in range(4):
        below = W._resonant_kernel(np.array([15 * (1 - 1e-13) / c]), sigma, order)[0]
        above = W._resonant_kernel(np.array([15 * (1 + 1e-13) / c]), sigma, order)[0]
        assert below == pytest.approx(above, rel=1e-9, abs=0.0)


def test_the_voigt_kernel_tends_to_the_gaussian_one():
    sigma = 1.0 / W.HCM
    delta = np.array([-20.0, -3.0, 0.7, 5.0, 20.0]) / W.HCM
    for order in range(4):
        np.testing.assert_allclose(
            W._voigt_kernel(delta, sigma, 1e-9 / W.HCM, order),
            W._resonant_kernel(delta, sigma, order),
            rtol=1e-6,
        )


def test_the_absorption_is_the_gaussian_line_sum(lines):
    nu, strength = lines
    size = math.ceil(SETTINGS.top_cm / SETTINGS.fine_cm) + 2
    absorption = W.smoothed_absorption(nu, strength, SETTINGS, size)
    points = np.array([1000.0, 1500.5, 2999.0, 3200.0])
    g = np.exp(-((points[:, None] - nu[None, :]) ** 2) / 2) / math.sqrt(2 * math.pi)
    density = 101325.0 / (W.K_B * 296.0) * 1e-6
    expected = density * np.sum(strength * g, axis=1) * 100
    got = absorption[np.rint(points / SETTINGS.fine_cm).astype(int)]
    # Linear deposition errs by (h / sigma)^2 / 8 = 3e-6 relative near the lines.
    np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-12)
    assert np.all(absorption >= 0)


def test_runs_and_pieces():
    nu = np.arange(10.0)
    mask = np.array([0, 1, 1, 0, 0, 1, 1, 1, 0, 1], bool)
    assert W.runs(mask, nu) == [(1.0, 2.0), (5.0, 7.0), (9.0, 9.0)]
    assert W._pieces((0.0, 10.0), [(2.0, 3.0), (6.0, 7.0)]) == [
        (0.0, 2.0),
        (3.0, 6.0),
        (7.0, 10.0),
    ]


@pytest.mark.parametrize("value", [1.23456e-7, 9.99999e3, 1.0, 0.000123401])
def test_a_stored_bound_is_rounded_up(value):
    stored = float(W._up(value))
    assert stored >= value
    assert stored <= value * (1 + 1e-3)


def test_bands_are_split_at_the_windows():
    lines = W.Lines(
        nu=np.array([100.0, 900.0, 1100.0, 1500.0, 2100.0]),
        strength=np.full(5, 1e-20),
        gamma_self=np.full(5, 0.1),
        band=np.array([0, 1, 1, 1, 1]),
        keys=[("1", "0 1", "0 1"), ("1", "1", "0")],
    )
    systems, inside = W.partition(lines, [(1000.0, 2000.0)])
    kinds = sorted((s.treatment, s.where, len(s.index)) for s in systems)
    assert kinds == [
        ("rotational", "gap 0", 1),
        ("rule", "gap 0", 1),
        ("rule", "gap 1", 1),
    ]
    assert [len(p.index) for p in inside[0]] == [2]

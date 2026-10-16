# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import itertools
import math

import numpy as np
import pytest
from scipy import constants

from optmat.constants import AU_POLARIZABILITY
from optmat.datasets import datasets
from optmat.exceptions import DataError
from optmat.gas.infrared import (
    InfraredBand,
    band_systems,
    fold_poles,
    gauss_rule,
    infrared_bands,
    infrared_polarizability_au,
    static_infrared_polarizability_au,
)
from optmat.gas.models.formulas import published_models
from optmat.gas.models.registry import default_model

HARTREE_WAVENUMBER = (
    constants.physical_constants["hartree-inverse meter relationship"][0] / 100
)


def test_oscillator_strength_is_the_integrated_cross_section_over_pi_r_e():
    band = InfraredBand("test", strength=1e-17, strength_over_nu2=1e-17 / 1000.0**2)
    # pi r_e = e^2 / (4 epsilon_0 m_e c^2), here in cm.
    pi_r_e = (
        constants.e**2
        / (4 * constants.epsilon_0 * constants.m_e * constants.c**2)
        * 100
    )
    assert band.oscillator_strength == pytest.approx(1e-17 / pi_r_e, rel=1e-12)
    assert band.wavenumber == pytest.approx(1000.0, rel=1e-12)


@pytest.mark.parametrize("gas", ["CH4", "N2O", "SF6"])
def test_static_value_is_the_sum_of_strengths_over_squared_wavenumbers(gas):
    bands = infrared_bands(gas)
    total = sum(band.strength_over_nu2 for band in bands)
    pi_r_e = (
        math.pi * constants.physical_constants["classical electron radius"][0] * 100
    )
    expected = total / pi_r_e * HARTREE_WAVENUMBER**2
    # The poles keep each band's sum S / nu^2 by construction: round-off only.
    assert static_infrared_polarizability_au(gas) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("gas", ["CH4", "N2O", "SF6"])
def test_optical_contribution_is_small_and_negative(gas):
    omega = np.array([0.114, 0.228])  # 400 and 200 nm, above every band (< 11000 cm^-1)
    alpha = infrared_polarizability_au(gas, omega)
    total_f = sum(band.oscillator_strength for band in infrared_bands(gas))
    # Far above the bands f / (E^2 - omega^2) tends to -f / omega^2; at 400 nm the
    # fundamentals (below 3100 cm^-1) are 1.5 % from that limit.
    np.testing.assert_allclose(alpha, -total_f / omega**2, rtol=0.03)


def test_gases_without_a_band_table_raise():
    with pytest.raises(DataError, match="no infrared bands"):
        infrared_bands("N2")


@pytest.mark.parametrize(("gas", "tolerance"), [("SF6", 0.05), ("CH4", 0.25)])
def test_shelton_ratios_at_two_wavelengths_confirm_the_band_strengths(gas, tolerance):
    # Shelton's N2/gas phase-matching ratio at 1319 and 1064 nm fixes
    # Delta alpha = alpha(2 omega) - alpha(omega) of the gas at both wavelengths. The
    # infrared share of it differs by a factor of 2.4 between them, so a single scale
    # on the electronic formula's dispersion fits both only if the band strengths are
    # right: solve for the scale k on HITRAN's strengths that makes them consistent.
    # k = 1.007 (SF6, infrared share 45 %) and 0.88 (CH4, share 6 %, so +-0.08).
    (ratio,) = [
        d
        for d in datasets(quantity="phase_match_density_ratio")
        if d.key == f"shelton_1990/density_N2_over_{gas}"
    ]
    # He et al.'s formula: electronic only (no infrared term); nitrogen as recommended.
    electronic = published_models()[f"he_2021/{gas}"]
    nitrogen = default_model("N2")
    normalised, infrared = [], []
    for wavelength, value in zip(
        ratio.vacuum_wavelength(), ratio.scaled_values(), strict=True
    ):
        points = np.array([wavelength, wavelength / 2])
        omega = 1e-2 / points / HARTREE_WAVENUMBER
        d_nitrogen = np.diff(nitrogen.polarizability(points, warn=False))[0]
        d_formula = np.diff(electronic.polarizability(points, warn=False))[0]
        d_infrared = np.diff(infrared_polarizability_au(gas, omega))[0]
        normalised.append(value * d_nitrogen / d_formula)
        infrared.append(d_infrared * AU_POLARIZABILITY / d_formula)
    scale = (normalised[0] - normalised[1]) / (infrared[0] - infrared[1])
    assert scale == pytest.approx(1.0, abs=tolerance)


@pytest.mark.parametrize("gas", ["CH4", "N2O", "CO2"])
def test_pure_rotational_bands_are_left_out(gas):
    # A band with the same vibrational state above and below is a rotational band: for
    # a polar molecule its static contribution is orientation, not polarization
    # (N2O: 1.45 a.u. at 296 K, TC-156).
    # Lines HITRAN leaves unassigned (no quanta: CH4 between 900 and 10900 cm^-1) are
    # kept; they are vibrational.
    for band in infrared_bands(gas):
        upper, separator, lower = band.label.partition(" <- ")
        if separator:
            assert upper.split(" ", 1)[1] != lower
    # The derived rules (OEP 3) keep them, as one rotational system of their own: a
    # term class a fit with band_nodes keeps apart from the infrared terms, and leaves
    # out of electric observations (tests/test_infrared_intervals.py).
    rotational = [s for s in band_systems(gas) if s.treatment == "rotational"]
    assert len(rotational) == 1


def _band(seed=1, count=200):
    # A synthetic band: lines spread over 2000-2400 cm^-1 with random strengths.
    rng = np.random.default_rng(seed)
    energies = np.sort(rng.uniform(2000, 2400, count)) / HARTREE_WAVENUMBER
    return rng.uniform(0.1, 1.0, count) * 1e-5, energies


def _moments(f, e, count):
    # mu_k of the measure f / E^2 at x = E^2, k = 0 .. count - 1.
    return np.array([np.sum(f / e**2 * (e**2) ** k) for k in range(count)])


@pytest.mark.parametrize("n", [1, 2, 4, 8])
def test_a_rule_reproduces_its_measures_moments(n):
    f, e = _band()
    strengths, nodes = gauss_rule(f, e, n)
    assert len(nodes) == n
    assert np.all(strengths > 0)
    assert e.min() < nodes.min() and nodes.max() < e.max()
    assert np.all(np.diff(nodes) > 0)
    # mu_0 .. mu_{2n-1}: the Lanczos process with full reorthogonalisation is exact to
    # a few ulps of the largest moment's scale.
    np.testing.assert_allclose(
        _moments(strengths, nodes, 2 * n), _moments(f, e, 2 * n), rtol=1e-12
    )


def test_a_rule_converges_to_the_line_sum_outside_the_band():
    f, e = _band()
    omega = np.array([1000.0, 1900.0, 2500.0, 4000.0]) / HARTREE_WAVENUMBER
    exact = np.sum(f / (e**2 - omega[:, None] ** 2), axis=1)
    errors = []
    for n in (1, 2, 4, 8):
        strengths, nodes = gauss_rule(f, e, n)
        rule = np.sum(strengths / (nodes**2 - omega[:, None] ** 2), axis=1)
        errors.append(np.max(np.abs(rule / exact - 1)))
    assert all(b < a for a, b in itertools.pairwise(errors))
    # Far from the band (1000 and 4000 cm^-1) eight nodes reach roundoff; 100 cm^-1
    # from its edges (1900, 2500 cm^-1) convergence is slower, as OEP 3 section 4 says.
    assert errors[-1] < 1e-6


def test_fewer_distinct_lines_than_nodes_give_the_exact_rule():
    # Three distinct energies, one of them twice: at most three nodes, exact.
    e = np.array([1.0, 1.5, 1.5, 2.0]) * 1e-2
    f = np.array([1.0, 0.5, 0.5, 2.0])
    strengths, nodes = gauss_rule(f, e, 8)
    assert len(nodes) == 3
    np.testing.assert_allclose(nodes, [1e-2, 1.5e-2, 2e-2], rtol=1e-12)
    np.testing.assert_allclose(strengths, [1.0, 1.0, 2.0], rtol=1e-12)
    single, energy = gauss_rule(np.array([3.0]), np.array([0.01]), 4)
    assert (float(single[0]), float(energy[0])) == pytest.approx((3.0, 0.01))


def test_a_band_split_at_a_window_sums_to_the_whole():
    f, e = _band()
    cut = e < 2200 / HARTREE_WAVENUMBER
    parts = [gauss_rule(f[m], e[m], 4) for m in (cut, ~cut)]
    for k in range(2):
        whole = _moments(f, e, 2)[k]
        split = sum(_moments(s, x, 2)[k] for s, x in parts)
        assert split == pytest.approx(whole, rel=1e-12)


@pytest.mark.parametrize(
    ("strengths", "energies", "n"),
    [([1.0], [1.0], 0), ([], [], 1), ([-1.0], [1.0], 1), ([1.0], [0.0], 1)],
)
def test_a_rule_refuses_bad_lines(strengths, energies, n):
    with pytest.raises(ValueError, match="positive strength"):
        gauss_rule(np.array(strengths, dtype=float), np.array(energies, dtype=float), n)


def test_a_fold_keeps_the_bands_two_moments():
    f, e = _band()
    m0, m_minus2 = float(f.sum()), float(np.sum(f / e**2))
    below, above = 1500 / HARTREE_WAVENUMBER, 3000 / HARTREE_WAVENUMBER
    low, high = fold_poles(m0, m_minus2, below, above)
    assert low > 0 and high > 0
    assert low + high == pytest.approx(m0, rel=1e-12)
    assert low / below**2 + high / above**2 == pytest.approx(m_minus2, rel=1e-12)


@pytest.mark.parametrize(("below", "above"), [(2.5, 3.0), (1.0, 1.5), (2.0, 1.0)])
def test_a_fold_needs_the_mean_between_its_poles(below, above):
    # The band's mean energy is 2 (M0 = 4, M_-2 = 1).
    with pytest.raises(ValueError, match="not between"):
        fold_poles(4.0, 1.0, below, above)

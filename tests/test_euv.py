# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

from optmat.datasets import datasets
from optmat.exceptions import DataError
from optmat.gas.euv import (
    EV_PER_HARTREE,
    MB_PER_INVERSE_EV,
    OscillatorDistribution,
    constrained,
    euv_model,
)
from optmat.gas.gases import Gas
from optmat.gas.models.fitted import fitted_models
from optmat.gas.models.oscillator import angular_frequency_au


def _single_line(f=0.5, e=1.0):
    # A line plus a flat continuum from 2 to 4 hartree (and its tail).
    return OscillatorDistribution(
        gas=Gas.He,
        name="toy",
        line_energies=np.array([e]),
        line_strengths=np.array([f]),
        continuum_energies=np.array([2.0, 3.0, 4.0]),
        continuum_density=np.array([0.3, 0.3, 0.3]),
    )


def test_cross_section_conversion_is_the_known_constant():
    # pi e^2 hbar / (2 eps0 m_e c) = 109.761 Mb eV.
    assert pytest.approx(109.761, abs=1e-3) == MB_PER_INVERSE_EV


def test_below_threshold_the_real_part_is_the_static_sum_series():
    toy = _single_line()
    omega = np.array([0.0, 0.1])
    real, imaginary = toy.complex_polarizability_au(omega)
    # At omega = 0 the real part is S(-2) exactly (the same quadrature).
    assert real[0] == pytest.approx(toy.sum_rule(-2), rel=1e-12, abs=0.0)
    assert np.all(imaginary == 0)


def test_above_threshold_the_imaginary_part_is_the_absorption():
    toy = _single_line()
    omega = np.array([2.5, 3.5])
    _, imaginary = toy.complex_polarizability_au(omega)
    np.testing.assert_allclose(imaginary, np.pi / (2 * omega) * 0.3, rtol=1e-12)


def test_principal_value_is_continuous_across_a_tabulated_energy():
    toy = _single_line()
    omega = np.array([2.999999, 3.0, 3.000001])
    real, _ = toy.complex_polarizability_au(omega)
    # The subtraction makes the integrand regular: a smooth function of omega.
    assert abs(real[1] - (real[0] + real[2]) / 2) < 1e-5


def test_constraint_meets_its_target_and_scales_only_below_the_split():
    toy = _single_line()
    scaled, factor = constrained(toy, 3.5, (-2, toy.sum_rule(-2) * 1.01))
    assert scaled.sum_rule(-2) == pytest.approx(
        toy.sum_rule(-2) * 1.01, rel=1e-10, abs=0.0
    )
    assert scaled.continuum_density[-1] == toy.continuum_density[-1]
    assert factor > 1


def _analytic_static_target(factor):
    # The synthetic continuum is 0.3*f on [2,3], log-log interpolated to 0.3 at 4,
    # then proportional to E^-3.5. Integrate g(E)/E^2 over each interval exactly.
    exponent = -np.log(factor) / np.log(4 / 3)
    join = 0.3 * factor / 3 * ((4 / 3) ** (exponent - 1) - 1) / (exponent - 1)
    return 0.5 + 0.3 * factor / 6 + join + 0.3 / 18


def test_strong_continuum_constraint_matches_an_independent_integral():
    toy = _single_line()
    target = _analytic_static_target(0.2)
    scaled, factor = constrained(toy, 3.5, (-2, target))
    # The 20,000-point log grid has spacing ~5e-4: its trapezoidal error is O(h^2).
    assert factor == pytest.approx(0.2, rel=1e-6, abs=0.0)
    np.testing.assert_allclose(
        scaled.continuum_density, [0.06, 0.06, 0.3], rtol=1e-6, atol=0
    )
    # Independently integrate the returned curve, without sum_rule's quadrature.
    assert _analytic_static_target(factor) == pytest.approx(target, rel=1e-6, abs=0.0)
    np.testing.assert_array_equal(scaled.line_strengths, toy.line_strengths)


def test_constraint_rejects_a_nonphysical_trial_scale_before_taking_its_log():
    # This valid synthetic target needs more robust root bracketing than the local
    # secant solver provides. It must report failure, not return a NaN distribution.
    with pytest.raises(DataError, match="invalid continuum scale"):
        constrained(_single_line(), 3.5, (-2, _analytic_static_target(0.05)))


def test_constraint_cannot_change_a_sum_without_an_adjustable_continuum():
    toy = _single_line()
    with pytest.raises(DataError, match="does not respond"):
        constrained(toy, 1.5, (-2, _analytic_static_target(0.5)))


def test_helium_distribution_meets_the_sum_rules():
    model = euv_model("He")
    # Scaled within Samson's 1-2 %; TRK within 2 % (the high-energy data).
    assert 1.0 < model.factor < 1.02
    assert model.distribution.sum_rule(0) == pytest.approx(2.0, rel=0.02, abs=0.0)


def test_helium_below_threshold_agrees_with_the_ab_initio_fit():
    model = euv_model("He")
    fit = fitted_models()["fit/He"]
    wavelength = np.array([1e-6, 200e-9, 100e-9, 70e-9])
    omega = angular_frequency_au(wavelength)
    dosd = np.real(model.polarizability_au(omega))
    # The Kramers-Kronig transform of measured photoabsorption reproduces the ab initio
    # refractive polarizability to 0.02 % in the infrared, 0.11 % at 70 nm.
    np.testing.assert_allclose(dosd, fit.polarizability_au(omega), rtol=1.5e-3)


def test_helium_absorption_agrees_with_liu():
    (liu,) = [d for d in datasets("He") if d.key == "liu_1997/alpha_imaginary"]
    omega = np.array([float(x) for x in liu.abscissae]) / EV_PER_HARTREE
    imaginary = np.imag(euv_model("He").polarizability_au(omega))
    # Samson's cross sections (1-2 %) against Liu's ab initio (a few 1e-3).
    np.testing.assert_allclose(imaginary, liu.scaled_values(), rtol=0.02)


@pytest.mark.parametrize(
    ("gas", "electrons", "scale"),
    # Chan et al.'s Ar-Xe continua are absolute to 5 %: the constraint moves them
    # less. Their neon continuum is TRK-normalised, not absolute: it moves 5.4 %.
    [("Ne", 10, 0.06), ("Ar", 18, 0.05), ("Kr", 36, 0.05), ("Xe", 54, 0.05)],
)
def test_heavier_noble_gases_meet_their_sums_with_small_scales(gas, electrons, scale):
    model = euv_model(gas)
    assert abs(model.factor - 1) < scale
    # TRK within 7 %: the high-energy (Henke) and autoionising-region approximations.
    assert model.distribution.sum_rule(0) == pytest.approx(electrons, rel=0.07, abs=0.0)


@pytest.mark.parametrize("gas", ["Ne", "Ar", "Kr", "Xe"])
def test_heavier_noble_gases_below_threshold_agree_with_the_fits(gas):
    model = euv_model(gas)
    fit = fitted_models()[f"fit/{gas}"]
    omega = angular_frequency_au(np.array([2e-6, 500e-9]))
    # The transform of the photoabsorption reproduces the fitted refraction in the
    # visible and infrared to 0.15 % (krypton's worst, 0.13 %).
    np.testing.assert_allclose(
        np.real(model.polarizability_au(omega)),
        fit.polarizability_au(omega),
        rtol=1.5e-3,
    )


def test_argon_dispersion_moment_matches_ab_initio():
    # S(-4) of the constrained distribution against Lesiuk & Jeziorski's 27.976(15).
    assert euv_model("Ar").distribution.sum_rule(-4) == pytest.approx(27.976, abs=0.05)


def test_neon_dispersion_moment_from_photoabsorption_is_below_ab_initio():
    # S(-4) of the constrained distribution, 2.826, against Lesiuk et al.'s 2.850(7):
    # the photoabsorption gives less dispersion than theory (TC-126), where
    # refractometry gives more (TC-92). A regression bound, not an agreement test.
    assert euv_model("Ne").distribution.sum_rule(-4) == pytest.approx(2.826, abs=0.005)

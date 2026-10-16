# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The measured-absorption term: Kramers-Kronig transform, moments, derivatives."""

import numpy as np
import pytest

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_POLARIZABILITY, EPS0, C
from optmat.exceptions import DataError
from optmat.gas.absorption import Piece, assemble, measured_absorption


@pytest.fixture(scope="module")
def ozone():
    return measured_absorption("O3")


def _omega(wavelength_nm):
    return 2 * np.pi * C / (np.asarray(wavelength_nm) * 1e-9) / AU_ANGULAR_FREQUENCY


def test_the_transform_agrees_with_a_direct_principal_value_sum(ozone):
    # Off-grid points half a step from the samples: the plain sum is then a
    # principal value (the midpoint rule), independent of the FFT convolution.
    grid = ozone.wavenumber_cm * 100
    for wavelength in (300.0, 600.0, 800.0):
        index = int(np.argmin(np.abs(ozone.wavenumber_cm - 1e7 / wavelength)))
        nu = grid[index] + 0.5 * ozone.step_cm * 100
        direct = (
            EPS0
            / np.pi**2
            * np.sum(ozone.sigma_m2 / (grid**2 - nu**2))
            * ozone.step_cm
            * 100
            / AU_POLARIZABILITY
        )
        omega = nu * 2 * np.pi * C / AU_ANGULAR_FREQUENCY
        real = ozone.complex_polarizability_au(np.array([omega]))[0][0]
        assert real == pytest.approx(direct, rel=2e-4)


def test_the_static_limit_is_the_inverse_square_moment(ozone):
    static = ozone.complex_polarizability_au(np.array([0.0]))[0][0]
    assert static == pytest.approx(ozone.sum_rule(-2), rel=1e-3)


def test_the_imaginary_part_follows_the_cross_section(ozone):
    omega = _omega([255.5])
    imaginary = ozone.complex_polarizability_au(omega)[1][0]
    sigma = np.interp(1e7 / 255.5, ozone.wavenumber_cm, ozone.sigma_m2)
    expected = EPS0 * C * sigma / (omega[0] * AU_ANGULAR_FREQUENCY) / AU_POLARIZABILITY
    assert imaginary == pytest.approx(expected, rel=1e-12)


def test_above_the_measured_range_the_term_is_real_and_regular(ozone):
    real, imaginary = ozone.complex_polarizability_au(_omega([50.0, 80.0]))
    assert np.all(imaginary == 0)
    assert np.all(np.isfinite(real))


def test_the_hartley_band_strength_matches_kumar_and_thakkar(ozone):
    # Their 4.352-5.391 eV interval holds 0.0787 (distribution C, Table S of the
    # supplement, as the ozone project extracted it): within the cross sections' few %.
    assert ozone.strength(4.352, 5.391) == pytest.approx(0.0787, rel=0.04)


def test_the_first_derivative_matches_a_difference_of_the_real_part(ozone):
    omega = _omega(np.array([450.0, 700.0]))
    step = 1e-4
    real_plus = ozone.complex_polarizability_au(omega + step)[0]
    real_minus = ozone.complex_polarizability_au(omega - step)[0]
    numeric = (real_plus - real_minus) / (2 * step)
    analytic = ozone.real_derivative_au(omega, 1)
    np.testing.assert_allclose(analytic, numeric, rtol=2e-3)


def test_smoothing_leaves_the_static_value_nearly_unchanged():
    wide = measured_absorption("O3", fwhm_cm=200.0)
    narrow = measured_absorption("O3", fwhm_cm=10.0)
    assert wide.sum_rule(-2) == pytest.approx(narrow.sum_rule(-2), rel=1e-3)


def test_overlapping_pieces_are_refused():
    source = "mainz_uvvis/O3_Bogumil2003_293K_230-1070nm.txt.gz"
    with pytest.raises(DataError, match="overlap"):
        assemble(
            "O3",
            (Piece(source, 300.0, 400.0), Piece(source, 350.0, 500.0)),
            temperature_k=293.0,
        )


def test_the_temperature_variant_replaces_the_serdyuchenko_window():
    cold = measured_absorption("O3", variant="T243")
    assert cold.temperature_k == 243.0
    assert any("243K" in p.source for p in cold.pieces)

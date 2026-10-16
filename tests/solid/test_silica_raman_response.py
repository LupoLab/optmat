# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import numpy as np
import pytest

from optmat.constants import EPS0, C
from optmat.exceptions import DataError
from optmat.solid.kerr import kerr_models
from optmat.solid.raman_response import silica_raman_response, silica_raman_responses

WAVENUMBER = 2 * math.pi * C * 100  # rad/s per cm^-1


@pytest.mark.parametrize("name", ["hollenbeck_2002", "lin_2006", "blow_1989"])
def test_each_shape_has_unit_integral(name):
    response = silica_raman_response(name)
    # 0.02 fs steps to 8 ps: the trapezoid's error on the fastest period (27 fs) is
    # (dt omega)^2/12 = 2e-7 of the largest term; Lin & Agrawal's h_b jumps to 2/tau_b
    # at zero, which the trapezoid on a response that is 0 at t = 0 misses by
    # f_b dt h_b(0+)/2 = 0.21 x 0.02/96 = 4.4e-5.
    t = np.arange(0.0, 8e-12, 0.02e-15)
    assert np.trapezoid(response.response(t), t) == pytest.approx(1.0, abs=5e-5)
    assert np.real(response.spectrum(0.0)) == pytest.approx(1.0, abs=1e-12)
    assert response.response(np.array([-1e-15, 0.0]))[0] == 0.0


@pytest.mark.parametrize("name", ["hollenbeck_2002", "lin_2006", "blow_1989"])
def test_the_spectrum_is_the_transform_of_the_response(name):
    response = silica_raman_response(name)
    dt = 0.05e-15
    t = np.arange(0.0, 6e-12, dt)
    h = response.response(t)
    for shift in (100.0, 440.0, 800.0):
        omega = shift * WAVENUMBER
        numeric = np.trapezoid(h * np.exp(1j * omega * t), t)
        # The trapezoid on a 0.05 fs grid: (dt omega)^2/12 < 1e-7, plus Lin &
        # Agrawal's jump at zero (dt h_b(0+)/2 = 1e-4 of the scale).
        assert numeric == pytest.approx(response.spectrum(omega), abs=3e-4)


def test_the_gain_of_hollenbecks_shape_peaks_at_440_wavenumbers():
    response = silica_raman_response()
    shift = np.linspace(1.0, 1300.0, 13000)
    gain = np.imag(response.spectrum(shift * WAVENUMBER))
    # Stolen et al.'s measured peak is at about 440 cm-1; Hollenbeck & Cantrell's
    # fit puts it at 439-448 cm-1 (TC-382).
    assert 430 < shift[np.argmax(gain)] < 450
    assert np.all(gain > -1e-12)


def test_lin_agrawals_parts():
    response = silica_raman_response("lin_2006")
    t = np.linspace(0.0, 2e-12, 1001)
    parallel = response.response(t)
    a = response.response(t, "a")
    b = response.response(t, "b")
    assert a + b == pytest.approx(parallel, rel=1e-12, abs=1e-3)
    # f_a = 0.75 isotropic, f_b + f_c = 0.25 anisotropic (their Eqs. (1)-(4)).
    assert np.real(response.spectrum(0.0, "a")) == pytest.approx(0.75)
    assert np.real(response.spectrum(0.0, "b")) == pytest.approx(0.25)


def test_a_copolarised_model_has_no_parts():
    with pytest.raises(DataError, match="copolarised only"):
        silica_raman_response().response(np.zeros(3), "a")
    with pytest.raises(DataError, match="no silica Raman model"):
        silica_raman_response("nope")


def test_the_default_is_paired_with_the_kerr_fit():
    response = silica_raman_response()
    kerr = kerr_models()["kerr/SiO2"]
    assert response.paired
    assert response.raman_fraction == kerr.raman_fraction
    assert not silica_raman_response("lin_2006").paired


def test_the_field_and_the_envelope_forms_add_up_to_n2():
    # A kernel slower than the cycle follows half of E^2: n2 = 3 chi_e / (4 n0^2 eps0
    # c) + K / (2 n0^2 eps0 c), the long-pulse index of the module docstring.
    response = silica_raman_response()
    wavelength = 1030e-9
    n0 = response._index_and_n2(wavelength)[0]
    n2 = kerr_models()["kerr/SiO2"].n2(wavelength)
    electronic = (
        3 * response.electronic_susceptibility(wavelength) / (4 * n0**2 * EPS0 * C)
    )
    raman = response.raman_coefficient(wavelength) / (2 * n0**2 * EPS0 * C)
    assert electronic + raman == pytest.approx(float(n2), rel=1e-14)
    assert raman / (electronic + raman) == pytest.approx(response.raman_fraction)


def test_blow_and_woods_share_of_the_envelope():
    # alpha = 0.7 on E^2, and 2/3 of 1 - alpha in the envelope (their Eqs. (2), (8)).
    assert silica_raman_response("blow_1989").raman_fraction == pytest.approx(0.2 / 0.9)


def test_lin_and_agrawal_reproduce_their_peak_gain():
    # Their f_R = 0.245 was set with n2 = 2.6e-20 m^2/W to 1.2e-11 cm/W at 795.5 nm and
    # 1.81e-11 cm/W at 526 nm (p. 3087), with gamma = n2 omega_0 / c at the pump: this
    # gain uses the Stokes frequency, lower by the shift (5.5 % and 3.5 % at 434 cm-1).
    response = silica_raman_response("lin_2006")
    shift = np.linspace(300.0, 600.0, 3001) * WAVENUMBER
    for pump, printed in ((795.5e-9, 1.2e-13), (526e-9, 1.81e-13)):
        gain = response.gain(pump, shift, n2=2.6e-20)
        peak = int(np.argmax(gain))
        at_pump = (
            gain[peak]
            * (2 * math.pi * C / pump)
            / (2 * math.pi * C / pump - shift[peak])
        )
        # Printed to two and three figures: within their rounding.
        assert at_pump == pytest.approx(printed, rel=0.01)


def test_every_model_builds():
    assert set(silica_raman_responses()) == {"hollenbeck_2002", "lin_2006", "blow_1989"}

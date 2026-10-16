# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import json
import math

import numpy as np
import pytest

from optmat.exceptions import DataError
from optmat.solid import kerr

STORED = kerr.KERR_DIRECTORY / "SiO2.json"


def test_the_raman_response_is_hollenbecks_normalised():
    # Normalised to unit integral analytically (TC-691); the trapezoid on 0.1 fs, to
    # 4 ps where it has decayed, is 4e-6 low (its error on the 27 fs period of the
    # fastest mode); the first maximum at 16 fs is the paper's (checked at
    # extraction, TC-393).
    t = np.arange(0, 4e-12, 0.1e-15)
    h = kerr.raman_response(t)
    assert np.trapezoid(h, t) == pytest.approx(1.0, abs=1e-5)
    assert t[np.argmax(h)] == pytest.approx(16.3e-15, abs=0.5e-15)
    assert kerr.raman_response(np.array([-1e-15]))[0] == 0.0


def test_eta_matches_the_direct_time_domain_integral():
    # The spectral formula against the defining double integral, for a 60 fs pulse on
    # a 0.2 fs grid (quadrature: 1e-4).
    fwhm = 60e-15
    sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))
    dt = 0.2e-15
    t = np.arange(-6 * sigma, 6 * sigma, dt)
    s = np.arange(0, 2e-12, dt)
    h = kerr.raman_response(s)
    intensity = np.exp(-(t**2) / (2 * sigma**2))
    delayed = np.convolve(intensity, h)[: len(t)] * dt
    direct = np.sum(intensity * delayed) / np.sum(intensity**2)
    assert kerr.raman_weight(fwhm) == pytest.approx(direct, abs=1e-4)


def test_eta_goes_from_zero_for_short_pulses_to_one_for_long():
    assert kerr.raman_weight(2e-15) < 0.05
    assert kerr.raman_weight(100e-15) == pytest.approx(0.957, abs=0.002)
    assert kerr.raman_weight(20e-12) == pytest.approx(1.0, abs=1e-12)


def test_the_stored_fit_is_current_and_reproduces():
    config = kerr.load_config("SiO2")
    record = json.loads(STORED.read_text())
    assert record["config_sha256"] == config.digest
    result = kerr.fit(config)
    fitted = result.model.parameters
    stored = record["parameters"]
    # The fit is deterministic; 1e-6 relative allows for the optimiser's tolerance.
    assert fitted[0] == pytest.approx(stored["N_m2_per_W"], rel=1e-6)
    assert fitted[3] == pytest.approx(stored["raman_fraction"], rel=1e-6)


def test_the_silica_model_is_physical():
    model = kerr.kerr_models()["kerr/SiO2"]
    wavelength = np.linspace(360e-9, 3400e-9, 50)
    n2 = model.n2(wavelength)
    assert np.all(n2 > 0)
    # Normal dispersion of n2 below the two-photon edge: it falls with wavelength.
    assert np.all(np.diff(n2) < 0)
    assert 0 < model.raman_fraction < 0.3
    electronic = model.n2_electronic(800e-9)
    assert electronic == pytest.approx((1 - model.raman_fraction) * model.n2(800e-9))
    assert model.n2_effective(800e-9, 1e-9) == pytest.approx(model.n2(800e-9))
    assert np.all(model.n2_sigma(wavelength) > 0)


def test_the_fitted_gain_is_the_silica_responses():
    # The Kerr fit predicts Stolen & Ippen's gain as 2 omega_s n2 f_R Im H_R / c; the
    # stored fit's f_R and n2 in optmat.solid.raman_response give the same (TC-696).
    from optmat.solid.raman_response import silica_raman_response  # noqa: PLC0415

    (options,) = kerr.load_config("SiO2").raman_gains
    gain = kerr.gain_points(options)
    model = kerr.kerr_models()["kerr/SiO2"]
    predicted = gain.factor() * model.raman_fraction * float(model.n2(gain.wavelength))
    response = silica_raman_response()
    assert float(response.gain(gain.wavelength, gain.shift)) == pytest.approx(
        predicted, rel=1e-12
    )
    # One standard uncertainty (10 %) below the measurement, where the n2 data pull.
    assert predicted / gain.value == pytest.approx(0.9, abs=0.03)


def test_a_config_needs_a_prior_or_a_gain(tmp_path):
    text = (kerr.KERR_DIRECTORY / "SiO2.toml").read_text()
    head = text.split("[[raman_gains]]")[0]
    path = tmp_path / "bare.toml"
    path.write_text(head)
    with pytest.raises(DataError, match="needs a prior"):
        kerr.load_config(path)

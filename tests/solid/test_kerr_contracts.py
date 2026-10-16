# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Controlled Kerr observations separate units, pulse weighting and calibration."""

import dataclasses
import json

import numpy as np
import pytest

from optmat.constants import E_CHARGE, HBAR, C
from optmat.datasets import Uncertainty, load_source
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.solid import kerr


@pytest.fixture
def sample(tmp_path, monkeypatch):
    wavelength = np.array([400, 500, 600, 800, 1000, 1200], dtype=float)
    energy = 2 * np.pi * HBAR * C / (wavelength * 1e-9 * E_CHARGE)
    points = [
        [str(i), str(w), str(3e-20 * (1 + 0.03 * e**2))]
        for i, (w, e) in enumerate(zip(wavelength, energy, strict=True))
    ]
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[datasets]]\nid="n2"\n'
        'material="SiO2"\nquantity="nonlinear_index"\nunit="m2/W"\nabscissa_unit="nm"\n'
        f"points={json.dumps(points)}\nuncertainty={{relative=0.01}}\n"
        'parameters={part="total",coefficient="self"}\n'
    )
    sample = load_source(path).datasets[0]
    monkeypatch.setattr(kerr, "datasets", lambda: (sample,))
    return sample


def test_pulse_metadata_is_required_only_for_pulse_weighting(sample):
    with pytest.raises(DataError, match="no pulse length"):
        kerr.points(kerr.KerrOptions(sample.key))
    no_raman = kerr.points(kerr.KerrOptions(sample.key, raman="none"))
    full_raman = kerr.points(kerr.KerrOptions(sample.key, raman="full"))
    np.testing.assert_array_equal(no_raman.eta, np.zeros(6))
    np.testing.assert_array_equal(full_raman.eta, np.ones(6))
    # At 20 ps the response's documented long-pulse limit is unity.
    long = kerr.points(kerr.KerrOptions(sample.key, pulse_fs=20000))
    np.testing.assert_allclose(long.eta, np.ones(6), atol=1e-12, rtol=0)
    np.testing.assert_array_equal(no_raman.value, full_raman.value)


def test_kerr_observations_require_self_coefficients_and_uncertainty(
    sample, monkeypatch
):
    changed = dataclasses.replace(
        sample, parameters={"part": "total", "coefficient": "probe"}
    )
    monkeypatch.setattr(kerr, "datasets", lambda: (changed,))
    with pytest.raises(DataError, match="only self coefficients"):
        kerr.points(kerr.KerrOptions(sample.key, raman="full"))
    changed = dataclasses.replace(sample, uncertainty=Uncertainty())
    with pytest.raises(DataError, match="points without an uncertainty"):
        kerr.points(kerr.KerrOptions(sample.key, raman="full"))
    supplied = kerr.points(kerr.KerrOptions(sample.key, raman="full", sigma=0.02))
    np.testing.assert_allclose(
        supplied.sigma, supplied.value * 0.02, rtol=2e-15, atol=0
    )
    changed = dataclasses.replace(
        sample,
        quantity="nonlinear_index_field",
        unit="m2/V2",
        parameters=sample.parameters | {"definition": "synthetic"},
    )
    with pytest.raises(DataError, match="is not converted"):
        kerr.points(kerr.KerrOptions(sample.key, raman="full"))


@pytest.mark.parametrize(
    ("entry", "message"),
    [
        ({"key": "test", "typo": 1}, "unknown keys"),
        ({"key": "test", "raman": "automatic"}, "unknown raman treatment"),
    ],
)
def test_kerr_options_do_not_silently_ignore_unknown_settings(entry, message):
    with pytest.raises(DataError, match=message):
        kerr._kerr_options(entry, "synthetic")


@pytest.mark.parametrize("quartic", [False, True])
def test_synthetic_kerr_fit_with_fraction_prior_round_trips(
    sample, tmp_path, monkeypatch, quartic
):
    path = tmp_path / "synthetic.toml"
    text = (
        'material="SiO2"\nrange_nm=[399,1201]\nraman_fraction=[0.2,0.01]\n'
        f"quartic={str(quartic).lower()}\n"
        'fitted=[{key="synthetic/n2",raman="full",group="synthetic",scale_sigma=0.01}]\n'
    )
    path.write_text(text)
    result = kerr.fit(kerr.load_config(path))
    record = json.loads(kerr.dumps(result))
    # Exact noiseless quadratic-in-energy data: solver rounding, including log scale.
    assert np.array(result.model.parameters)[[0, 1, 3]] == pytest.approx(
        [3e-20, 0.03, 0.2], rel=1e-6, abs=0
    )
    # |b| < 1e-10 changes n2 by less than 1e-8 over this energy range.
    assert result.model.parameters[2] == pytest.approx(0, abs=1e-10)
    np.testing.assert_allclose(
        result.model.n2(sample.vacuum_wavelength()),
        sample.scaled_values(),
        rtol=1e-7,
        atol=0,
    )
    assert record["config_sha256"] == result.config.digest
    (tmp_path / "synthetic.json").write_text(kerr.dumps(result))
    monkeypatch.setattr(kerr, "KERR_DIRECTORY", tmp_path)
    kerr.kerr_models.cache_clear()
    try:
        rebuilt = kerr.kerr_models()["kerr/synthetic"]
        np.testing.assert_array_equal(rebuilt.n2(800e-9), result.model.n2(800e-9))
        with pytest.warns(ModelRangeWarning, match="outside"):
            rebuilt.n2(200e-9)
    finally:
        kerr.kerr_models.cache_clear()
    for invalid, message in [
        ("typo=1\n" + text, "unknown keys"),
        (text.replace("synthetic/n2", "missing"), "unknown dataset"),
        (text + 'raman_gains=[{key="synthetic/n2",typo=1}]\n', "unknown keys"),
    ]:
        path.write_text(invalid)
        with pytest.raises(DataError, match=message):
            kerr.load_config(path)


def test_gain_observations_require_gain_units_pump_shift_and_errors(
    sample, monkeypatch
):
    with pytest.raises(DataError, match="one raman_gain"):
        kerr.gain_points(kerr.GainOptions(sample.key))
    gain = dataclasses.replace(
        sample,
        quantity="raman_gain",
        unit="cm/W",
        rows=("0",),
        abscissae=("500",),
        values=("1e-11",),
        sigmas=("",),
        parameters={"transition": "synthetic"},
    )
    monkeypatch.setattr(kerr, "datasets", lambda: (gain,))
    with pytest.raises(DataError, match=r"parameters\.pump_nm"):
        kerr.gain_points(kerr.GainOptions(gain.key))
    gain = dataclasses.replace(
        gain, parameters={"transition": "synthetic", "pump_nm": 500}
    )
    with pytest.raises(DataError, match=r"parameters\.shift_cm1"):
        kerr.gain_points(kerr.GainOptions(gain.key))
    gain = dataclasses.replace(
        gain,
        parameters={"transition": "synthetic", "pump_nm": 500, "shift_cm1": 440},
        uncertainty=Uncertainty(),
    )
    with pytest.raises(DataError, match="no uncertainty"):
        kerr.gain_points(kerr.GainOptions(gain.key))
    point = kerr.gain_points(kerr.GainOptions(gain.key, sigma=0.1))
    assert point.value == pytest.approx(1e-13, rel=2e-15, abs=0)
    assert point.sigma == pytest.approx(1e-14, rel=2e-15, abs=0)

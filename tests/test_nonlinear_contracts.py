# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Synthetic third-order data exercise conversions, fit priors and diagnostics."""

import dataclasses
import hashlib
import json

import numpy as np
import pytest

from optmat.conditions import State
from optmat.datasets import Uncertainty, load_source
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import nonlinear
from optmat.gas.gases import Gas


@pytest.fixture
def sample(tmp_path):
    # Invented Kerr dispersion, gamma = 10 * (1 + 2 W2*1e-10 + 3 W2^2*1e-20).
    wavelengths = np.array([400, 500, 600, 800, 1000, 1500], dtype=float)
    x = 4 * (1e7 / wavelengths) ** 2 * 1e-10
    values = 10 * (1 + 2 * x + 3 * x * x)
    points = [
        [str(i), str(w), str(v)]
        for i, (w, v) in enumerate(zip(wavelengths, values, strict=True))
    ]
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[datasets]]\nid="gamma"\n'
        'material="He"\nquantity="hyperpolarizability"\nunit="au"\nabscissa_unit="nm"\n'
        f"points={json.dumps(points)}\nuncertainty={{relative=0.001}}\n"
        'parameters={process="kerr",convention="T",part="electronic"}\n',
        encoding="utf-8",
    )
    return load_source(path).datasets[0]


@pytest.mark.parametrize("process", ["cars", "pump_probe"])
def test_two_colour_processes_require_both_frequencies(process):
    with pytest.raises(DataError, match="needs a second frequency"):
        nonlinear.frequency_sums(process, np.array([1.0]))


def test_unknown_third_order_process_is_not_treated_as_kerr():
    with pytest.raises(DataError, match="unknown third-order process"):
        nonlinear.frequency_sums("missing", np.array([1.0]))


def test_nonlinear_index_unit_validation_precedes_conversion(sample):
    changes = {
        "quantity": "nonlinear_index",
        "state": State(101325, 300),
        "parameters": {"part": "electronic", "coefficient": "self"},
    }
    for unit in ("cm2/W", "m2/W"):
        assert dataclasses.replace(sample, unit=unit, **changes).unit == unit
    with pytest.raises(DataError, match=r"unit 'esu'.*nonlinear_index"):
        dataclasses.replace(sample, unit="esu", **changes)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (
            {
                "unit": "SI",
                "parameters": {
                    "process": "kerr",
                    "convention": "ward_1969",
                    "part": "electronic",
                },
            },
            "no conversion",
        ),
        ({"quantity": "nonlinear_susceptibility"}, "only per-amagat CARS"),
        ({"quantity": "polarizability", "parameters": {}}, "not converted to gamma"),
        ({"abscissae": ("not stated",) * 6}, "wavelength is not stated"),
    ],
)
def test_unsupported_conversions_do_not_assume_units(sample, changes, message):
    with pytest.raises(DataError, match=message):
        nonlinear.gamma_observation(
            dataclasses.replace(sample, **changes),
            nonlinear.NonlinearOptions(sample.key),
        )


def test_row_selection_and_standard_air_probe_conversion(sample):
    options = nonlinear.NonlinearOptions(sample.key, rows=("1", "3"), factor=2.0)
    obs = nonlinear.gamma_observation(sample, options)
    # Wavelengths round-trip through SI metres.
    np.testing.assert_allclose(obs.wavelength_nm, [500, 800], rtol=2e-15, atol=0)
    expected = np.array([float(sample.values[i]) for i in (1, 3)]) * 2
    np.testing.assert_allclose(obs.gamma, expected, rtol=2e-15, atol=0)
    np.testing.assert_allclose(obs.sigma, expected * 0.001, rtol=2e-15, atol=0)
    with pytest.raises(DataError, match="unknown rows"):
        nonlinear.gamma_observation(
            sample, dataclasses.replace(options, rows=("missing",))
        )
    air = dataclasses.replace(sample, wavelength_medium="standard_air")
    assert nonlinear._wavenumber(air, 500) == pytest.approx(
        1e7 / (500 * 1.000277), rel=2e-15
    )
    with pytest.raises(DataError, match="no formula"):
        nonlinear._formula("missing/formula")


@pytest.mark.parametrize(
    ("entry", "message"),
    [
        (1, "must be a key or a table"),
        ({"key": "missing", "typo": 1}, "unknown keys"),
        ("missing", "unknown dataset"),
    ],
)
def test_nonlinear_dataset_options_are_validated(entry, message):
    with pytest.raises(DataError, match=message):
        nonlinear._options(entry, "synthetic")


def test_synthetic_fit_recovers_dispersion_with_a_systematic_scale_prior(
    sample, monkeypatch, tmp_path
):
    monkeypatch.setattr(nonlinear, "datasets", lambda: (sample,))
    path = tmp_path / "fit.toml"
    text = (
        'gas="He"\nrange_nm=[400,1500]\n'
        'fitted=[{key="synthetic/gamma",systematic=0.02}]\n'
        'validators=[{key="synthetic/gamma",rows=["missing"]}]\n'
    )
    path.write_text(text, encoding="utf-8")
    config = nonlinear.load_nonlinear_config(path)
    result = nonlinear.fit_hyperpolarizability(config)
    # Exact synthetic observations; solver and finite-difference conditioning set 1e-6.
    np.testing.assert_allclose(result.parameters, [10, 2, 3, 0], rtol=1e-6, atol=1e-8)
    assert result.scales[sample.key] == pytest.approx(1.0, abs=1e-7)
    record = json.loads(nonlinear.dumps(result))
    assert record["config_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["parameters"]["B_prime_1e-20_cm4"] == 0
    assert "unknown rows" in record["validation"][sample.key]["note"]
    assert np.shape(record["covariance"]) == (4, 4)
    assert record["max_w2_fitted_cm-2"] == pytest.approx(
        4 * (1e7 / 400) ** 2, rel=2e-15
    )
    unstated = dataclasses.replace(sample, uncertainty=Uncertainty())
    monkeypatch.setattr(nonlinear, "datasets", lambda: (unstated,))
    with pytest.raises(DataError, match="fitted value needs an uncertainty"):
        nonlinear.fit_hyperpolarizability(config)
    path.write_text("typo=1\n" + text, encoding="utf-8")
    with pytest.raises(DataError, match="unknown keys"):
        nonlinear.load_nonlinear_config(path)


def test_range_endpoints_and_static_limit_do_not_warn_but_extrapolation_does():
    model = nonlinear.HyperpolarizabilityModel(
        Gas.He, (10.0, 0.0, 0.0, 0.0), np.zeros((4, 4)), wavelength_range=(400e-9, 1e-6)
    )
    np.testing.assert_array_equal(
        model.gamma_au("kerr", np.array([400e-9, 1e-6, np.inf])), [10, 10, 10]
    )
    with pytest.warns(ModelRangeWarning, match="outside"):
        assert model.gamma_au("kerr", 399e-9) == 10

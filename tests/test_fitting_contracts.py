# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Small synthetic fitting problems with independently calculable answers."""

import dataclasses
import hashlib
import json
import math

import numpy as np
import pytest

import optmat.gas.models.fitted as stored
from optmat.constants import AU_ANGULAR_FREQUENCY, C
from optmat.datasets import Uncertainty, load_source
from optmat.exceptions import DataError
from optmat.gas import fitting
from optmat.gas.fitting import DatasetOptions, FitConfig, Line, Observation, fit
from optmat.gas.gases import Gas
from optmat.gas.models import registry


@pytest.fixture
def config():
    # The only free physical parameter is the strength of a pole at one hartree.
    return FitConfig(
        name="synthetic",
        gas=Gas.He,
        description="synthetic",
        lines=(Line("one", fitting.HARTREE_NM, 2.0, 0.5),),
        effective=0,
        effective_min_energy=2.0,
        sum_priors={},
        chi0=0.0,
        fitted=(),
        validators=(),
        range_nm=(100.0, math.inf),
        digest="a" * 64,
    )


@pytest.fixture
def sample(tmp_path):
    path = tmp_path / "synthetic.toml"
    path.write_text(
        'citation="synthetic"\nkind="theory"\n[[datasets]]\nid="alpha"\n'
        'material="He"\nquantity="polarizability"\nunit="au"\n'
        'abscissa_unit="nm"\npoints=[["a","500","2.2","0.1"],'
        '["b","1000","2.3","0.2"]]\nuncertainty={column=true}\n',
        encoding="utf-8",
    )
    return load_source(path).datasets[0]


def test_weighted_fit_and_covariance_equal_the_linear_normal_equations(config):
    omega = np.array([0.0, 0.2, 0.4, 0.6])
    y = np.array([2.03, 2.13, 2.32, 3.10])
    sigma = np.array([0.03, 0.05, 0.02, 0.04])
    data = Observation(
        "synthetic/linear", omega, y, sigma, None, DatasetOptions("synthetic/linear")
    )
    result = fit(config, [data], starts=1)
    x = 1 / (1 - omega**2)
    information = np.sum((x / sigma) ** 2) + 1 / 0.5**2
    expected = (np.sum(x * y / sigma**2) + 2 / 0.5**2) / information
    # Solver convergence and its finite-difference Jacobian: relative error < 1e-7.
    assert result.model.strengths == pytest.approx((expected,), rel=1e-7, abs=0.0)
    assert result.covariance[0, 0] == pytest.approx(1 / information, rel=1e-7, abs=0.0)
    expected_sigma = x / np.sqrt(information)
    np.testing.assert_allclose(
        result.polarizability_sigma_au(omega, birge=False),
        expected_sigma,
        rtol=1e-7,
        atol=0,
    )
    scaled = dataclasses.replace(result, chi2=12.0, dof=3)
    np.testing.assert_allclose(
        scaled.polarizability_sigma_au(omega), 2 * expected_sigma, rtol=1e-7, atol=0
    )
    assert result.dataset_chi2[data.key][0] == 4


def test_systematic_scale_prior_and_fixed_sum_terms_enter_residuals(config):
    # Test the residual vector directly: this isolates the likelihood, before fitting.
    config = dataclasses.replace(
        config, sum_priors={-2: (2.25, 0.1)}, fixed=(("fixed", 1.0, 2.0),)
    )
    obs = Observation(
        "test",
        np.array([0.0]),
        np.array([2.5]),
        np.array([0.2]),
        0.1,
        DatasetOptions("test"),
    )
    residual = fitting._residuals(
        np.array([2.0, 0.2]), config, fitting._Layout(1, 0, (0,)), [obs]
    )
    # Model = 2 + 1/4; scaled by 1.2; the line and sum priors are exactly met.
    np.testing.assert_allclose(residual, [1.0, 0.0, 2.0, 0.0], rtol=2e-15, atol=2e-15)
    free = dataclasses.replace(obs, options=DatasetOptions("test", free_scale=True))
    np.testing.assert_allclose(
        fitting._residuals(
            np.array([2.0, 0.2]), config, fitting._Layout(1, 0, (0,)), [free]
        ),
        [1.0, 0.0, 0.0],
        rtol=2e-15,
        atol=2e-15,
    )


def test_observation_window_shape_uncertainty_and_missing_reference(
    sample, config, monkeypatch
):
    monkeypatch.setattr(fitting, "datasets", lambda: (sample,))
    options = DatasetOptions(
        sample.key,
        range_nm=(750, 1250),
        free_scale=True,
        shape_relative=0.01,
        extra_relative=0.02,
    )
    (obs,) = fitting.observations(dataclasses.replace(config, fitted=(options,)))
    assert obs.alpha == pytest.approx([2.3], rel=2e-15, abs=0.0)
    assert obs.sigma == pytest.approx(
        [2.3 * math.hypot(0.01, 0.02)], rel=2e-15, abs=0.0
    )
    with pytest.raises(DataError, match="shape_relative needs free_scale"):
        fitting._observation(sample, dataclasses.replace(options, free_scale=False))
    with pytest.raises(DataError, match="reference_temperature_k"):
        fitting._observation(
            sample, DatasetOptions(sample.key, temperature=(300, 0.01, 0.001))
        )
    unstated = dataclasses.replace(sample, uncertainty=Uncertainty())
    with pytest.raises(DataError, match="needs a stated uncertainty"):
        fitting._observation(unstated, DatasetOptions(sample.key))


def test_validator_report_retains_unsupported_quantities_and_filters_moments(
    sample, config, monkeypatch
):
    moment = dataclasses.replace(
        sample,
        id="moments",
        quantity="cauchy_moment",
        abscissa="order",
        abscissa_unit="1",
        abscissae=("-2", "0"),
        values=("2.2", "99"),
    )
    unsupported = dataclasses.replace(
        sample,
        id="unsupported",
        quantity="hyperpolarizability",
        parameters={"process": "kerr", "convention": "T", "part": "electronic"},
    )
    config = dataclasses.replace(
        config,
        validators=tuple(DatasetOptions(d.key) for d in (sample, moment, unsupported)),
    )
    model = fitting.OscillatorModel(
        gas=Gas.He,
        name="synthetic",
        citation="synthetic",
        wavelength_range=(100e-9, math.inf),
        kind="fitted",
        strengths=(2.0,),
        energies=(1.0,),
        labels=("one",),
    )
    result = fitting.FitResult(config, model, np.array([2.0]), np.array([[0.01]]))
    monkeypatch.setattr(fitting, "datasets", lambda: (sample, moment, unsupported))
    report = fitting.validation(result)
    assert report[moment.key]["n"] == 1
    assert report[moment.key]["chi2"] == pytest.approx(4.0, rel=3e-15, abs=0.0)
    assert report[moment.key]["mean_relative_deviation"] == pytest.approx(
        0.1, abs=2e-15
    )
    assert "hyperpolarizability" in str(report[unsupported.key]["note"])
    assert report[sample.key]["n"] == 2
    with pytest.raises(DataError, match="fitted sum needs a stated uncertainty"):
        fitting._moments(
            dataclasses.replace(moment, uncertainty=Uncertainty()),
            DatasetOptions(moment.key),
        )


CONFIG = (
    'gas="He"\nrange_nm=[100,"inf"]\n[model]\n'
    "effective=1\neffective_min_energy_hartree=1\n"
)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("typo=1\n" + CONFIG, "unknown top-level keys"),
        (CONFIG + '[[fitted]]\nkey="synthetic/alpha"\ntypo=1\n', "unknown keys"),
        (CONFIG + '[[fitted]]\nkey="missing"\n', "unknown dataset"),
        ('validators=["missing"]\n' + CONFIG, "unknown validator"),
        (CONFIG.replace('[100,"inf"]', "[100]"), r"expected \[low, high\]"),
        (
            CONFIG + '[[fitted]]\nkey="synthetic/alpha"\ntemperature={}\n',
            "temperature needs exactly",
        ),
        (CONFIG + "fixed_intervals={}\n", "fixed_intervals needs a dataset"),
        (
            CONFIG + 'fixed_intervals={dataset="synthetic/alpha"}\n',
            "not an interval table",
        ),
    ],
)
def test_invalid_fit_configuration_is_rejected_before_optimisation(
    tmp_path, monkeypatch, sample, text, message
):
    monkeypatch.setattr(fitting, "datasets", lambda: (sample,))
    path = tmp_path / "fit.toml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DataError, match=message):
        fitting.load_config(path)


def test_serialised_fit_preserves_provenance_prediction_and_uncertainty(
    tmp_path, sample, monkeypatch
):
    monkeypatch.setattr(fitting, "datasets", lambda: (sample,))
    path = tmp_path / "synthetic.toml"
    path.write_text(CONFIG + '[[fitted]]\nkey="synthetic/alpha"\n', encoding="utf-8")
    config = fitting.load_config(path)
    result = fit(config, starts=1)
    record = json.loads(fitting.dumps(result))
    assert record["config_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["gas"] == "He"
    assert record["scales"] == {sample.key: 1.0}
    assert record["n_fixed"] == 0
    assert record["dataset_chi2"][sample.key][0] == 2
    root = tmp_path / "models" / "fitted"
    root.mkdir(parents=True)
    (root / "synthetic.json").write_text(fitting.dumps(result), encoding="utf-8")
    monkeypatch.setattr(stored, "GAS_DATA_DIRECTORY", tmp_path)
    stored.fitted_models.cache_clear()
    try:
        rebuilt = stored.fitted_models()["fit/synthetic"]
        omega = np.array([0.0, 0.02, 0.05])
        np.testing.assert_array_equal(
            rebuilt.polarizability_au(omega), result.model.polarizability_au(omega)
        )
        np.testing.assert_array_equal(
            rebuilt.polarizability_sigma_au(omega),
            result.polarizability_sigma_au(omega),
        )
    finally:
        stored.fitted_models.cache_clear()


def test_ratio_validators_keep_orientation_and_require_uncertainty(
    sample, config, monkeypatch
):
    ratio = dataclasses.replace(
        sample,
        quantity="phase_match_density_ratio",
        unit="1",
        parameters={"reference_material": "Ar", "harmonic": 2},
    )
    with pytest.raises(DataError, match="neither gas is"):
        fitting._ratio(ratio, DatasetOptions(ratio.key), Gas.Ne)
    with pytest.raises(DataError, match="ratio needs a stated uncertainty"):
        fitting._ratio(
            dataclasses.replace(ratio, uncertainty=Uncertainty()),
            DatasetOptions(ratio.key),
            Gas.He,
        )
    monkeypatch.setattr(fitting, "datasets", lambda: (ratio,))
    result = fit(
        config,
        [
            Observation(
                "synthetic",
                np.array([0.0]),
                np.array([2.0]),
                np.array([0.1]),
                None,
                DatasetOptions("synthetic"),
            )
        ],
        starts=1,
    )
    result = dataclasses.replace(
        result,
        config=dataclasses.replace(config, validators=(DatasetOptions(ratio.key),)),
    )
    partner = dataclasses.replace(
        result.model, gas=Gas.Ar, strengths=(5.0,), energies=(2.0,)
    )
    monkeypatch.setattr(registry, "default_model", {Gas.Ar: partner}.__getitem__)
    # Independently evaluate the two invented poles: alpha_He=2/(1-w^2),
    # alpha_Ar=5/(4-w^2). Reverse the recorded gases to check both ratio orientations.
    w = 2 * np.pi * C / np.array([500e-9, 1000e-9]) / AU_ANGULAR_FREQUENCY
    partner_delta = 5 / (4 - (2 * w) ** 2) - 5 / (4 - w**2)
    fitted_delta = 2 / (1 - (2 * w) ** 2) - 2 / (1 - w**2)
    for reverse in (False, True):
        ratio = dataclasses.replace(
            ratio,
            material="Ar" if reverse else "He",
            parameters={"reference_material": "He" if reverse else "Ar", "harmonic": 2},
        )
        report = fitting.validation(result)[ratio.key]
        assert report["n"] == 2
        predicted = (
            fitted_delta / partner_delta if reverse else partner_delta / fitted_delta
        )
        expected = np.array([2.2, 2.3]) / predicted - 1
        # Pole subtraction and the one-parameter fit introduce rounding only.
        assert report["mean_relative_deviation"] == pytest.approx(
            np.mean(expected), rel=1e-12, abs=0.0
        )


def test_rayleigh_observation_dispatch_checks_uncertainty(sample):
    from optmat.conditions import State  # noqa: PLC0415

    rayleigh = dataclasses.replace(
        sample,
        quantity="rayleigh_cross_section",
        unit="cm2",
        state=State(101325, 300),
        parameters={"king_factor": 1.0},
        values=("1e-27", "2e-27"),
    )
    obs = fitting._observation(rayleigh, DatasetOptions(rayleigh.key))
    assert obs.options.electric
    assert np.all(obs.alpha > 0)
    with pytest.raises(DataError, match="cross section needs a stated uncertainty"):
        fitting._observation(
            dataclasses.replace(rayleigh, uncertainty=Uncertainty()),
            DatasetOptions(rayleigh.key),
        )


def test_measured_and_fixed_electronic_terms_enter_moments_and_serialisation(config):
    from optmat.constants import (  # noqa: PLC0415
        AU_ENERGY,
        E_CHARGE,
        EPS0,
        HBAR,
        M_E,
        C,
    )
    from optmat.gas.absorption import MeasuredAbsorption, Piece  # noqa: PLC0415

    grid = np.array([0.0, 10000.0, 20000.0])
    cross_section = np.array([0.0, 1e-22, 1e-22])
    measured = MeasuredAbsorption(
        Gas.He,
        300.0,
        (Piece("synthetic", 500, 1000, "test"),),
        0.0,
        10000.0,
        grid,
        cross_section,
        np.zeros(3),
    )
    # S(0) is the rectangular df/dE area between the positive grid endpoints.
    density = 1e-22 * 2 * EPS0 * M_E * C * AU_ENERGY / (np.pi * E_CHARGE**2 * HBAR)
    width = 10000 * 100 * 2 * np.pi * C * HBAR / AU_ENERGY
    area = density * width
    config = dataclasses.replace(
        config,
        absorption=measured,
        fixed=(("fixed", 1.0, 2.0),),
        sum_priors={0: (3 + area, 0.1)},
    )
    obs = Observation(
        "moment",
        np.array([0.0]),
        np.array([3 + area]),
        np.array([0.1]),
        None,
        DatasetOptions("moment"),
        orders=np.array([0.0]),
    )
    result = fit(config, [obs], starts=1)
    assert result.model.strengths == pytest.approx((2, 1), rel=1e-7, abs=0.0)
    assert result.chi2 < 1e-12  # exact moment and prior, up to optimiser rounding
    record = fitting.result_record(result)
    assert record["absorption"] == {
        "fwhm_cm": 0,
        "step_cm": 10000,
        "temperature_k": 300,
        "pieces": [["synthetic", 500, 1000]],
    }

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Stored nonlinear fits agree at the solver's numerical convergence scale."""

from collections.abc import Mapping
from dataclasses import replace
from typing import cast

import numpy as np
import pytest

from optmat.gas.gases import Gas
from optmat.gas.nonlinear import (
    NONLINEAR_DIRECTORY,
    NonlinearFit,
    _dispersion,
    _observations,
    _scale_groups,
    fit_hyperpolarizability,
    frequency_sums,
    hyperpolarizability_models,
    load_nonlinear_config,
)

GASES = sorted(p.stem for p in NONLINEAR_DIRECTORY.glob("*.toml"))
# least_squares uses default ftol/xtol/gtol=1e-8. Near a least-squares minimum,
# parameter resolution scales as sqrt(objective resolution), in units of sigma.
# Finite-difference, three-point and analytic-Jacobian convergence experiments gave
# <7e-6 sigma parameter displacement and <1.4e-8 relative prediction differences
# at fitted points. Use 1e-4 sigma and 1e-7 predictions, also over the declared range.
# These are numerical regression bounds, not estimates of scientific uncertainty.
PARAMETER_SIGMA_TOL = 1e-4
PREDICTION_RTOL = 1e-7
OBJECTIVE_TOL = 1e-8


def _assert_parameters(actual, expected, covariance):
    sigma = np.sqrt(np.diag(covariance))
    assert np.all(np.isfinite(actual)), "nonfinite parameters"
    assert np.all(np.abs(actual - expected) <= PARAMETER_SIGMA_TOL * sigma), (
        "parameter displacement exceeds numerical convergence scale"
    )
    # A fixed parameter has zero variance, hence must agree exactly above.


def _objective(config, parameters, scales):
    """Reconstruct data and prior contributions, counting each scale prior once."""
    data = _observations(config.fitted)
    groups, membership, group_prior = _scale_groups(data)
    assert set(scales) == {obs.key for obs in data}, "dataset selection"
    group_scales = {}
    data_chi2 = 0.0
    for obs, group in zip(data, membership, strict=True):
        scale = scales[obs.key]
        if group < 0:
            assert scale == 1.0, "fixed dataset scale"
        else:
            previous = group_scales.setdefault(group, scale)
            assert scale == previous, "shared dataset scale"
        residuals = (
            scale * _dispersion(parameters, obs.w2, obs.w4) - obs.gamma
        ) / obs.sigma
        data_chi2 += float(residuals @ residuals)
    assert len(group_scales) == len(groups)
    priors = sum(
        ((parameters[position] - prior[0]) / prior[1]) ** 2
        for position, prior in (
            (1, config.a_prior),
            (2, config.b_prior),
            (3, config.b_prime_prior),
        )
        if prior is not None
    )
    priors += sum(
        ((group_scales[group] - 1) / sigma) ** 2 for group, sigma in group_prior.items()
    )
    return data_chi2, data_chi2 + priors


def _assert_objective(actual, expected):
    assert np.isfinite(actual), "nonfinite objective"
    assert abs(actual - expected) <= OBJECTIVE_TOL * max(1.0, expected), (
        "full objective differs beyond solver convergence"
    )


def _assert_predictions(config, actual, expected):
    probes = [(obs.w2, obs.w4) for obs in _observations(config.fitted)]
    # Independent process geometries constrain B and B', even where the fitted
    # observations alone cannot distinguish them. Include static and range endpoints.
    low, high = config.range_nm
    wavenumbers = 1e7 / np.geomspace(low, high, 41)  # cm^-1
    probes.extend(
        frequency_sums(process, wavenumbers)
        for process in ("static", "dc_kerr", "kerr", "eshg", "thg")
    )
    for w2, w4 in probes:
        np.testing.assert_allclose(
            _dispersion(actual, w2, w4),
            _dispersion(expected, w2, w4),
            rtol=PREDICTION_RTOL,
            atol=0,
            err_msg="model predictions",
        )


def _assert_fit(config, stored, fresh):
    record = stored.record
    assert record["config_sha256"] == config.digest, "configuration provenance"
    assert fresh.config.digest == config.digest, "fresh configuration provenance"
    assert fresh.dof == record["dof"], "degrees of freedom"
    assert fresh.dataset_chi2.keys() == record["dataset_chi2"].keys(), (
        "dataset selection"
    )
    for key, (count, _) in fresh.dataset_chi2.items():
        assert count == record["dataset_chi2"][key][0], "dataset point count"
    _assert_parameters(fresh.parameters, stored.parameters, stored.covariance)
    _assert_predictions(config, fresh.parameters, stored.parameters)
    assert fresh.scales.keys() == record["scales"].keys(), "dataset selection"
    np.testing.assert_allclose(
        [fresh.scales[k] for k in fresh.scales],
        [record["scales"][k] for k in fresh.scales],
        rtol=PREDICTION_RTOL,
        atol=0,
        err_msg="dataset scales",
    )
    old_data_chi2, old_objective = _objective(
        config, stored.parameters, record["scales"]
    )
    new_data_chi2, new_objective = _objective(config, fresh.parameters, fresh.scales)
    # Statistics must still describe the parameters recorded beside them; 1e-12
    # allows only arithmetic rounding, not solver drift, in these recomputations.
    np.testing.assert_allclose(old_data_chi2, record["chi2"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(new_data_chi2, fresh.chi2, rtol=1e-12, atol=1e-12)
    _assert_objective(new_objective, old_objective)


@pytest.mark.parametrize("gas", GASES)
def test_stored_hyperpolarizability_matches_a_fresh_fit(gas):
    config = load_nonlinear_config(gas)
    stored = hyperpolarizability_models()[Gas(gas)]
    _assert_fit(config, stored, fit_hyperpolarizability(config))


@pytest.fixture
def reference_fit():
    config = load_nonlinear_config("CH4")
    stored = hyperpolarizability_models()[Gas.CH4]
    record = stored.record
    return (
        config,
        stored,
        NonlinearFit(
            config=config,
            parameters=np.asarray(stored.parameters),
            covariance=stored.covariance,
            scales=cast(Mapping[str, float], record["scales"]),
            chi2=cast(float, record["chi2"]),
            dof=cast(int, record["dof"]),
            dataset_chi2={
                key: (int(values[0]), values[1])
                for key, values in cast(
                    Mapping[str, list[float]], record["dataset_chi2"]
                ).items()
            },
        ),
    )


def test_small_parameter_roundoff_is_accepted(reference_fit):
    config, stored, fresh = reference_fit
    parameters = np.nextafter(fresh.parameters, np.inf)
    _assert_fit(config, stored, replace(fresh, parameters=parameters))


@pytest.mark.parametrize("position", range(4))
def test_sub_uncertainty_parameter_regressions_are_rejected(reference_fit, position):
    config, stored, fresh = reference_fit
    parameters = fresh.parameters.copy()
    parameters[position] += 1e-3 * np.sqrt(stored.covariance[position, position])
    with pytest.raises(AssertionError, match="parameter displacement"):
        _assert_fit(config, stored, replace(fresh, parameters=parameters))


def test_fixed_parameter_cannot_move():
    with pytest.raises(AssertionError, match="parameter displacement"):
        _assert_parameters(np.array([1e-30]), np.array([0.0]), np.zeros((1, 1)))


def test_process_predictions_detect_a_fit_degeneracy(reference_fit):
    config, stored, _ = reference_fit
    # CH4's ESHG observations constrain B + B'/2. Moving along that degeneracy
    # preserves those predictions but must fail in the other sampled processes.
    parameters = np.asarray(stored.parameters).copy()
    parameters[2] -= 0.005
    parameters[3] += 0.01
    for obs in _observations(config.fitted):
        np.testing.assert_allclose(
            _dispersion(parameters, obs.w2, obs.w4),
            _dispersion(stored.parameters, obs.w2, obs.w4),
            rtol=1e-12,
            atol=0,
        )
    with pytest.raises(AssertionError, match="model predictions"):
        _assert_predictions(config, parameters, stored.parameters)


def test_dataset_scale_regression_is_rejected(reference_fit):
    config, stored, fresh = reference_fit
    scales = dict(fresh.scales)
    scales[next(iter(scales))] *= 1.001
    with pytest.raises(AssertionError, match="dataset scales"):
        _assert_fit(config, stored, replace(fresh, scales=scales))


def test_changed_dataset_selection_is_rejected(reference_fit):
    config, stored, fresh = reference_fit
    with pytest.raises(AssertionError, match="dataset selection"):
        _assert_fit(config, stored, replace(fresh, dataset_chi2={}))


def test_parameter_prior_enters_full_objective(reference_fit):
    config, stored, fresh = reference_fit
    data, total = _objective(config, stored.parameters, fresh.scales)
    assert total - data == pytest.approx(
        (stored.parameters[3] / config.b_prime_prior[1]) ** 2, abs=1e-15
    )
    # Changing the prior penalty alone must affect the regression, even when the
    # model predictions and data residuals are untouched.
    changed = replace(config, b_prime_prior=(1.0, 1.3))
    new_data, new_total = _objective(changed, stored.parameters, fresh.scales)
    assert new_data == data
    with pytest.raises(AssertionError, match="full objective"):
        _assert_objective(new_total, total)


def test_shared_scale_prior_is_counted_once():
    config = load_nonlinear_config("Ar")
    stored = hyperpolarizability_models()[Gas.Ar]
    # Ar's three calculation datasets share one freely fitted scale. Give that
    # group a prior in this test only: it must add one penalty, not three.
    group = [o for o in config.fitted if o.scale_group == "hattig_ccsd"]
    assert len(group) == 3
    constrained = replace(
        config,
        fitted=tuple(
            replace(o, free_scale=False, systematic=0.1) if o in group else o
            for o in config.fitted
        ),
    )
    scales = cast(Mapping[str, float], stored.record["scales"])
    data, total = _objective(config, stored.parameters, scales)
    constrained_data, constrained_total = _objective(
        constrained, stored.parameters, scales
    )
    assert constrained_data == data
    penalty = ((scales[group[0].key] - 1) / 0.1) ** 2
    assert constrained_total - total == pytest.approx(penalty, rel=1e-12, abs=1e-12)
    with pytest.raises(AssertionError, match="full objective"):
        _assert_objective(constrained_total, total)


def test_shared_scales_cannot_diverge():
    config = load_nonlinear_config("Ar")
    stored = hyperpolarizability_models()[Gas.Ar]
    scales = dict(cast(Mapping[str, float], stored.record["scales"]))
    key = next(o.key for o in config.fitted if o.scale_group == "hattig_ccsd")
    scales[key] *= 1.001
    with pytest.raises(AssertionError, match="shared dataset scale"):
        _objective(config, stored.parameters, scales)


def test_fixed_scales_cannot_move(reference_fit):
    config, stored, fresh = reference_fit
    scales = dict(fresh.scales)
    scales[next(iter(scales))] *= 1.001
    with pytest.raises(AssertionError, match="fixed dataset scale"):
        _objective(config, stored.parameters, scales)


def test_changed_configuration_digest_is_rejected(reference_fit):
    config, stored, fresh = reference_fit
    with pytest.raises(AssertionError, match="configuration provenance"):
        _assert_fit(replace(config, digest="0" * 64), stored, fresh)

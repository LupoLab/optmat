# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Every stored fit is what its configuration gives today."""

import copy
import json
import math
from dataclasses import replace

import numpy as np
import pytest

from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.fitting import (  # pyright: ignore[reportPrivateUsage]
    FIT_DIRECTORY,
    HARTREE_NM,
    _Layout,
    _residuals,
    fit,
    load_config,
    observations,
    result_record,
)
from optmat.gas.models.fitted import fitted_models

CONFIGS = sorted(p.stem for p in FIT_DIRECTORY.glob("*.toml"))


# A fresh refit against the stored record (TC-1090). The effective oscillators are
# nearly degenerate (a far one enters only through f / E^2), so the solver's stopping
# point depends on its path: refits of the committed code move single parameters by
# up to 3e7 of their stored sigma (SF6) while alpha within the fit's intervals moves
# by at most 1.1 % of its own uncertainty band (He) or 1.5e-5 relative (Ne), and the
# full objective, priors included, by at most 1.1e-3 (Ne). The checks are on those
# predictions and on the objective, at 4-9 times those values; a change of 0.01 in
# chi-squared is a 0.1-sigma move along any direction.
BAND_FRACTION = 0.05
ALPHA_RELATIVE = 1e-4
OBJECTIVE_CHI2 = 1e-2


def _objective(config, record):
    """Return the fit's full objective (data and priors) at a record's parameters."""
    data = observations(config)
    scaled = tuple(
        i for i, obs in enumerate(data) if obs.options.free_scale or obs.systematic
    )
    layout = _Layout(len(config.lines), config.effective, scaled)
    n_fitted = len(record["strengths"]) - len(config.fixed_terms)
    strengths = np.asarray(record["strengths"])[:n_fitted]
    energies = np.asarray(record["energies_hartree"])[:n_fitted]
    n = layout.n_lines
    excess = np.maximum(energies[n:] - config.effective_min_energy, 1e-12)
    offsets = [record["scales"].get(data[i].key, 1.0) - 1.0 for i in layout.scaled]
    theta = np.concatenate(
        [strengths[:n], np.log(strengths[n:]), np.log(excess), offsets]
    )
    residual = _residuals(theta, config, layout, data)
    return float(residual @ residual)


def _assert_refit(config, stored, record, fresh, fresh_record):
    """Fail unless a fresh fit agrees with the stored one within the solver's reach."""
    assert fresh_record["dof"] == record["dof"]
    assert {k: v[0] for k, v in fresh_record["dataset_chi2"].items()} == {
        k: v[0] for k, v in record["dataset_chi2"].items()
    }
    for interval in stored.intervals:
        omega = _interval_grid(interval, 200)
        old = stored.polarizability_au(omega)
        new = fresh.polarizability_au(omega)
        band = stored.polarizability_sigma_au(omega)
        assert np.max(np.abs(new - old) / band) <= BAND_FRACTION
        np.testing.assert_allclose(new, old, rtol=ALPHA_RELATIVE, atol=0)
    difference = _objective(config, fresh_record) - _objective(config, record)
    assert abs(difference) <= OBJECTIVE_CHI2, difference


@pytest.mark.fits
@pytest.mark.parametrize("name", CONFIGS)
def test_stored_fit_matches_a_fresh_fit(name):
    stored = fitted_models()[f"fit/{name}"]
    config = load_config(name)
    result = fit(config)
    _assert_refit(config, stored, stored.record, result.model, result_record(result))


def _perturbed(record, index, factor):
    moved = copy.deepcopy(record)
    strengths = list(moved["strengths"])
    strengths[index] *= factor
    moved["strengths"] = strengths
    return moved


@pytest.fixture(scope="module")
def nitrogen():
    stored = fitted_models()["fit/N2"]
    return load_config("N2"), stored, dict(stored.record)


def test_a_refit_equal_to_the_record_is_accepted(nitrogen):
    config, stored, record = nitrogen
    _assert_refit(config, stored, record, stored, record)


@pytest.mark.parametrize(
    "change",
    [
        lambda r: _perturbed(r, 0, 1 + 1e-4),  # alpha moves by about 1e-4
        lambda r: {**r, "dof": r["dof"] + 1},
        lambda r: {**r, "dataset_chi2": dict(list(r["dataset_chi2"].items())[1:])},
    ],
    ids=["strength", "dof", "selection"],
)
def test_a_real_change_in_a_refit_is_rejected(nitrogen, change):
    config, stored, record = nitrogen
    moved = change(record)
    fresh = replace(stored, strengths=tuple(moved["strengths"]))
    with pytest.raises(AssertionError):
        _assert_refit(config, stored, record, fresh, moved)


@pytest.mark.parametrize("name", CONFIGS)
def test_stored_fit_records_its_configuration(name):
    stored = fitted_models()[f"fit/{name}"]
    assert stored.record["config_sha256"] == load_config(name).digest


@pytest.mark.parametrize("name", CONFIGS)
def test_stored_fit_is_a_physical_polarizability_over_its_range(name):
    # Per qualified interval (OEP 3): every stored fit has one today, its range.
    stored = fitted_models()[f"fit/{name}"]
    for interval in stored.intervals:
        _check_interval(stored, _interval_grid(interval, 400))


def _check_interval(stored, omega):
    """Check positivity, passivity and monotony between poles over one interval."""
    alpha = stored.polarizability_au(omega)
    if stored.absorbing:
        # Inside its absorption bands the index has anomalous dispersion: only
        # positivity and passivity (Im alpha >= 0) hold.
        wavelength = HARTREE_NM * 1e-9 / omega[omega > 0]
        complex_alpha = stored.complex_polarizability(wavelength, warn=False)
        assert np.all(np.real(complex_alpha) > 0)
        assert np.all(np.imag(complex_alpha) >= 0)
        return
    # Between its poles a sum of positive-strength oscillators rises monotonically
    # (the infrared bands of CH4, SF6, N2O lie below the declared range, the electronic
    # terms above it); a pole inside the range would break it. Positive because the
    # electronic terms dominate the small negative infrared tail.
    energies = np.asarray(stored.energies)
    assert np.all(np.asarray(stored.strengths) > 0)
    assert not np.any((energies > omega[0]) & (energies < omega[-1]))
    assert np.all(alpha > 0)
    assert np.all(np.diff(alpha) > 0)


def _interval_grid(interval, points):
    """Photon energies (hartree) across one qualified interval, m."""
    shortest, longest = (x * 1e9 for x in interval)
    low = 0.0 if math.isinf(longest) else 1.001 * HARTREE_NM / longest
    return np.linspace(low, 0.999 * HARTREE_NM / shortest, points)


def _range_grid(stored, points):
    """Photon energies (hartree) across a model's declared wavelength range."""
    return _interval_grid(stored.wavelength_range, points)


@pytest.mark.parametrize("name", CONFIGS)
def test_stored_fit_uncertainty_is_positive_over_its_range(name):
    stored = fitted_models()[f"fit/{name}"]
    covariance = np.asarray(stored.covariance)
    scale = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    # A parameter at its bound (SF6's lower effective energy) has zero variance in the
    # log parametrisation: its row and column are zero, symmetric by construction.
    scale = np.where(scale > 0, scale, 1.0)
    # Symmetric to round-off: a far effective oscillator enters only through f / E^2,
    # so (f, E) is degenerate and the covariance singular to 1e-17 of its scale;
    # its asymmetry (at most 4e-6 of the scale, helium) is that round-off.
    assert np.max(np.abs(covariance - covariance.T) / scale) < 1e-4
    # The degenerate direction does not move alpha: the propagated variance is
    # positive wherever the model is valid.
    omega = _range_grid(stored, 200)
    sigma = stored.polarizability_sigma_au(omega)
    assert np.all(np.isfinite(sigma))
    assert np.all(sigma > 0)


@pytest.mark.parametrize("name", CONFIGS)
def test_robustness_record_is_current_and_complete(name):
    path = GAS_DATA_DIRECTORY / "models" / "robustness" / f"{name}.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    config = load_config(name)
    # Regenerate with tools/fit/robustness.py after changing a fit.
    assert record["config_sha256"] == config.digest
    assert list(record["dropped"]) == [options.key for options in config.fitted]

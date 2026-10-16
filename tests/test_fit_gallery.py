# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The gallery compares data in the same conventions and windows as the fitter."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from matplotlib.figure import Figure

from optmat.constants import AU_POLARIZABILITY
from optmat.datasets import datasets
from optmat.gas import reduction
from optmat.gas.fitting import load_config, observations
from optmat.gas.models.fitted import fitted_models
from optmat.gas.models.oscillator import angular_frequency_au


@pytest.fixture(scope="module")
def helpers():
    path = Path(__file__).parents[1] / "docs" / "gallery" / "_helpers.py"
    spec = importlib.util.spec_from_file_location("fit_gallery_helpers", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(params=[name.removeprefix("fit/") for name in fitted_models()])
def plotted(request, helpers):
    model = fitted_models()[f"fit/{request.param}"]
    config = load_config(request.param)
    figure = Figure()
    top, bottom = figure.subplots(2)
    helpers._fit_data(top, bottom, model, config, (100e-9, 20e-6))
    yield model, config, top, bottom
    figure.clear()


def test_fitted_residuals_reproduce_stored_dataset_chi_squared(plotted, helpers):
    model, config, top, bottom = plotted
    plotted_keys = {
        container.get_label().removesuffix(" (fitted)")
        for container in top.containers
        if container.get_label().endswith(" (fitted)")
    }
    residuals = {line.get_label(): line for line in bottom.lines}
    eligible = {d.key for d in helpers.plottable(config.gas.value)}
    expected_keys = eligible & {option.key for option in config.fitted}
    assert residuals.keys() == plotted_keys == expected_keys
    for key, line in residuals.items():
        count, chi_squared = model.record["dataset_chi2"][key]
        assert len(line.get_ydata()) == count
        # Stored statistics and plotting differ only in evaluation round-off, which
        # subtraction amplifies for the most precise data (helium's metrology).
        assert np.sum(line.get_ydata() ** 2) == pytest.approx(
            chi_squared, rel=1e-6, abs=1e-7
        )
        options = next(option for option in config.fitted if option.key == key)
        shown = np.asarray(line.get_xdata())
        finite = shown != 18000.0  # static observations are displayed at 0.9 x xlim[1]
        assert np.all(shown[finite] >= options.range_nm[0])
        assert np.all(shown[finite] <= options.range_nm[1])


def test_fitted_error_bars_and_percentages_use_the_fit_uncertainties(plotted):
    model, config, top, bottom = plotted
    data = {obs.key: obs for obs in observations(config)}
    residuals = {line.get_label(): line.get_ydata() for line in bottom.lines}
    for container in top.containers:
        label = container.get_label()
        if not label.endswith(" (fitted)"):
            continue
        key = label.removesuffix(" (fitted)")
        obs = data[key]
        scale = model.record["scales"][key]
        reference = model.polarizability_au(obs.omega)
        percent_sigma = obs.sigma / scale / np.abs(reference) * 100
        marker, _, bars = container.lines
        segments = np.asarray(bars[0].get_segments())
        errors = (segments[:, 1, 1] - segments[:, 0, 1]) / 2
        # Artist endpoints add/subtract the error from the percentage deviation.
        np.testing.assert_allclose(errors, percent_sigma, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(
            marker.get_ydata(),
            residuals[key] * obs.sigma / scale / reference * 100,
            rtol=1e-10,
            atol=1e-12,
        )


@pytest.mark.parametrize(
    ("gas", "key", "role"),
    [
        ("Ar", "gaiser_2018/molar_polarizability_Ar", "fitted"),
        ("He", "puchalski_2016/alpha", "validator"),
        ("CH4", "reinsch_1985/mcscf_CH4", "validator"),
    ],
)
def test_electric_data_are_compared_with_electric_predictions(
    helpers, monkeypatch, gas, key, role
):
    model = fitted_models()[f"fit/{gas}"]
    config = load_config(gas)
    dataset = next(d for d in datasets(gas) if d.key == key)
    # Exercise the comparison independently of the gallery's eligibility filter:
    # helium's finite-frequency electric check also tests the retardation convention.
    monkeypatch.setattr(helpers, "plottable", lambda gas: [dataset])
    figure = Figure()
    top, bottom = figure.subplots(2)
    try:
        helpers._fit_data(top, bottom, model, config, (100e-9, 20e-6))
        container = next(
            c for c in top.containers if c.get_label() == f"{key} ({role})"
        )
        omega = angular_frequency_au(dataset.vacuum_wavelength())
        # These datasets are electric; the methane theory is also electronic only.
        # Evaluate the electronic oscillator sum independently of the fitter.
        count = len(model.strengths) - model.n_fixed
        strengths = np.asarray(model.strengths)[:count]
        energies = np.asarray(model.energies)[:count]
        prediction = np.sum(
            strengths[:, None] / (energies[:, None] ** 2 - omega**2), axis=0
        )
        values = reduction.polarizability(dataset) / AU_POLARIZABILITY
        reference = model.polarizability_au(omega)
        expected = (values - prediction) / reference * 100
        # Round-off from the equivalent SI/atomic-unit evaluation routes.
        np.testing.assert_allclose(
            container.lines[0].get_ydata(), expected, rtol=1e-7, atol=1e-12
        )
        if gas == "Ar":
            residual = next(line for line in bottom.lines if line.get_label() == key)
            # The stored fit agrees within one sigma; comparing to the refractive
            # value instead gave about nine sigma and hid the point outside the axes.
            assert abs(residual.get_ydata()[0]) < 1
    finally:
        figure.clear()

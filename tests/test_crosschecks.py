# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

from optmat.constants import C
from optmat.datasets import datasets
from optmat.gas.crosschecks import (
    Comparison,
    group_delay_dispersions,
    phase_match_ratios,
)
from optmat.gas.models.registry import default_model


def test_deviation_is_relative_to_the_measurement():
    comparison = Comparison("a/b", "r", 1e-6, 2.0, 0.1, 2.02)
    assert comparison.deviation == pytest.approx(0.01)


def test_ratio_prediction_is_the_inverse_ratio_of_the_dispersions():
    comparison = next(
        c for c in phase_match_ratios() if c.key == "shelton_1990/density_He_over_Ar"
    )
    points = np.array([comparison.wavelength, comparison.wavelength / 2])
    helium = default_model("He").polarizability(points, warn=False)
    argon = default_model("Ar").polarizability(points, warn=False)
    expected = (argon[1] - argon[0]) / (helium[1] - helium[0])
    # The same arithmetic on the same models: round-off only.
    assert comparison.predicted == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_fitted_ratios_agree_within_their_uncertainties():
    # The hydrogen isotopes and the He/H2 ratio are fitted (TC-96); each should lie
    # within 1.5 sigma (the largest pull among them is 1.0).
    fitted = {
        "shelton_1986/ratio_D2_H2",
        "shelton_1990/density_H2_over_D2",
        "shelton_1990/density_He_over_H2",
    }
    pulls = [
        (c.predicted - c.measured) / c.sigma
        for c in phase_match_ratios()
        if c.key in fitted
    ]
    assert len(pulls) == 13
    assert np.max(np.abs(pulls)) < 1.5


def test_group_delay_dispersion_is_the_second_difference_at_the_dataset_state():
    comparison = next(
        c for c in group_delay_dispersions() if c.key == "borzsonyi_2008/sgdd_He"
    )
    model = default_model("He")
    (measured,) = datasets("He", quantity="group_delay_dispersion")
    state = measured.state
    assert state is not None
    lam = comparison.wavelength
    step = lam * 1e-4
    n = model.refractivity(np.array([lam - step, lam, lam + step]), state, warn=False)
    second = (n[0] - 2 * n[1] + n[2]) / step**2
    expected = lam**3 / (2 * np.pi * C**2) * second * 1e30
    # Two step sizes of a centred difference (1e-3 and 1e-4 of lambda) differ at the
    # 1e-6 level; the prediction must also be at the dataset's state (1 bar, 19.5 C).
    assert comparison.predicted == pytest.approx(expected, rel=1e-4, abs=0.0)

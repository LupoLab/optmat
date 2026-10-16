# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The frequency derivative of every gas model's polarizability (OEP 1, 4.1 and 6.4).

Every registry model gives ``d alpha / d omega`` analytically. The tests compare it with
a four-point central difference of ``complex_polarizability`` in ``omega``, whose error
falls as the fourth power of the step, away from poles; check the measured absorption
(O3) branch by branch, including at a knot and above the grid; and check the static
limit, helium's retardation term and the temperature metadata of section 4.9.
"""

import math
import tomllib
import warnings

import numpy as np
import pytest

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_POLARIZABILITY, C
from optmat.exceptions import ModelRangeWarning
from optmat.gas import GAS_DATA_DIRECTORY
from optmat.gas.lorentz_lorenz import (
    dpolarizability_dsusceptibility,
    polarizability_from_susceptibility,
)
from optmat.gas.models.fitted import fitted_models
from optmat.gas.models.oscillator import OscillatorModel
from optmat.gas.models.registry import _all_models, model

MODELS = sorted(_all_models())


def _difference(m, omega, step):
    """Return the four-point central difference of alpha in omega (rad/s)."""

    def alpha(w):
        return m.complex_polarizability(2 * math.pi * C / w, warn=False)

    h = step * omega
    return (
        alpha(omega - 2 * h)
        - 8 * alpha(omega - h)
        + 8 * alpha(omega + h)
        - alpha(omega + 2 * h)
    ) / (12 * h)


@pytest.mark.parametrize("name", MODELS)
def test_the_derivative_is_the_difference_of_the_polarizability(name):
    m = model(name)
    # Inside the qualified intervals, away from their edges: between two of them a
    # molecule's bands are poles (OEP 3). One interval: five points across it, 5 %
    # inside its ends; several: the quarter points of each, in log wavelength.
    if len(m.intervals) == 1:
        low, high = m.wavelength_range
        wavelength = np.geomspace(low * 1.05, min(high, 20e-6) * 0.95, 5)
    else:
        wavelength = np.concatenate(
            [np.geomspace(low, high, 5)[[1, 3]] for low, high in m.intervals]
        )
    omega = 2 * math.pi * C / wavelength
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        analytic = m.dpolarizability_domega(wavelength, warn=False)
        coarse = _difference(m, omega, 1e-4)
        fine = _difference(m, omega, 5e-5)
    scale = np.maximum(np.abs(analytic), 1e-300)
    error_coarse = np.abs(coarse - analytic) / scale
    error_fine = np.abs(fine - analytic) / scale
    # Within 1e-5 at the finer step. Where the difference is not yet at its rounding
    # floor (up to about 2e-7 here), halving the step must cut its error by nearly 16:
    # the difference converges to the analytic value, not merely near it. O3's spline
    # is tested between its knots below, since these steps cross them.
    if name == "fit/O3":
        return
    assert np.all(error_fine < 1e-5), (name, error_fine)
    converging = error_coarse > 1e-6
    assert np.all(error_fine[converging] < error_coarse[converging] / 8), name


@pytest.mark.parametrize("name", MODELS)
def test_the_static_limit_is_zero(name):
    value = complex(model(name).dpolarizability_domega(math.inf, warn=False))
    # The seven refractiveindex.info entries of form 2 are NaN in the static limit,
    # as their polarizability is; every other model is exactly 0 there.
    assert value == 0 or (math.isnan(value.real) and name.startswith("refractiveindex"))


def test_the_derivative_warns_outside_the_range():
    m = model("fit/Ar")
    with pytest.warns(ModelRangeWarning):
        m.dpolarizability_domega(m.wavelength_range[0] / 2)


@pytest.fixture(scope="module")
def ozone():
    return model("fit/O3")


class TestOzone:
    """The measured absorption's derivative, branch by branch (OEP 1, section 4.1)."""

    def _at(self, ozone, nu_cm, step_cm):
        omega = 2 * math.pi * C * np.asarray(nu_cm) * 100
        analytic = ozone.dpolarizability_domega(1e-2 / np.asarray(nu_cm), warn=False)
        h = 2 * math.pi * C * 100 * step_cm

        def alpha(w):
            return ozone.complex_polarizability(2 * math.pi * C / w, warn=False)

        numeric = (
            alpha(omega - 2 * h)
            - 8 * alpha(omega - h)
            + 8 * alpha(omega + h)
            - alpha(omega + 2 * h)
        ) / (12 * h)
        return analytic, numeric

    def test_between_the_knots(self, ozone):
        # Midway between knots (the grid's step is 0.5 cm^-1) with steps well inside
        # one segment: the spline is a cubic there and the interpolated cross section
        # linear, so the difference is exact to rounding.
        analytic, numeric = self._at(ozone, [20000.25, 40000.25, 52000.75], 0.01)
        np.testing.assert_allclose(analytic.real, numeric.real, rtol=1e-7)
        np.testing.assert_allclose(analytic.imag, numeric.imag, rtol=1e-7)

    def test_at_a_knot_the_imaginary_part_takes_the_mean_slope(self, ozone):
        absorption = ozone.absorption
        per_hartree = AU_ANGULAR_FREQUENCY / (2 * math.pi * C) / 100
        nu = 40000.0  # a knot
        omega = np.array([nu / per_hartree])
        below = absorption.complex_derivative_au(omega - 0.1 / per_hartree)[1]
        above = absorption.complex_derivative_au(omega + 0.1 / per_hartree)[1]
        at = absorption.complex_derivative_au(omega)[1]
        # The slope part of the derivative is the mean; the 1/omega^2 part is smooth.
        assert at == pytest.approx(0.5 * (below + above), rel=1e-3, abs=0.0)

    def test_above_the_grid(self, ozone):
        top = ozone.absorption.max_wavenumber_cm
        analytic, numeric = self._at(ozone, [top * 1.2, top * 2], 1.0)
        np.testing.assert_allclose(analytic.real, numeric.real, rtol=1e-7)
        assert np.all(analytic.imag == 0)

    def test_the_first_order_dispersion_uses_the_splines_derivative(self, ozone):
        omega = np.array([0.1, 0.2])
        first = ozone.dispersion_au(omega, 1)
        derivative = ozone.dpolarizability_domega_au(omega)
        # d(omega alpha)/d omega = alpha + omega alpha', with the measured term's alpha
        # from the grid (as dispersion_au has it) and alpha' from the spline.
        spline = ozone.absorption.complex_derivative_au(omega)[0]
        grid = ozone.absorption.real_derivative_au(omega, 1)
        assert not np.allclose(spline, grid, rtol=0, atol=0)
        expected = ozone.dispersion_au(omega, 0) / omega + omega * derivative
        np.testing.assert_allclose(first, expected, rtol=1e-12)


def test_heliums_retardation_term_enters_the_derivative():
    helium = model("fit/He")
    assert isinstance(helium, OscillatorModel)
    omega = 0.2  # hartree
    wavelength = 2 * math.pi * C / (omega * AU_ANGULAR_FREQUENCY)
    import dataclasses  # noqa: PLC0415

    plain = dataclasses.replace(helium, retardation_au=0.0)
    difference = complex(
        helium.dpolarizability_domega(wavelength, warn=False)
        - plain.dpolarizability_domega(wavelength, warn=False)
    )
    expected = (
        2 * helium.retardation_au * omega * AU_POLARIZABILITY / AU_ANGULAR_FREQUENCY
    )
    assert difference.real == pytest.approx(expected, rel=1e-9, abs=0.0)


def test_the_jacobian_of_inverse_lorentz_lorenz():
    chi, density, virial = 6e-4, 2.6e25, 1e-6
    h = 1e-10
    numeric = (
        polarizability_from_susceptibility(chi + h, density, virial)
        - polarizability_from_susceptibility(chi - h, density, virial)
    ) / (2 * h)
    assert dpolarizability_dsusceptibility(chi, density, virial) == pytest.approx(
        numeric, rel=1e-7, abs=0.0
    )


class TestTemperatures:
    """The reference temperatures and coefficients as metadata (OEP 1, section 4.9)."""

    def test_the_curated_coefficients_are_the_fits(self):
        with (GAS_DATA_DIRECTORY / "temperatures.toml").open("rb") as handle:
            table = tomllib.load(handle)["coefficients"]
        for gas, entry in table.items():
            with (GAS_DATA_DIRECTORY / "fits" / f"{gas}.toml").open("rb") as handle:
                fit = tomllib.load(handle)
            options = [
                d["temperature"]
                for d in fit.get("fitted", []) + fit.get("validators", [])
                if isinstance(d, dict) and d.get("key") == entry["dataset"]
            ]
            assert options, gas
            assert options[0]["coefficient"] == entry["coefficient_per_k"]
            assert options[0]["sigma"] == entry["sigma_per_k"]

    def test_the_metadata_of_the_fits(self):
        fits = fitted_models()
        assert fits["fit/N2"].reference_temperature == 293.15
        assert fits["fit/N2"].infrared_temperature is None
        assert fits["fit/CO2"].infrared_temperature == 296.0
        assert fits["fit/SF6"].infrared_temperature == 298.1
        assert fits["fit/CH4"].reference_temperature is None
        assert fits["fit/Ar"].temperature_coefficient is None
        coefficient = fits["fit/CO2"].temperature_coefficient
        assert coefficient is not None
        assert coefficient.wavelength == pytest.approx(632.9908e-9)
        identity = fits["fit/CO2"].identity
        assert (identity.reference_temperature, identity.infrared_temperature) == (
            293.15,
            296.0,
        )

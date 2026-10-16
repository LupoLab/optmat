# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The complex-index interface of the models: transparent and absorbing."""

import dataclasses
import warnings

import numpy as np
import pytest

from optmat.conditions import ONE_ATM_0C
from optmat.constants import AU_POLARIZABILITY
from optmat.exceptions import ModelAbsorptionWarning
from optmat.gas.gases import Gas
from optmat.gas.models.base import Model
from optmat.gas.models.registry import default_model, refractive_index, refractivity

WAVELENGTHS = np.array([400e-9, 800e-9, 1800e-9])


@pytest.mark.parametrize("gas", [Gas.Ar, Gas.N2, Gas.O2])
def test_a_transparent_model_has_a_zero_imaginary_part(gas):
    chosen = default_model(gas)
    alpha = chosen.complex_polarizability(WAVELENGTHS)
    assert np.all(np.imag(alpha) == 0)
    np.testing.assert_array_equal(np.real(alpha), chosen.polarizability(WAVELENGTHS))


@pytest.mark.parametrize("gas", [Gas.Ar, Gas.N2, Gas.O2])
def test_the_complex_index_agrees_with_the_refractivity(gas):
    n = refractive_index(gas, WAVELENGTHS, ONE_ATM_0C)
    expected = refractivity(gas, WAVELENGTHS)
    np.testing.assert_allclose(np.real(n) - 1, expected, rtol=1e-12)
    assert np.all(np.imag(n) == 0)


@dataclasses.dataclass(frozen=True)
class _Absorber(Model):
    """A Lorentzian line at 500 nm with a small width, for the interface tests."""

    @property
    def absorbing(self):
        return True

    def _polarizability(self, wavelength):
        return np.real(self._complex_polarizability(wavelength))

    def _complex_polarizability(self, wavelength):
        omega = 1 / np.asarray(wavelength)
        omega0, gamma = 1 / 500e-9, 1 / 5e-3
        denominator = omega0**2 - omega**2 - 1j * gamma * omega
        return 5 * AU_POLARIZABILITY * omega0**2 / denominator

    def _dcomplex_polarizability_domega(self, wavelength):
        # The same line, differentiated in its variable 1/lambda, times d(1/lambda)/d
        # omega = 1 / (2 pi c).
        omega = 1 / np.asarray(wavelength)
        omega0, gamma = 1 / 500e-9, 1 / 5e-3
        denominator = omega0**2 - omega**2 - 1j * gamma * omega
        slope = 5 * AU_POLARIZABILITY * omega0**2 * (2 * omega + 1j * gamma)
        return slope / denominator**2 / (2 * np.pi * 299792458.0)


ABSORBER = _Absorber(Gas.CO, "test/absorber", "none", (100e-9, 10e-6), "theory")


def test_an_absorbing_model_warns_where_the_real_part_alone_is_used():
    with pytest.warns(ModelAbsorptionWarning, match="absorbs"):
        ABSORBER.polarizability(np.array([500e-9]))


def test_an_absorbing_model_is_silent_far_from_its_band():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ABSORBER.polarizability(np.array([5e-6]))


def test_absorption_gives_a_positive_imaginary_index():
    n = ABSORBER.refractive_index(np.array([500e-9]))
    assert np.imag(n)[0] > 0


def test_the_complex_index_accepts_a_single_wavelength():
    # A scalar wavelength reaches the measured-absorption term as a 0-d array.
    from optmat.conditions import State  # noqa: PLC0415
    from optmat.gas.models.registry import refractive_index  # noqa: PLC0415

    n = refractive_index("O3", 600e-9, State(101325.0, 293.15))
    assert np.iscomplexobj(n)
    assert np.imag(n) > 0

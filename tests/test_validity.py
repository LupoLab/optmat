# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The wavelengths a model accepts: positive ones, and inf for a gas's static limit."""

import math

import numpy as np
import pytest

from optmat import _validity
from optmat.exceptions import ArgumentError
from optmat.gas.air import air_model
from optmat.gas.models.registry import model as gas_model
from optmat.solid.models.registry import model as solid_model

SOLID_METHODS = (
    "permittivity",
    "dpermittivity_domega",
    "refractive_index",
    "dindex_domega",
    "group_index",
    "gvd",
)
GAS_METHODS = (
    "polarizability",
    "complex_polarizability",
    "dpolarizability_domega",
    "refractivity",
    "refractive_index",
)


@pytest.mark.parametrize("wavelength", [0.0, -800e-9, -math.inf])
@pytest.mark.parametrize("static", [True, False])
def test_a_wavelength_at_or_below_zero_is_refused(wavelength, static):
    with pytest.raises(ArgumentError, match="test: wavelengths must be positive"):
        _validity.checked_wavelength([800e-9, wavelength], "test", static=static)


def test_inf_is_the_static_limit_only_where_it_is_allowed():
    assert _validity.checked_wavelength(math.inf, "test", static=True) == math.inf
    with pytest.raises(ArgumentError, match="positive and finite"):
        _validity.checked_wavelength(math.inf, "test", static=False)


def test_a_nan_wavelength_passes_as_a_missing_value():
    checked = _validity.checked_wavelength([np.nan, 1e-6], "test", static=False)
    assert np.isnan(checked[0])
    assert checked[1] == 1e-6


@pytest.mark.parametrize("method", SOLID_METHODS)
@pytest.mark.parametrize("name", ["fit/SiO2", "malitson_1965/sellmeier", "table/Ag"])
@pytest.mark.parametrize("wavelength", [0.0, -1e-6, math.inf])
def test_a_solid_model_refuses_an_unphysical_wavelength(method, name, wavelength):
    # Before the check, fit/SiO2 gave n = 0.9998 at -1 um: a Sellmeier term is even
    # in lambda, so a negative wavelength returned finite nonsense, not an error.
    with pytest.raises(ArgumentError, match=name):
        getattr(solid_model(name), method)(np.array([800e-9, wavelength]), warn=False)


@pytest.mark.parametrize("method", GAS_METHODS)
@pytest.mark.parametrize("wavelength", [0.0, -1e-6])
def test_a_gas_model_refuses_an_unphysical_wavelength(method, wavelength):
    with pytest.raises(ArgumentError, match="fit/Ar"):
        getattr(gas_model("fit/Ar"), method)(wavelength, warn=False)


@pytest.mark.parametrize("method", GAS_METHODS)
def test_air_refuses_an_unphysical_wavelength(method):
    with pytest.raises(ArgumentError, match="air/dry air"):
        getattr(air_model(), method)(-1e-6, warn=False)


def test_a_gas_and_air_keep_their_static_limit():
    static = gas_model("fit/Ar").refractivity(math.inf, warn=False, ideal=True)
    assert np.isfinite(static)
    assert np.isfinite(air_model().refractivity(math.inf, warn=False, ideal=True))

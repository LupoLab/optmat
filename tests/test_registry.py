# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import subprocess
import sys

import numpy as np
import pytest

from optmat.conditions import ONE_ATM_0C
from optmat.exceptions import OptmatError
from optmat.gas.gases import Gas
from optmat.gas.models.registry import (
    _all_models,
    default_model,
    model,
    models,
    refractivity,
)


@pytest.mark.parametrize("gas", list(Gas))
def test_every_gas_has_a_recommended_model_of_its_own(gas):
    chosen = default_model(gas)
    assert chosen.gas is gas
    assert chosen.name in models(gas)


def test_unknown_model_name_raises():
    with pytest.raises(OptmatError, match="no model named"):
        model("fit/Og")


def test_a_name_is_looked_up_in_the_registry_that_holds_it():
    # model() builds only the registry its name's prefix names; it must find every
    # model the full listing holds (the same model: a test may rebuild a registry).
    for name, listed in _all_models().items():
        assert model(name).identity.digest == listed.identity.digest


def test_a_fit_is_looked_up_without_the_source_catalogue():
    # The published formulas read every source record (seconds); the recommended
    # models are fits, whose lookup does not need them (TC-810).
    code = (
        "import optmat.datasets, optmat.gas.models.registry as r; "
        "r.default_model('Ar'); r.model('fit/O3'); "
        "print(optmat.datasets.sources.cache_info().currsize)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "0"


def test_refractivity_uses_the_default_or_a_named_model():
    wavelength = np.array([500e-9, 1e-6])
    default = refractivity("Ar", wavelength)
    named = refractivity("Ar", wavelength, ONE_ATM_0C, model_name="fit/Ar")
    np.testing.assert_array_equal(default, named)
    other = refractivity(
        "Ar", wavelength, ONE_ATM_0C, model_name="borzsonyi_2008/sellmeier_Ar"
    )
    # Borzsonyi's argon is 0.11-0.15 % above the fit in the visible
    # (docs/gas_models.md).
    assert np.all((other / default - 1 > 5e-4) & (other / default - 1 < 3e-3))


def test_scalar_wavelength_gives_a_zero_dimensional_array():
    value = refractivity("He", 800e-9)
    assert isinstance(value, np.ndarray)
    assert value.shape == ()

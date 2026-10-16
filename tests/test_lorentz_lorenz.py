# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from optmat.constants import AU_POLARIZABILITY, EPS0, N_A
from optmat.gas.lorentz_lorenz import (
    molar_refractivity,
    polarizability_from_refractivity,
    polarizability_from_susceptibility,
    refractivity,
    susceptibility,
    susceptibility_dilute,
)

ALPHAS = st.floats(min_value=0.1, max_value=100.0).map(lambda a: a * AU_POLARIZABILITY)
DENSITIES = st.floats(min_value=1e20, max_value=1e28)


@given(ALPHAS, DENSITIES)
def test_susceptibility_round_trip(alpha, density):
    chi = float(susceptibility(alpha, density))
    back = float(polarizability_from_susceptibility(chi, density))
    # Two rational maps of a y < 0.2: a few ulps each.
    assert back == pytest.approx(alpha, rel=1e-13, abs=0.0)


@given(ALPHAS, DENSITIES)
def test_refractivity_round_trip(alpha, density):
    dn = float(refractivity(alpha, density))
    back = float(polarizability_from_refractivity(dn, density))
    assert back == pytest.approx(alpha, rel=1e-13, abs=0.0)


@given(ALPHAS, st.floats(min_value=1e20, max_value=1e25))
def test_dilute_limit_agrees_to_first_order(alpha, density):
    chi = float(susceptibility(alpha, density))
    dilute = float(susceptibility_dilute(alpha, density))
    # chi / dilute = 1 / (1 - y) with y = dilute / 3: the first-order correction.
    assert chi / dilute - 1 == pytest.approx(dilute / 3, rel=1e-3, abs=0.0)


def test_refractivity_has_no_cancellation_for_tiny_values():
    alpha = 1.0 * AU_POLARIZABILITY
    dn = float(refractivity(alpha, 1e10))
    assert dn == pytest.approx(1e10 * alpha / (2 * EPS0), rel=1e-12, abs=0.0)


def test_arrays_keep_shape():
    alpha = np.full((3, 2), 11.0 * AU_POLARIZABILITY)
    assert susceptibility(alpha, 2.7e25).shape == (3, 2)
    assert refractivity(alpha, 2.7e25).shape == (3, 2)


def test_molar_refractivity_definition():
    alpha = 2.0 * AU_POLARIZABILITY
    assert float(molar_refractivity(alpha)) == pytest.approx(
        N_A * alpha / (3 * EPS0), rel=1e-15, abs=0.0
    )


@given(ALPHAS, DENSITIES, st.floats(min_value=-1e-5, max_value=1e-5))
def test_round_trip_with_the_refractivity_virial(alpha, density, virial):
    dn = float(refractivity(alpha, density, virial))
    back = float(polarizability_from_refractivity(dn, density, virial))
    assert back == pytest.approx(alpha, rel=1e-12, abs=0.0)


def test_virial_term_is_the_pair_correction():
    # Argon-like: B_R/A_R = 1.7 cm^3/mol / 4.14 cm^3/mol... as m^3/mol, at 1 atm
    # (4.46e-5 mol/cm^3): y grows by b rho, 1.6e-5 here.
    alpha, density, virial = 11.08 * AU_POLARIZABILITY, 2.687e25, 0.41e-6
    plain = float(susceptibility(alpha, density))
    paired = float(susceptibility(alpha, density, virial))
    y = density * alpha / (3 * EPS0)
    expected = (
        3 * y * (1 + virial * density / N_A) / (1 - y * (1 + virial * density / N_A))
    )
    assert paired == pytest.approx(expected, rel=1e-14, abs=0.0)
    assert paired / plain - 1 == pytest.approx(
        virial * density / N_A, rel=2e-3, abs=0.0
    )

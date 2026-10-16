# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import numpy as np
import pytest

from optmat.constants import AU_ENERGY, AU_POLARIZABILITY, E_CHARGE
from optmat.datasets import load_source
from optmat.exceptions import DataError
from optmat.gas import _tunnelling as tunnelling
from optmat.gas.gases import Gas
from optmat.gas.ionisation import (
    _parameters,  # pyright: ignore[reportPrivateUsage]
    adk_coefficient_squared,
    electronic_static_polarizability,
    hartree_coefficient_squared,
)
from optmat.gas.models.fitted import fitted_models

RECORD = """citation = "demo"
kind = "theory"

[[datasets]]
id = "ip"
material = "He"
quantity = "ionization_energy"
unit = "eV"
abscissa = "none"
points = [["1", "", "24.587389"]]
[datasets.parameters]
kind = "adiabatic"

[[datasets]]
id = "c_l"
material = "He"
quantity = "tunnelling_structure_coefficient"
unit = "1"
abscissa = "label"
points = [["1", "0", "3.13"]]
[datasets.parameters]
model = "moadk"
orbital = "1s"
m = 0
normalisation = "Tong et al. 2002 eq. (2)"

[[datasets]]
id = "ion"
material = "He+"
quantity = "polarizability"
unit = "au"
abscissa_unit = "nm"
points = [["1", "inf", "0.28125"]]
"""


@pytest.fixture
def by_key(tmp_path):
    path = tmp_path / "demo_2000.toml"
    path.write_text(RECORD)
    return {d.key: d for d in load_source(path).datasets}


def test_a_gas_entry_gives_its_parameters_with_their_sources(by_key):
    entry = {
        "ionization_energy": "demo_2000/ip",
        "orbital": "1s",
        "orbital_l": 0,
        "structure": {"test": "demo_2000/c_l"},
        "default_structure": "test",
        "ion_polarizability": "demo_2000/ion",
        "reasons": {"ionization_energy": "a test"},
    }
    parameters = _parameters(Gas.He, entry, by_key)
    assert parameters.ionization_energy == pytest.approx(24.587389 * E_CHARGE)
    assert parameters.ionization_energy_kind == "adiabatic"
    assert parameters.structure["test"].coefficients == {0: 3.13}
    # He+ is hydrogenic: 9 / (2 Z^4) = 9/32 au exactly.
    assert parameters.ion_polarizability == pytest.approx(0.28125 * AU_POLARIZABILITY)
    assert parameters.polarizability_difference == pytest.approx(
        parameters.polarizability - 0.28125 * AU_POLARIZABILITY
    )
    # n* = 1 / sqrt(2 I_p) with I_p in hartree.
    assert parameters.n_star == pytest.approx(
        1 / math.sqrt(2 * 24.587389 * E_CHARGE / AU_ENERGY)
    )
    assert parameters.l_star == pytest.approx(parameters.n_star - 1)
    assert parameters.reasons == {"ionization_energy": "a test"}


def test_a_wrong_material_or_default_is_refused(by_key):
    entry = {"ionization_energy": "demo_2000/ip", "orbital": "1s"}
    with pytest.raises(DataError, match="not Ne's"):
        _parameters(Gas.Ne, entry, by_key)
    with pytest.raises(DataError, match="not a tabulation"):
        _parameters(Gas.He, {**entry, "default_structure": "none"}, by_key)
    with pytest.raises(DataError, match="unknown dataset"):
        _parameters(Gas.He, {**entry, "ionization_energy": "demo_2000/nope"}, by_key)
    with pytest.raises(DataError, match="not He\\+'s"):
        _parameters(Gas.He, {**entry, "ion_polarizability": "demo_2000/ip"}, by_key)


def test_the_electronic_static_polarizability_is_the_fitted_static_limit():
    # Helium's fit rests on Puchalski et al.'s ab initio 4He static value, 1.38376077
    # au (with the relativistic, QED and finite-mass corrections); the fit reproduces
    # it to its 1e-7 backbone, and the magnetic term is left out (it is not electric).
    alpha, sigma, source = electronic_static_polarizability("He")
    assert source == "fit/He"
    assert alpha / AU_POLARIZABILITY == pytest.approx(1.38376077, rel=2e-7)
    assert sigma < 1e-6 * alpha
    # A molecule's static limit leaves out its infrared bands: CO2's full static
    # polarizability, with its vibrational part (about 2 au of 19.6, the dielectric
    # value), is 12 % larger.
    co2, _, _ = electronic_static_polarizability("CO2")
    full = float(fitted_models()["fit/CO2"].polarizability(np.inf, warn=False))
    assert 0.08 < (full - co2) / co2 < 0.16


def test_the_hartree_coefficient_and_its_stirling_form():
    # Hydrogen: C_0 = 2 exactly (psi = exp(-r) / sqrt(pi) = 2 exp(-r) Y_00).
    assert hartree_coefficient_squared(1.0, 0.0) == pytest.approx(4.0, rel=1e-15)
    # The corrected eq. (20) is Stirling's formula for the two Gamma functions of
    # eq. (19), whose error is about 1/(12 (n* + l*)) + 1/(12 (n* - l*)): 17.6 % at
    # n* = 1, l* = 0; 1.7 % and 2.0 % at n* = 10 with l* = 0 and 4; 8.9 % at l* = 9,
    # where Gamma(n* - l*) = Gamma(1) is Stirling's worst case.
    assert adk_coefficient_squared(1.0, 0.0) / 4.0 == pytest.approx(1.176, abs=1e-3)
    for l_star in (0.0, 4.0, 9.0):
        ratio = adk_coefficient_squared(10.0, l_star) / hartree_coefficient_squared(
            10.0, l_star
        )
        assert ratio == pytest.approx(1.0, abs=0.1 if l_star == 9.0 else 0.025)


def test_the_hydrogen_static_rate_is_exact_asymptotically():
    # Landau & Lifshitz: w = (4 / F) exp(-2 / 3F) for H(1s) in a weak static field.
    field = np.array([0.02, 0.04, 0.06])
    factors = tunnelling.aligned_structure_factors({0: 2.0}, 0)
    rates = [tunnelling.static_rate_au(f, 0.5, factors) for f in field]
    np.testing.assert_allclose(rates, 4 / field * np.exp(-2 / (3 * field)), rtol=1e-13)


def test_the_orientation_average_matches_a_quadrature():
    # A p_sigma orbital (l = 1, m = 0): B(m') = C Q(1, m') d^1_{m'0}(beta), with
    # d^1_00 = cos(beta) and d^1_{+-1,0} = -+sin(beta)/sqrt(2); averaging over
    # orientations (sin(beta) d(beta) / 2) must give C^2 Q^2 / 3 for every m'.
    beta = np.linspace(0, np.pi, 20001)
    weight = np.sin(beta) / 2
    c = 1.7
    averages = {
        0: np.trapezoid(
            (c * tunnelling.q_factor(1, 0) * np.cos(beta)) ** 2 * weight, beta
        ),
        1: np.trapezoid(
            (c * tunnelling.q_factor(1, 1) * np.sin(beta) / np.sqrt(2)) ** 2 * weight,
            beta,
        ),
    }
    factors = tunnelling.isotropic_structure_factors({1: c})
    # The trapezoid rule on 20001 points is exact to about 1e-8 here.
    assert factors[0] == pytest.approx(averages[0], rel=1e-7)
    assert factors[1] == pytest.approx(averages[1], rel=1e-7)
    assert factors[-1] == pytest.approx(factors[1])
    with pytest.raises(ValueError, match="exceeds"):
        tunnelling.q_factor(1, 2)


def test_the_noble_gases_have_nists_ionisation_energies():
    from optmat.gas.ionisation import ionisation_parameters  # noqa: PLC0415

    argon = ionisation_parameters("Ar")
    assert argon.ionization_energy / E_CHARGE == pytest.approx(15.7596119, rel=1e-12)
    assert argon.ionization_energy_source == "nist_asd_2026/ie_Ar_I"
    assert argon.orbital_l == 1
    # kappa = sqrt(2 I_p) = 1.0762 and n* = 1 / kappa = 0.9292 for argon.
    assert argon.n_star == pytest.approx(0.92915, abs=1e-5)
    assert argon.polarizability_source == "fit/Ar"
    with pytest.raises(DataError, match="no ionisation parameters"):
        ionisation_parameters("O3")


def test_the_model_constants_of_zhang_and_talebpour():
    from optmat.gas.ionisation import ionisation_parameters  # noqa: PLC0415

    argon = ionisation_parameters("Ar")
    assert argon.barrier_suppression_zhang == (0.16178, -1.50441, 0.32127)
    assert argon.barrier_suppression_zhang_source == "zhang_2014/tableXI_Ar"
    assert ionisation_parameters("N2").effective_charge == 0.9
    assert ionisation_parameters("O2").effective_charge == 0.53
    assert ionisation_parameters("Ar").effective_charge is None
    assert ionisation_parameters("N2").barrier_suppression_zhang is None

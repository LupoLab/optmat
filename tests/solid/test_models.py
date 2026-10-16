# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import dataclasses
import math

import numpy as np
import pytest

from optmat.conditions import State, air_refractivity
from optmat.constants import C
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.solid.models.base import complex_index
from optmat.solid.models.formulas import FormulaModel, UdmModel, published_models

# A three-pole Sellmeier with Malitson's (1965) printed coefficients, as a test fixture
# of the form (the curated record is malitson_1965).
SELLMEIER = FormulaModel(
    material="SiO2",
    name="test/sellmeier",
    citation="",
    wavelength_range=(210e-9, 3710e-9),
    kind="published",
    form="n2_sellmeier_wavelength",
    terms={
        "B": (0.6961663, 0.4079426, 0.8974794),
        "L": (0.0684043, 0.1162414, 9.896161),
    },
)

# The same formula read as an index relative to room air at 20 degrees C and 1 atm.
IN_AIR = dataclasses.replace(SELLMEIER, air=State(101325.0, 293.15))

# A temperature-dependent Sellmeier whose polynomials are Malitson's constants plus a
# linear term, so that at T = 0 it is the fixture above.
IN_TEMPERATURE = FormulaModel(
    material="SiO2",
    name="test/temperature",
    citation="",
    wavelength_range=(400e-9, 2600e-9),
    kind="published",
    form="n2_sellmeier_temperature",
    terms={
        "S1": (0.6961663, 1e-5),
        "S2": (0.4079426, 0.0),
        "S3": (0.8974794, 0.0),
        "L1": (0.0684043, 0.0),
        "L2": (0.1162414, 0.0),
        "L3": (9.896161, 0.0),
    },
    temperature=295.0,
    temperature_range=(30.0, 300.0),
)

_LITHOSIL = published_models()["franta_2016/udm_lithosil_q2"]
assert isinstance(_LITHOSIL, UdmModel)
UDM = UdmModel(
    material="SiO2",
    name="test/udm",
    citation="",
    wavelength_range=(12.4e-9, 124e-6),
    kind="published",
    udm=_LITHOSIL.udm,
)


def test_complex_index_keeps_a_tiny_extinction_coefficient():
    # numpy's sqrt gives k = 0 here; k = eps_i / (2n) keeps it to full precision.
    eps = np.array([2.1 + 2.7e-195j, 2.25 + 1e-12j, -4.0 + 1e-10j, 1e-3 + 2.0j])
    index = complex_index(eps)
    np.testing.assert_allclose(index**2, eps, rtol=1e-14)
    assert np.imag(index[0]) == pytest.approx(
        2.7e-195 / (2 * math.sqrt(2.1)), rel=1e-14, abs=0.0
    )
    assert np.all(np.real(index) >= 0)
    assert np.all(np.imag(index) >= 0)


def test_sellmeier_reproduces_the_printed_formula():
    # n at 0.5 um from the formula written out.
    lam2 = 0.5**2
    b, lam = SELLMEIER.terms["B"], SELLMEIER.terms["L"]
    n2 = 1 + sum(bi * lam2 / (lam2 - li**2) for bi, li in zip(b, lam, strict=True))
    assert np.real(SELLMEIER.refractive_index(0.5e-6)) == pytest.approx(
        math.sqrt(n2), rel=1e-15, abs=0.0
    )


@pytest.mark.parametrize(
    "model", [SELLMEIER, IN_AIR, UDM], ids=["sellmeier", "in_air", "udm"]
)
@pytest.mark.parametrize("wavelength", [250e-9, 800e-9, 1550e-9, 3000e-9])
def test_the_omega_derivative_matches_a_finite_difference(model, wavelength):
    # Fourth-order central difference in omega with step 1e-4 omega: truncation below
    # 1e-10 relative for these smooth, transparent regions.
    omega = 2 * math.pi * C / wavelength
    h = 1e-4 * omega
    w = omega + h * np.array([-2.0, -1.0, 1.0, 2.0])
    f = model.permittivity(2 * math.pi * C / w)
    numeric = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * h)
    exact = model.dpermittivity_domega(wavelength)
    assert abs(exact - numeric) <= 1e-8 * abs(exact)


@pytest.mark.parametrize(
    "name",
    ["malitson_1965/sellmeier", "daimon_2002/sellmeier_20C", "dodge_1984/sellmeier_o"],
)
@pytest.mark.parametrize("wavelength", [600e-9, 1000e-9, 2000e-9])
def test_the_derivative_of_an_index_made_absolute_is_smooth(name, wavelength):
    # Air- and nitrogen-relative formulas (standard-air and vacuum wavelengths): the
    # derivative at 41 frequencies 1e-5 apart departs from a quartic through them by
    # 1e-14 to 9e-14 of itself (rms), the rounding of the formula itself; differencing
    # the air's squared index instead of its refractivity gave 2e-9 to 1.3e-8, a
    # jitter of alternating sign that spoiled differenced higher derivatives (TC-838).
    model = published_models()[name]
    omega = 2 * math.pi * C / wavelength * (1 + 1e-5 * np.arange(-20, 21))
    x = np.arange(-20, 21) / 20
    derivative = np.real(model.dindex_domega(2 * math.pi * C / omega))
    residual = derivative - np.polyval(np.polyfit(x, derivative, 4), x)
    assert np.std(residual) < 3e-13 * abs(np.mean(derivative))


@pytest.mark.parametrize("model", [SELLMEIER, UDM], ids=["sellmeier", "udm"])
def test_group_index_and_gvd_follow_from_n(model):
    # n_g = n - lambda dn/dlambda and beta_2 = lambda^3/(2 pi c^2) d2n/dlambda2 from
    # finite differences of n in wavelength (step 1 nm at 800 nm: 1e-7 relative).
    lam = 800e-9
    h = 1e-9
    n = np.real(model.refractive_index(lam + h * np.array([-2.0, -1.0, 0.0, 1.0, 2.0])))
    d1 = (n[0] - 8 * n[1] + 8 * n[3] - n[4]) / (12 * h)
    d2 = (-n[0] + 16 * n[1] - 30 * n[2] + 16 * n[3] - n[4]) / (12 * h * h)
    assert model.group_index(lam) == pytest.approx(n[2] - lam * d1, rel=1e-9, abs=0.0)
    assert model.gvd(lam) == pytest.approx(
        lam**3 / (2 * math.pi * C**2) * d2, rel=1e-5, abs=0.0
    )


def test_malitson_gvd_at_800_nm_is_stable():
    # A regression check of the units (s^2/m): 36.16 fs^2/mm is this implementation's
    # value from Malitson's printed coefficients, recorded 2026-09-26.
    assert SELLMEIER.gvd(800e-9) * 1e27 == pytest.approx(36.16, abs=0.01)


def test_outside_its_range_a_model_warns():
    with pytest.warns(ModelRangeWarning, match="outside"):
        SELLMEIER.refractive_index(5e-6)
    SELLMEIER.refractive_index(5e-6, warn=False)


def test_the_udm_absorbs_in_the_vacuum_ultraviolet_and_is_transparent_in_the_visible():
    index = UDM.refractive_index(np.array([100e-9, 600e-9]))
    assert np.imag(index[0]) > 0.5
    assert np.imag(index[1]) < 1e-12


def test_an_index_relative_to_air_is_multiplied_by_the_air_index():
    wavelength = np.array([300e-9, 800e-9, 1500e-9])
    air = 1 + air_refractivity(1 / wavelength, State(101325.0, 293.15))
    relative = np.real(SELLMEIER.refractive_index(wavelength))
    # Exact up to rounding: the model's epsilon is the relative one times n_air^2.
    assert np.real(IN_AIR.refractive_index(wavelength)) == pytest.approx(
        relative * air, rel=1e-15, abs=0.0
    )


def test_a_temperature_sellmeier_follows_its_polynomials():
    at_zero = dataclasses.replace(IN_TEMPERATURE, temperature=0.0)
    wavelength = np.array([500e-9, 1000e-9])
    # At T = 0 the polynomials are Malitson's constants (rounding only).
    assert at_zero.permittivity(wavelength) == pytest.approx(
        SELLMEIER.permittivity(wavelength), rel=1e-15, abs=0.0
    )
    warmer = IN_TEMPERATURE.at_temperature(300.0)
    # S1 grows by 1e-5 per K: n rises with temperature.
    assert np.all(
        np.real(warmer.refractive_index(wavelength))
        > np.real(IN_TEMPERATURE.refractive_index(wavelength))
    )
    with pytest.warns(ModelRangeWarning, match="outside"):
        IN_TEMPERATURE.at_temperature(350.0)
    with pytest.raises(DataError, match="no temperature"):
        SELLMEIER.at_temperature(300.0)

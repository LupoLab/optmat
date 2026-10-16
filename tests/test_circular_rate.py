# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The rate in a circularly polarised field (Barth & Smirnova 2013; TC-820)."""

import math

import numpy as np
import pytest

from optmat.gas import ionisation_validation as validation
from optmat.gas.ionisation_rates import (
    _log_circular_terms,  # pyright: ignore[reportPrivateUsage]
    circular_enhancement,
    circular_saddle_point,
    ionisation_rate_model,
    keldysh_parameter,
)


def test_the_saddle_point_in_the_adiabatic_limit():
    # Appendix F: zeta_0 = gamma^2/3 (1 - 28 gamma^2/45) for gamma << 1.
    gamma = 0.05
    expansion = gamma**2 / 3 * (1 - 28 * gamma**2 / 45)
    assert float(circular_saddle_point(gamma)) == pytest.approx(expansion, rel=1e-5)


def test_the_saddle_point_solves_its_equation_at_any_keldysh_parameter():
    gamma = np.array([0.3, 1.0, 3.0, 30.0, 1e4, 1e7])
    zeta = circular_saddle_point(gamma)
    u = np.sqrt((zeta**2 + gamma**2) / (1 + gamma**2))
    # Eq. (94) in a form that stays finite where u rounds to 1.
    artanh = np.log1p(u) - 0.5 * (np.log1p(-(zeta**2)) - np.log1p(gamma**2))
    assert artanh == pytest.approx(u / (1 - zeta), rel=1e-9)
    assert np.all((zeta > 0) & (zeta < 1))


def test_the_enhancement_in_the_adiabatic_limit():
    # Eqs. (106) and (109): g = 1 - gamma^2/15, h_+ + h_- = 3 + 2 gamma^2/3; p_0 adds
    # a term of order F/F_0 to both numerator and denominator.
    field, omega, energy = 0.1, 0.005, 0.5792
    gamma = float(keldysh_parameter(np.array(field), omega, energy))
    ratio = (2 * energy) ** 1.5 / field
    expected = math.exp(2 * ratio / 3 * gamma**2 / 15) * (1 + 2 * gamma**2 / 9)
    assert float(circular_enhancement(field, omega, energy, 1)) == pytest.approx(
        expected, rel=2e-3
    )
    # And it tends to 1.
    assert float(circular_enhancement(field, 1e-6, energy, 1)) == pytest.approx(
        1.0, abs=1e-6
    )


def _approximate(field, omega, energy):
    """Barth & Smirnova's approximate rates, eqs. (89) and (90), summed over p."""
    gamma = omega * math.sqrt(2 * energy) / field
    n0 = energy / omega * (1 + 1 / gamma**2)
    photons = np.arange(math.ceil(n0), math.ceil(n0) + 4000)
    zeta = 2 * n0 / photons - 1
    root = np.sqrt((zeta**2 + gamma**2) / (1 + gamma**2))
    exponent = -4 * n0 / (1 + zeta) * (np.arctanh(root) - root)
    common = (1 + 1 / gamma**2) ** 1.5 * np.exp(exponent)
    p0 = (
        3
        * energy
        / (16 * math.sqrt(2 * math.pi) * n0**2.5)
        * np.sum(
            common
            * (1 - zeta**2)
            * np.sqrt(1 - zeta)
            * ((1 + gamma**2) / (zeta**2 + gamma**2)) ** 1.25
        )
    )
    plus_minus = sum(
        3
        * energy
        / (8 * math.sqrt(2 * math.pi) * n0**1.5)
        * np.sum(
            common
            * (root - sign * zeta) ** 2
            / np.sqrt(1 - zeta)
            * ((1 + gamma**2) / (zeta**2 + gamma**2)) ** 0.75
        )
        for sign in (1, -1)
    )
    return p0 + plus_minus


def test_the_simple_rates_match_the_approximate_ones():
    # Their example, Ip = 0.5, F = 0.06, omega = 0.057 (800 nm, gamma = 0.95), and
    # Fig. 2's range: the simple and approximate rates 'coincide within graphical
    # resolution'.
    for field in (0.03, 0.045, 0.06):
        simple = (
            0.5
            * field
            / (2 * 1.0)
            * math.exp(float(_log_circular_terms(np.array(field), 0.057, 0.5, 1)))
        )
        assert simple == pytest.approx(_approximate(field, 0.057, 0.5), rel=0.08)


def test_the_circular_rate_of_the_model():
    static = ionisation_rate_model("Ar", math.inf)
    field = (2e18 / (8.8541878128e-12 * 299792458.0)) ** 0.5
    # A static field: the rate at the magnitude.
    assert float(static.circular_rate(field)) == pytest.approx(
        float(static.instantaneous_rate(field)), rel=1e-12
    )
    # A long wavelength tends to it.
    far = ionisation_rate_model("Ar", 1e-3)
    assert float(far.circular_rate(field)) == pytest.approx(
        float(static.instantaneous_rate(field)), rel=0.01
    )
    # At 800 nm and 1e14 W/cm^2, circular ionises less than linear at the same
    # intensity, whose peak field is larger by the square root of 2.
    model = ionisation_rate_model("Ar", 800e-9)
    amplitude = float(validation.amplitude(1e18))
    assert model.circular_rate(amplitude / math.sqrt(2)) < 0.1 * (
        model.cycle_averaged_rate(amplitude)
    )
    with pytest.raises(ValueError, match="s and p orbitals"):
        _log_circular_terms(np.array(0.05), 0.057, 0.5, 2)


@pytest.mark.parametrize(("gas", "linear"), [("Ar", 1.01), ("Kr", 0.96), ("Xe", 0.89)])
def test_the_circular_rate_against_nambus_measurement(gas, linear):
    # At 800 nm the circular data need the same intensity scale as the linear, to
    # 0.07 to 0.10: 1.08, 1.00, 0.98 against 1.01, 0.96, 0.89 (TC-820).
    circular = validation.phase_shift_test(
        f"nambu_2026/fig6_{gas}_CP", "default_static"
    )
    assert circular.scale == pytest.approx(linear, abs=0.1)
    # At 400 nm the default is too high in circular as in linear polarisation.
    assert (
        validation.phase_shift_test(f"nambu_2026/fig7_{gas}_CP", "default_static").scale
        < 0.75
    )


def test_the_s_shell_zero_field_and_molecules():
    # Helium's 1s: the s orbital's factor, eq. (103), also tends to 1 and exceeds it.
    assert float(circular_enhancement(0.1, 1e-6, 0.9036, 0)) == pytest.approx(
        1.0, abs=1e-6
    )
    assert float(circular_enhancement(0.1, 0.057, 0.9036, 0)) > 1.0
    model = ionisation_rate_model("He", 800e-9)
    assert float(model.circular_rate(0.0)) == 0.0
    # Molecules have no circular-polarisation rate.
    from optmat.exceptions import DataError  # noqa: PLC0415
    from optmat.gas.molecular_ionisation_rates import (  # noqa: PLC0415
        molecular_rate_model,
    )

    with pytest.raises(DataError, match="atoms only"):
        validation._envelope_yields(  # pyright: ignore[reportPrivateUsage]
            molecular_rate_model("N2", 800e-9), np.array([1e18]), 55e-15, circular=True
        )

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import numpy as np
import pytest
from scipy import special

from optmat.constants import EPS0, C
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import raman_lines, raman_parameters, stokes_gain
from optmat.gas.raman_response import (
    placzek_teller_q,
    raman_response,
    raman_response_status,
)

GASES = ("N2", "O2", "H2", "D2", "CH4", "SF6", "CO2", "N2O")


def _inside(response, factor=1.0):
    """A density inside the response's range (its geometric middle, times factor)."""
    low, high = response.density_range
    return factor * math.sqrt(max(low, 1e-3 * AMAGAT) * min(high, 1e3 * AMAGAT))


@pytest.mark.parametrize("gas", GASES)
def test_the_widths_are_raman_lines_widths_at_any_density(gas):
    response = raman_response(gas, 295.0, wavelength=800e-9)
    for factor in (0.5, 1.0, 1.7):
        density = _inside(response, factor)
        if not response.density_range[0] < density < response.density_range[1]:
            continue
        lines = raman_lines(gas, 295.0, density, wavelength=800e-9)
        # The three coefficients are solved from the same law at other densities:
        # the difference is the solve's rounding, 20 eps (the design's condition).
        assert response.fwhm(density) == pytest.approx(lines.fwhm, rel=1e-12)
        rotation = np.array(response.band) == "rotation"
        assert np.all(response.anisotropic[rotation] == lines.amplitude[rotation])
        assert np.all(response.isotropic[rotation] == 0)
        assert np.all(response.isotropic[~rotation] == lines.amplitude[~rotation])


@pytest.mark.parametrize("gas", ["N2", "H2", "dry air"])
def test_the_spectrum_is_the_transform_of_the_response(gas):
    response = raman_response(gas, 295.0, wavelength=800e-9)
    density = 1.2 * AMAGAT if gas != "H2" else 3 * AMAGAT
    # A step a tenth of the fastest period (H2's Q(1), 8 fs) and a window of 12
    # dephasing times of the slowest line: the tail left out is e^-12 = 6e-6.
    rate = response.decay_rates(density)
    dt = 0.1 * 2 * math.pi / np.max(response.frequency)
    n = int(12 / np.min(rate) / dt)
    t = np.arange(n) * dt
    h = response.response(t, density)
    for omega in (0.0, 2 * math.pi * C * 100 * 60.0, response.frequency[-1]):
        # The trapezoid rule on the oscillating kernel: its error is (dt Omega)^2/12
        # of the largest term, 3e-3 for dt Omega = 2 pi / 10; the comparison allows
        # 1e-2 of the transform's magnitude.
        integrand = h * np.exp(1j * omega * t)
        numeric = np.trapezoid(integrand, t)
        analytic = response.spectrum(omega, density)
        assert abs(numeric - analytic) < 1e-2 * abs(analytic)


def test_the_response_is_zero_before_the_impulse_and_finite():
    response = raman_response("H2", 295.0)
    h = response.response(np.array([-1e-12, 0.0, 1e-12]), 2 * AMAGAT)
    assert h[0] == 0.0
    assert h[1] == 0.0
    assert np.all(np.isfinite(h))


def test_the_tensor_components():
    response = raman_response("N2", 295.0)
    a = response.component_amplitudes("a")
    b = response.component_amplitudes("b")
    assert a + b == pytest.approx(response.amplitude, rel=1e-15)
    rotation = np.array(response.band) == "rotation"
    # A weak probe across the pump sees -1/2 of the rotational response (the
    # traceless drive), all of the isotropic vibrational one.
    assert a[rotation] / response.amplitude[rotation] == pytest.approx(-0.5, rel=1e-15)
    q = ~rotation
    ratio = a[q] / response.amplitude[q]
    # 1 for Q(0), which has no anisotropic share; N2's anisotropic share is 3 %.
    assert np.all(ratio > 0.9)
    assert np.all(ratio <= 1.0)


def _p2_matrix_element_squared(j, j_prime):
    """(1/(2J+1)) sum_M |<J'M| cos^2 - 1/3 |JM>|^2, by Gauss-Legendre quadrature."""
    x, w = np.polynomial.legendre.leggauss(64)
    total = 0.0
    for m in range(-j, j + 1):
        if abs(m) > j_prime:
            continue
        # |Y_JM|^2 integrated over phi: the normalised associated Legendre functions.
        y1 = special.sph_harm_y(j, m, np.arccos(x), 0.0).real
        y2 = special.sph_harm_y(j_prime, m, np.arccos(x), 0.0).real
        element = 2 * math.pi * np.sum(w * y1 * y2 * (x**2 - 1 / 3))
        total += element**2
    return total / (2 * j + 1)


@pytest.mark.parametrize("j", [1, 2, 3, 6, 11])
def test_the_placzek_teller_coefficients_by_quadrature(j):
    # The anisotropic share of a line is (4/45) b gamma^2 per molecule: the M average
    # of the squared element of P2's Cartesian part, cos^2 - 1/3.
    q = _p2_matrix_element_squared(j, j)
    s = _p2_matrix_element_squared(j, j + 2)
    o = _p2_matrix_element_squared(j, j - 2)
    assert q == pytest.approx(4 / 45 * float(placzek_teller_q(np.array([j]))[0]))
    b_s = 3 * (j + 1) * (j + 2) / (2 * (2 * j + 1) * (2 * j + 3))
    assert s == pytest.approx(4 / 45 * b_s)
    assert q + s + o == pytest.approx(4 / 45)


def test_the_anisotropic_q_share_matches_the_depolarization_ratio():
    response = raman_response("N2", 295.0)
    band = raman_parameters("N2").vibrational[0]
    rho = band.depolarization
    assert rho is not None
    q = np.array(response.branch) == "Q"
    ratio = np.sum(response.anisotropic[q]) / np.sum(response.isotropic[q])
    assert ratio == pytest.approx(4 * rho / (3 - 4 * rho), rel=1e-12)


def test_the_anisotropic_q_share_of_hydrogen_is_hunts():
    response = raman_response("H2", 295.0)
    band = raman_parameters("H2").vibrational[0]
    assert band.transition_anisotropy is not None
    q = np.array(response.branch) == "Q"
    j = response.j[q]
    gamma = band.transition_anisotropy(j)
    element = band.transition_polarizability(j)  # pyright: ignore[reportCallIssue]
    expected = (
        response.isotropic[q] * 4 / 45 * placzek_teller_q(j) * (gamma / element) ** 2
    )
    assert response.anisotropic[q] == pytest.approx(expected, rel=1e-13)
    # Q(0) has none; Q(1) about 2 %.
    assert response.anisotropic[q][j == 0] == 0.0
    share = response.anisotropic[q][j == 1] / response.isotropic[q][j == 1]
    assert 0.015 < float(share[0]) < 0.03


def test_the_damping_groups_share_a_width_at_every_density():
    response = raman_response("dry air", 295.0)
    densities = np.array([0.1, 0.7, 1.9]) * AMAGAT
    widths = response.fwhm(densities)
    frequencies = response.frequencies(densities) - response.frequency[:, None]
    for group in set(response.damping_group.tolist()):
        members = response.damping_group == group
        assert np.all(widths[members] == widths[members][0])
        assert np.all(frequencies[members] == frequencies[members][0])


def test_hydrogens_q_lines_are_one_group():
    response = raman_response("H2", 295.0)
    q = np.array(response.branch) == "Q"
    assert len(set(response.damping_group[q].tolist())) == 1


def test_hydrogens_q_lines_shift_with_density():
    response = raman_response("H2", 298.0)
    q = np.array(response.branch) == "Q"
    shift = (response.frequencies(10 * AMAGAT) - response.frequency) / (2 * math.pi)
    # Bischel & Dyer's Eq. (4) at 298 K: -96.1 MHz/amagat. The difference of two
    # frequencies of 1.2e14 Hz keeps 1e-16 x 1.2e14 / 1e9 = 1e-11 of the shift.
    assert shift[q] == pytest.approx(-961e6, rel=1e-9)
    assert np.all(shift[~q] == 0)


def test_the_range_warnings():
    response = raman_response("N2", 295.0)
    low, high = response.density_range
    with pytest.warns(ModelRangeWarning, match="Doppler"):
        response.decay_rates(0.5 * low)
    with pytest.warns(ModelRangeWarning, match="overlap"):
        response.decay_rates(2 * high)
    with pytest.warns(ModelRangeWarning, match="temperatures of the width"):
        raman_response("N2", 600.0)
    with pytest.raises(DataError, match="positive"):
        response.decay_rates(0.0)


def test_the_nitrogen_q_branch_overlaps_near_two_amagat():
    # Koszykowski et al. 1987: isolated lines at 1 atm, collapsed by 10 atm.
    low, high = raman_response("N2", 295.0).density_range
    assert 1 * AMAGAT < high < 3 * AMAGAT
    assert low < 0.1 * AMAGAT


def test_the_dicke_term_is_capped_at_the_doppler_width():
    response = raman_response("H2", 295.0)
    density = 0.01 * AMAGAT
    with pytest.warns(ModelRangeWarning):
        width = response.fwhm(density)
    q = np.array(response.branch) == "Q"
    expected = (
        response.width_constant
        + response.doppler_width
        + response.width_broadening * 0.01
    )
    assert width[q] == pytest.approx(expected[q], rel=1e-12)
    # Herring et al.'s Eq. (4) for Q(1) at 295 K: 7.15e-7 x 4155 cm-1 x (295/2)^1/2.
    assert response.doppler_width[q][0] == pytest.approx(
        7.15e-7
        * float(response.frequency[q][response.j[q] == 1][0])
        / (2 * math.pi)
        * math.sqrt(295 / 2),
        rel=1e-12,
    )


def test_the_mixture_carries_the_mole_fractions():
    air = raman_response("dry air", 295.0)
    fractions = dict(air.composition)
    assert sum(fractions.values()) == pytest.approx(1.0, rel=1e-12)
    assert set(fractions) == {"N2", "O2", "Ar"}
    nitrogen = raman_response("N2", 295.0)
    in_air = np.array(air.species) == "N2"
    assert air.isotropic[in_air] + air.anisotropic[in_air] == pytest.approx(
        fractions["N2"] * nitrogen.amplitude, rel=1e-14
    )
    # O2 and Ar broaden nitrogen's lines differently from nitrogen.
    assert not np.allclose(air.fwhm(AMAGAT)[in_air], nitrogen.fwhm(AMAGAT))


def test_a_mixture_by_fractions():
    mixture = raman_response({"N2": 0.8, "Ar": 0.2}, 295.0)
    assert set(mixture.species) == {"N2"}
    with pytest.raises(DataError, match="sum to one"):
        raman_response({"N2": 0.5, "Ar": 0.4}, 295.0)
    with pytest.raises(DataError, match="no broadening coefficient"):
        raman_response({"H2": 0.5, "He": 0.5}, 295.0)
    with pytest.raises(DataError, match="no Raman-active"):
        raman_response({"Ar": 1.0}, 295.0)


def test_the_gain_and_the_index_are_the_line_formulas():
    response = raman_response("H2", 295.0, wavelength=532e-9)
    density = 10 * AMAGAT
    lines = response.lines(density)
    assert response.stokes_gain(532e-9, density) == pytest.approx(
        stokes_gain(lines, 532e-9, density), rel=1e-15
    )
    nitrogen = raman_response("N2", 295.0)
    reference = raman_lines("N2", 295.0, AMAGAT)
    s = np.array(reference.branch) == "S"
    assert reference.fwhm is not None
    omega, rate = reference.frequency[s], math.pi * reference.fwhm[s]
    expected = AMAGAT * np.sum(reference.amplitude[s] * omega / (omega**2 + rate**2))
    assert nitrogen.rotational_nonlinear_index(AMAGAT) == pytest.approx(
        expected / (2 * EPS0 * C), rel=1e-12
    )


def test_the_truncation_of_the_ladder_is_converged():
    coarse = raman_response("N2", 295.0, population_floor=1e-8)
    fine = raman_response("N2", 295.0)
    t = np.linspace(0.0, 20e-12, 4001)
    difference = coarse.response(t, AMAGAT) - fine.response(t, AMAGAT)
    # The lines dropped hold below 1e-8 of the molecules each.
    assert np.max(np.abs(difference)) < 1e-6 * np.max(np.abs(fine.response(t, AMAGAT)))


def test_the_references_are_the_sources():
    response = raman_response("dry air", 295.0)
    assert {"herring_1986", "rosasco_1983", "holzer_2021", "millot_1992"} <= set(
        response.references
    )


def _solprop():
    import json  # noqa: PLC0415
    from pathlib import Path  # noqa: PLC0415

    data = Path(__file__).parent / "data"
    return (
        json.loads((data / "solprop_raman_response.json").read_text(encoding="utf-8")),
        json.loads((data / "solprop_raman.json").read_text(encoding="utf-8"))["raman"],
    )


@pytest.mark.parametrize("gas", ["N2", "O2", "H2", "D2", "CH4", "SF6"])
def test_the_convention_is_solprops(gas):
    # optmat's response on solprop's (2dc84d0) over the first 100 fs at 1 amagat, where
    # no line has dephased: the least-squares scale of each part is the ratio of the
    # strengths, so no factor of the convention differs between the two (TC-690).
    responses, parameters = _solprop()
    t = np.array(responses["delays_s"])
    early = t <= 100e-15
    response = raman_response(gas, 295.0)
    rotation = np.array(response.band) == "rotation"
    omega = response.frequencies(AMAGAT)
    rate = response.decay_rates(AMAGAT)
    for part, select in (("rotation", rotation), ("vibration", ~rotation)):
        if part not in responses["gases"][gas]:
            continue
        reference = np.array(responses["gases"][gas][part]["1.0"])[early]
        ours = np.array(
            [
                np.sum(
                    response.amplitude[select]
                    * np.sin(omega[select] * x)
                    * np.exp(-rate[select] * x)
                )
                for x in t[early]
            ]
        )
        scale = np.dot(ours, reference) / np.dot(reference, reference)
        if part == "vibration":
            v = parameters[gas]["vibrational"]
            k = (4 * math.pi * EPS0 * v["polarisability_derivative"]) ** 2 / (
                2 * v["reduced_mass"] * v["frequency"] * EPS0
            )
            # Within 100 fs the Q lines' spread (H2: 6 cm-1 between Q(0) and Q(1))
            # dephases them by at most 0.11 rad, cos 0.11 = 0.994: 1 %.
            assert scale == pytest.approx(
                np.sum(response.amplitude[select]) / k, rel=0.01
            )
        elif gas in ("N2", "O2"):
            # The same anisotropy and constants to their roundings.
            assert scale == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize(
    "gas",
    [
        "N2",
        "O2",
        "H2",
        "D2",
        "CH4",
        "SF6",
        "CO2",
        "N2O",
        "dry air",
        "dry air (N2, O2)",
        {"N2": 0.8, "Ar": 0.2},
    ],
)
def test_the_status_without_the_response(gas):
    # The accessor reads the status from the choices alone (solprop SEP 3, Q42): the
    # same status and reason as the response built at room temperature (TC-840).
    response = raman_response(gas, 295.0)
    assert raman_response_status(gas) == (response.status, response.status_reason)


def test_a_gas_without_a_raman_response_has_no_status():
    # The noble gases have none, alone or mixed; a molecule without curated choices,
    # and the full CIPM-2007 air, are refused as raman_response refuses them.
    assert raman_response_status("Ar") is None
    assert raman_response_status({"He": 0.5, "Ar": 0.5}) is None
    with pytest.raises(DataError, match="no choices for CO"):
        raman_response_status("CO")
    with pytest.raises(DataError, match="CIPM-2007"):
        raman_response_status("dry air (CIPM-2007)")

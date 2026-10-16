# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import numpy as np
import pytest

from optmat.constants import AU_POLARIZABILITY, EPS0, HBAR, K_B, C
from optmat.exceptions import DataError
from optmat.gas.raman import (
    RotationalConstants,
    rotational_lines,
    rotational_populations,
    spin_statistics,
    vibrational_lines,
)

ATOMIC_MASS = 1.66053906660e-27  # kg (CODATA 2018), for solprop's reduced mass


@pytest.mark.parametrize(
    ("gas", "even", "odd"),
    [
        ("H2", 1, 3),
        ("D2", 6, 3),
        ("N2", 6, 3),
        ("O2", 0, 1),
        ("CO2", 1, 0),
        ("N2O", 1, 1),
    ],
)
def test_the_nuclear_spin_weights(gas, even, odd):
    weights = spin_statistics(gas)
    assert (weights.even, weights.odd) == (even, odd)


def test_a_spherical_top_or_an_atom_has_no_rotational_spectrum():
    for gas in ("CH4", "SF6", "Ar"):
        with pytest.raises(DataError, match="not a linear molecule"):
            spin_statistics(gas)


def test_the_s_branch_is_solprops_form():
    # solprop: A_J = -K (J+1)(J+2)/(2J+3) [p_{J+2}/(2J+5) - p_J/(2J+1)] with
    # K = (2/15) Delta alpha^2 / (hbar eps0) and p_J the level populations.
    constants = RotationalConstants(1.99, 5.7e-6)
    spin = spin_statistics("N2")
    anisotropy = 4.6 * AU_POLARIZABILITY
    lines = rotational_lines(constants, spin, anisotropy, 294.0)
    j, p = rotational_populations(constants, spin, 294.0, 62)
    k = 2 / 15 * anisotropy**2 / (HBAR * EPS0)
    expected = {
        int(jj): -k
        * (jj + 1)
        * (jj + 2)
        / (2 * jj + 3)
        * (p[jj + 2] / (2 * jj + 5) - p[jj] / (2 * jj + 1))
        for jj in j[:-2]
    }
    for jj, amplitude in zip(lines.j, lines.amplitude, strict=True):
        assert amplitude == pytest.approx(expected[int(jj)], rel=1e-13)


def test_the_rotational_response_has_its_classical_limit():
    # At k_B T >> B the populations' differences are n_J hbar Omega_J / k_B T and
    # (J+1)(J+2)/((2J+1)(2J+3)) -> 1/4, so the response integrates to
    # sum_J A_J / Omega_J = Delta alpha^2 / (30 eps0 k_B T): three quarters of the
    # classical isothermal orientational response, the rest being the Delta J = 0
    # (M-sublevel) part that a coherent drive does not excite (TC-493). The corrections
    # are O(B / k_B T), 1 % for B = 2 cm^-1 at 295 K (k_B T = 205 cm^-1).
    anisotropy = 4.6 * AU_POLARIZABILITY
    temperature = 295.0
    lines = rotational_lines(
        RotationalConstants(2.0), spin_statistics("N2"), anisotropy, temperature
    )
    integral = np.sum(lines.amplitude / lines.frequency)
    classical = anisotropy**2 / (30 * EPS0 * K_B * temperature)
    assert integral / classical == pytest.approx(1.0, abs=0.02)


def test_the_vibrational_amplitude_is_solprops_for_a_harmonic_mode():
    # solprop's H2 vibrational coupling (its doctest: K = 1.803981e-37 m^5/(V^2 s))
    # is (4 pi eps0 dalpha/dq)^2 / (2 mu Omega eps0) with dalpha/dq = 1.3e-20 m^2,
    # mu = 1.00784 u / 2 and Omega = 2 pi 124.5669 THz; the transition polarizability
    # of the harmonic mode, <0|alpha|1>^2 = (dalpha/dq)^2 hbar / (2 mu Omega), gives the
    # same number here.
    derivative = 4 * math.pi * EPS0 * 1.3e-20
    mu = 1.00784 * ATOMIC_MASS / 2
    omega = 2 * math.pi * 124.5669e12
    transition = math.sqrt(derivative**2 * HBAR / (2 * mu * omega))
    band = omega / (2 * math.pi * C * 100)
    lines = vibrational_lines(band, transition)
    assert lines.amplitude[0] == pytest.approx(1.803981e-37, rel=1e-6)
    assert lines.frequency[0] == pytest.approx(omega, rel=1e-14)


def test_the_q_branch_shares_the_band_strength_by_population():
    lower, upper = RotationalConstants(60.0, 0.05), RotationalConstants(57.0, 0.05)
    lines = vibrational_lines(
        4160.0, 2e-41, lower=lower, upper=upper, spin=spin_statistics("H2")
    )
    # The populations sum to one, to the floor that drops the highest levels.
    assert np.sum(lines.amplitude) == pytest.approx(
        (2e-41) ** 2 / (HBAR * EPS0), rel=1e-10
    )
    # Normal hydrogen at 296 K: J = 1 (ortho) holds about 66 % of the molecules.
    share = lines.amplitude / np.sum(lines.amplitude)
    assert share[lines.j == 1][0] == pytest.approx(0.66, abs=0.03)
    with pytest.raises(DataError, match="spin weights"):
        vibrational_lines(4160.0, 2e-41, lower=lower, upper=upper)


def test_the_distortion_series_is_cut_where_it_turns_over():
    # B x - D x^2 turns over at x = B / 2D = 60 (J = 7): the ladder stops there, and a
    # temperature that populates it is refused.
    j, _ = rotational_populations(
        RotationalConstants(60.0, 0.5), spin_statistics("H2"), 296.0, 20
    )
    assert j[-1] == 7
    with pytest.raises(DataError, match="distortion series"):
        rotational_lines(
            RotationalConstants(60.0, 0.5), spin_statistics("H2"), 1e-41, 5000.0
        )
    with pytest.raises(DataError, match="raise j_max"):
        rotational_lines(
            RotationalConstants(0.4), spin_statistics("N2O"), 1e-40, 296.0, j_max=10
        )


def test_every_curated_gas_gives_lines_with_widths():
    from optmat.gas.linewidths import AMAGAT  # noqa: PLC0415
    from optmat.gas.raman import raman_lines  # noqa: PLC0415

    for gas in ("N2", "O2", "H2", "D2", "CH4", "SF6"):
        lines = raman_lines(gas, 295.0, AMAGAT)
        assert lines.fwhm is not None
        assert np.all(lines.fwhm > 0)
        assert np.all(lines.amplitude > 0)
        assert np.all(lines.frequency > 0)


def test_oxygen_has_only_odd_levels():
    from optmat.gas.linewidths import AMAGAT  # noqa: PLC0415
    from optmat.gas.raman import raman_lines  # noqa: PLC0415

    lines = raman_lines("O2", 295.0, AMAGAT)
    assert np.all(lines.j % 2 == 1)


def test_the_transition_polarizability_from_a_normal_coordinate_derivative():
    from optmat.constants import ATOMIC_MASS  # noqa: PLC0415
    from optmat.gas.raman import raman_parameters  # noqa: PLC0415

    # Murphy et al.'s N2 d alpha / dQ = 0.66 A^2 amu^(-1/2) (volume units) at Bendtsen's
    # nu_0 = 2329.9168 cm-1: <0|alpha|1> = 4 pi eps0 x 0.66e-20 m^2 / sqrt(amu)
    # x sqrt(hbar / 2 omega).
    band = raman_parameters("N2").vibrational[0]
    omega = 2 * math.pi * C * 100 * 2329.9168
    expected = 4 * math.pi * EPS0 * 0.66e-20 / math.sqrt(ATOMIC_MASS)
    expected *= math.sqrt(HBAR / (2 * omega))
    assert band.transition_polarizability == pytest.approx(expected, rel=1e-12)
    assert band.band_origin == pytest.approx(2329.9168, rel=1e-12)


@pytest.mark.parametrize("gas", ["H2", "D2"])
def test_stoicheffs_constants_reproduce_teals_lines(gas):
    from optmat.datasets import datasets  # noqa: PLC0415
    from optmat.gas.raman import raman_parameters  # noqa: PLC0415

    # Stoicheff's B0, D0 and H0 (Weber 1979) against Teal & MacWood's independent S
    # lines, stated to +-0.3 cm-1 but scattered more: their two readings of H2's S(1)
    # differ by 0.6. Every line is within 1 cm-1; without H0, H2's S(3) misses by 1.2.
    constants = raman_parameters(gas).rotational.constants  # pyright: ignore[reportOptionalMemberAccess]
    (lines,) = [d for d in datasets(gas) if d.key == f"teal_1935/table1_{gas}"]
    assert constants.h > 0
    for row, value in zip(lines.rows, lines.scaled_values(), strict=True):
        low, high = sorted(int(x) for x in row.split("-"))
        position = (
            constants.energy(np.array([high]))[0] - constants.energy(np.array([low]))[0]
        )
        assert position == pytest.approx(value, abs=1.0)


def test_the_gain_needs_widths():
    from optmat.gas.raman import stokes_gain  # noqa: PLC0415

    lines = vibrational_lines(4160.0, 1e-41)
    with pytest.raises(DataError, match="widths"):
        stokes_gain(lines, 532e-9, 1e25)


def test_hunts_components_give_the_mean_element_per_level():
    from optmat.gas.raman import AnisotropyTable, raman_parameters  # noqa: PLC0415

    # (par + 2 perp) / 3 per J; Mizrahi & Shelton quote Hunt et al.'s Q(1) element of
    # H2, 0.74027 au, and Q(2) of D2, 0.61765 au.
    for gas, j, value in (("H2", 1, 0.74027), ("D2", 2, 0.61765)):
        element = raman_parameters(gas).vibrational[0].transition_polarizability
        assert isinstance(element, AnisotropyTable)
        assert element(np.array([j]))[0] / AU_POLARIZABILITY == pytest.approx(
            value, abs=2e-5
        )


def test_the_element_from_normalised_cross_sections():
    from optmat.gas.raman import raman_parameters  # noqa: PLC0415

    # Schrötter & Klöckner's Sigma for CH4 nu1 averages 8.55 (thirteen values); with
    # rho = 0, a^2 = Sigma C / (16 pi^4), C = 5.05e-48 cm6/sr.
    band = raman_parameters("CH4").vibrational[0]
    expected = math.sqrt(8.554 * 5.05e-48 / (16 * math.pi**4)) * 1e-6  # m3
    expected *= 4 * math.pi * EPS0  # C m2/V
    assert band.transition_polarizability == pytest.approx(expected, rel=2e-3)


def test_the_single_resonance_dispersion():
    from optmat.gas.raman import raman_parameters  # noqa: PLC0415

    band = raman_parameters("H2").vibrational[0]  # static elements, nu_i = 84800
    assert band.reference_wavenumber == 0.0
    assert band.dispersion_factor(None) == 1.0
    pump = 1 / 532e-7  # cm-1
    assert band.dispersion_factor(532e-9) == pytest.approx(
        84800.0**2 / (84800.0**2 - pump**2)
    )
    with pytest.raises(DataError, match="beyond the resonance"):
        band.dispersion_factor(100e-9)
    # A band curated at a visible wavelength is unity there.
    methane = raman_parameters("CH4").vibrational[0]
    assert methane.dispersion_factor(488e-9) == pytest.approx(1.0, abs=1e-3)


def test_a_dephasing_time_gives_a_collisional_width():
    from optmat.gas.raman import raman_lines  # noqa: PLC0415

    # Chen et al.'s 23.8 ps at 2.4 atm and 295 K: FWHM 1 / (pi T2) there, and half at
    # half the density.
    density = 2.4 * 101325 / (K_B * 295.0)
    lines = raman_lines("N2O", 295.0, density)
    assert lines.fwhm is not None
    assert lines.fwhm == pytest.approx(1 / (math.pi * 23.8e-12))
    half = raman_lines("N2O", 295.0, density / 2)
    assert half.fwhm is not None
    assert half.fwhm == pytest.approx(lines.fwhm / 2)


def test_the_nitrogen_q_widths_follow_rahn_and_palmer():
    from optmat.datasets import datasets  # noqa: PLC0415
    from optmat.gas.linewidths import WAVENUMBER_TO_HZ  # noqa: PLC0415
    from optmat.gas.raman import raman_parameters  # noqa: PLC0415

    # Rosasco's 298 K widths with a power law per J fitted to Rahn & Palmer's FWHM at
    # 500-1500 K (per atm at T), including the levels above Rosasco's J = 18.
    law = raman_parameters("N2").vibrational[0].widths
    by_key = {d.key: d for d in datasets("N2")}
    for temperature in (500, 730, 1000, 1500):
        dataset = by_key[f"rahn_1986/table1_{temperature}K"]
        j = np.array([int(r.strip("Q()")) for r in dataset.rows])
        density = 101325 / (K_B * temperature)
        model = law.fwhm(j, float(temperature), density, {"self": 1.0})
        measured = dataset.scaled_values() * WAVENUMBER_TO_HZ
        ratio = model / measured - 1
        assert j.max() > 18
        assert np.sqrt(np.mean(ratio**2)) < 0.06


@pytest.mark.parametrize(
    ("wavelength_nm", "element"),
    [(911.28, 0.75327), (632.80, 0.76705), (435.96, 0.80185), (237.93, 0.98923)],
)
def test_the_dispersed_element_follows_rychlewskis(wavelength_nm, element):
    from optmat.gas.raman import AnisotropyTable, raman_parameters  # noqa: PLC0415

    # Rychlewski 1983 Table VII: <00|alpha|10> of H2, static 0.74006 au, at these
    # wavelengths (printed in angstrom). Bischel & Black's measured single-resonance law
    # follows its ab initio dispersion to 0.8 % (TC-542), and Hunt et al.'s static
    # Q(0), 0.73923, is 0.11 % below Rychlewski's.
    band = raman_parameters("H2").vibrational[0]
    assert band.dispersion_factor(wavelength_nm * 1e-9) == pytest.approx(
        element / 0.74006, rel=0.009
    )
    static = band.transition_polarizability
    assert isinstance(static, AnisotropyTable)
    assert static(np.array([0]))[0] / AU_POLARIZABILITY == pytest.approx(
        0.74006, rel=0.0012
    )

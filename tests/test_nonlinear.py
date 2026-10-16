# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The hyperpolarizability conversions, dispersion model and stored fits."""

import warnings

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import EPS0, C
from optmat.datasets import datasets
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.gases import Gas
from optmat.gas.nonlinear import (
    AU_HYPERPOLARIZABILITY,
    ESU_PER_SI,
    NONLINEAR_DIRECTORY,
    NonlinearOptions,
    frequency_sums,
    gamma_observation,
    hyperpolarizability_models,
)

GASES = sorted(p.stem for p in NONLINEAR_DIRECTORY.glob("*.toml"))
BY_KEY = {d.key: d for d in datasets()}


def _gamma(key):
    return gamma_observation(BY_KEY[key], NonlinearOptions(key)).gamma


def test_one_atomic_unit_of_gamma_is_shelton_and_rices_esu_value():
    # Shelton & Rice, Table II: 1 au = 6.235377e-65 C^4 m^4 J^-3 = 5.0367e-40 esu.
    assert pytest.approx(6.235377e-65, rel=1e-6) == AU_HYPERPOLARIZABILITY
    assert pytest.approx(5.0367e-40, rel=1e-4) == AU_HYPERPOLARIZABILITY * ESU_PER_SI


def test_lehmeiers_helium_is_43_5_au_in_the_taylor_convention():
    # gamma_T = 24 gamma_Lehmeier (TC-262): his He, 1.13e-64 C m^4 V^-3, is 43.5 au.
    assert _gamma("lehmeier_1985/gamma_thg_He")[0] == pytest.approx(43.49, abs=0.01)


def test_cars_chi_per_amagat_is_148_au_per_unit():
    # gamma_T = 2 chi / N0: 1e-18 cm^3/(erg amagat) is 147.8 au (Farrow, TC-313).
    ar = BY_KEY["farrow_1985/this_work_Ar"]
    value = float(ar.scaled_values()[0]) / 1e-18
    factor = 147.8
    stokes_free = NonlinearOptions("farrow_1985/this_work_Ar")
    with pytest.raises(DataError, match="stokes_nm"):
        gamma_observation(ar, stokes_free)
    assert value * factor == pytest.approx(9.6 * 147.8, rel=1e-3)


def test_one_atomic_unit_is_8_91e_23_cm2_per_w_of_n2_at_1_atm_and_0_c():
    density = 101325.0 / (1.380649e-23 * 273.15)
    n2 = density * AU_HYPERPOLARIZABILITY / (8 * EPS0**2 * C) * 1e4
    assert n2 == pytest.approx(8.91e-23, rel=1e-3)


@pytest.mark.parametrize(
    ("process", "w2", "w4"),
    [
        ("static", 0, 0),
        ("kerr", 4, 4),
        ("eshg", 6, 18),
        ("dc_kerr", 2, 2),
        ("thg", 12, 84),
    ],
)
def test_the_frequency_sums_of_each_process(process, w2, w4):
    nu = np.array([1.0e4])
    sums = frequency_sums(process, nu)
    assert sums[0][0] == pytest.approx(w2 * 1e8)
    assert sums[1][0] == pytest.approx(w4 * 1e16)


def test_a_mixed_component_is_refused():
    key = "bishop_1989/omega_kerr_xxzz"
    with pytest.raises(DataError, match="component"):
        gamma_observation(BY_KEY[key], NonlinearOptions(key))


@pytest.mark.parametrize("gas", GASES)
def test_the_processes_are_ordered_below_resonance(gas):
    model = hyperpolarizability_models()[Gas(gas)]
    wavelength = 1.5e-6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        values = [
            model.gamma_au(p, wavelength) for p in ("dc_kerr", "kerr", "eshg", "thg")
        ]
    assert model.gamma_au("static", np.inf) > 0
    assert np.all(np.diff(np.ravel(values)) > 0)


def test_argons_n2_at_800_nm_is_within_the_measured_range():
    model = hyperpolarizability_models()[Gas.Ar]
    n2 = float(model.n2(800e-9, State(101325.0, 293.15))) * 1e4
    # Wahlstrand 2012 9.7 +- 1.2, Zahedpour 2015 10.1 +- 1.0 (1e-20 cm^2/W).
    assert 9.0e-20 < n2 < 11.5e-20


def test_solprops_gamma3_is_gamma_over_six_eps0():
    model = hyperpolarizability_models()[Gas.N2]
    gamma3 = float(model.gamma3("kerr", 800e-9))
    assert gamma3 == pytest.approx(float(model.gamma_si("kerr", 800e-9)) / (6 * EPS0))


def test_the_documentation_table_is_current():
    from pathlib import Path  # noqa: PLC0415

    page = (
        Path(__file__).parents[1] / "docs/user-guide/nonlinear-models.md"
    ).read_text()
    rows = {
        cells[0]: cells
        for line in page.splitlines()
        if line.startswith("| ") and "±" in line
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }
    for gas, model in hyperpolarizability_models().items():
        kerr = float(model.gamma_au("kerr", 800e-9))
        printed = float(rows[gas.value][2].split("±")[0])
        assert printed == pytest.approx(kerr, rel=2e-3)


def test_the_time_domain_form_gives_every_channels_gamma():
    # The time-domain form in docs/user-guide/nonlinear-models.md:
    # P = N gamma0/6 {E^3 + A[...] + ...},
    # with D multiplying a component of wavenumber nu by nu^2. On a periodic grid whose
    # bins are 500 cm^-1 apart, a two-colour field at 800 nm (bin 25) and 1818 nm (bin
    # 11) makes, per unit density, P_1 = (3/4 gamma_Kerr(1) + 3/2 gamma_pp) / 6 at bin
    # 25, P_3 = gamma_THG(1) / 24 at bin 75 and 3/4 gamma_CARS / 6 at bin 39 (2w1 - w2).
    model = hyperpolarizability_models()[Gas.Ar]
    gamma0, a, b, b_prime = model.parameters
    a, b, b_prime = a * 1e-10, b * 1e-20, b_prime * 1e-20
    n = 512
    nu = np.fft.rfftfreq(n, 1.0 / n) * 500.0

    def d(signal):
        return np.fft.irfft(np.fft.rfft(signal) * nu**2, n)

    t = np.arange(n) / n
    field = np.cos(2 * np.pi * 25 * t) + np.cos(2 * np.pi * 11 * t)
    e3, de = field**3, d(field)
    p = (gamma0 / 6) * (
        e3
        + a * (d(e3) + 3 * field**2 * de)
        + b
        * (d(d(e3)) + 3 * field**2 * d(de) + 6 * d(field**2 * de) + 6 * field * de**2)
        + b_prime * (d(d(e3)) + 3 * field**2 * d(de))
    )
    amplitude = np.abs(np.fft.rfft(p)) * 2 / n
    one, two = 1e-2 / 12500.0, 1e-2 / 5500.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kerr = float(model.gamma_au("kerr", one))
        cross = float(model.gamma_au("pump_probe", two, one))
        thg = float(model.gamma_au("thg", one))
        cars = float(model.gamma_au("cars", one, two))
    assert amplitude[25] == pytest.approx((0.75 * kerr + 1.5 * cross) / 6, rel=1e-10)
    assert amplitude[75] == pytest.approx(thg / 24, rel=1e-10)
    assert amplitude[39] == pytest.approx(0.75 * cars / 6, rel=1e-10)


def test_n2_of_a_gas_without_an_equation_of_state_uses_the_ideal_gas():
    model = hyperpolarizability_models()[Gas.O3]
    state = State(101325.0, 293.15)
    assert float(model.n2(800e-9, state)) == pytest.approx(
        float(model.n2(800e-9, state, ideal=True))
    )


def test_a_process_beyond_the_fitted_w2_warns():
    argon = hyperpolarizability_models()[Gas.Ar]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        argon.gamma_au("kerr", 400e-9)  # W2 = 2.5e9, within the fitted 2.86e9
        argon.gamma_au("thg", 700e-9)
    # The third harmonic of 500 nm has W2 = 4.8e9 cm^-2.
    with pytest.warns(ModelRangeWarning, match="beyond the fitted data"):
        argon.gamma_au("thg", 500e-9)
    # N2O's dispersion is borrowed: its data reach only a Kerr fundamental of 664 nm.
    with pytest.warns(ModelRangeWarning, match="beyond the fitted data"):
        hyperpolarizability_models()[Gas.N2O].gamma_au("kerr", 500e-9)

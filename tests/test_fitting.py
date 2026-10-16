# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import math

import numpy as np
import pytest

from optmat.exceptions import DataError
from optmat.gas.fitting import (
    FIT_DIRECTORY,
    DatasetOptions,
    FitConfig,
    Line,
    Observation,
    fit,
    load_config,
)
from optmat.gas.gases import Gas
from optmat.gas.models.oscillator import OscillatorModel

TRUE_LINE = Line("resonance", 58.4334, f_prior=0.276, f_sigma=0.01)
TRUE_EFFECTIVE = [(0.9, 1.3), (0.8, 3.0)]  # (strength, energy in hartree)


def _true_alpha(omega):
    alpha = TRUE_LINE.f_prior / (TRUE_LINE.energy**2 - omega**2)
    for f, e in TRUE_EFFECTIVE:
        alpha = alpha + f / (e**2 - omega**2)
    return alpha


def _config(**changes):
    base = {
        "name": "synthetic",
        "gas": Gas.He,
        "description": "synthetic data from a known oscillator model",
        "lines": (TRUE_LINE,),
        "effective": 2,
        "effective_min_energy": 0.9,
        "sum_priors": {},
        "chi0": 0.0,
        "fitted": (),
        "validators": (),
        "range_nm": (60.0, math.inf),
        "digest": "0" * 64,
    }
    return FitConfig(**{**base, **changes})


def _observation(key, omega, relative_sigma, rng, scale=1.0, **options):
    alpha = _true_alpha(omega)
    sigma = relative_sigma * alpha
    noisy = scale * alpha + sigma * rng.standard_normal(omega.shape)
    return Observation(
        key=key,
        omega=omega,
        alpha=noisy,
        sigma=sigma,
        systematic=None,
        options=DatasetOptions(key=key, **options),
    )


def test_fit_recovers_a_known_model_within_its_uncertainty():
    rng = np.random.default_rng(1)
    omega = np.linspace(0.0, 0.72, 40)
    data = [_observation("theory", omega, 1e-5, rng)]
    result = fit(_config(), data, starts=3)
    check = np.linspace(0.0, 0.75, 30)
    error = result.model.polarizability_au(check) - _true_alpha(check)
    sigma = result.polarizability_sigma_au(check)
    # A correct fit with correct errors lands within a few sigma everywhere.
    assert np.all(np.abs(error) < 5 * sigma)
    assert result.reduced_chi2 == pytest.approx(1.0, abs=0.5)


def test_free_scale_absorbs_a_dataset_offset():
    rng = np.random.default_rng(2)
    omega = np.linspace(0.0, 0.7, 30)
    data = [
        _observation("theory", omega, 1e-5, rng),
        _observation("offset", omega[::3], 1e-5, rng, scale=1.0025, free_scale=True),
    ]
    result = fit(_config(), data, starts=3)
    # The model scale is fixed by the unscaled dataset, so the offset dataset's
    # fitted scale is its true 1/1.0025 ... expressed as model x scale = data.
    assert result.scales["offset"] == pytest.approx(1.0025, abs=3e-5)
    assert result.scales["theory"] == 1.0


def test_line_prior_is_respected_when_data_do_not_constrain_it():
    rng = np.random.default_rng(3)
    omega = np.linspace(0.0, 0.2, 10)  # far from the line: little information on f
    data = [_observation("static", omega, 1e-3, rng)]
    result = fit(_config(effective=1), data, starts=3)
    assert result.model.strengths[0] == pytest.approx(TRUE_LINE.f_prior, abs=0.03)


def test_covariance_is_positive_definite():
    rng = np.random.default_rng(4)
    omega = np.linspace(0.0, 0.72, 40)
    result = fit(_config(), [_observation("theory", omega, 1e-5, rng)], starts=3)
    assert np.all(np.linalg.eigvalsh(result.covariance) > 0)


def test_moments_and_chi0_enter_as_electric_quantities():
    rng = np.random.default_rng(5)
    omega = np.linspace(0.0, 0.7, 25)
    chi0 = -2e-4
    refractive = _observation("refractive", omega, 1e-6, rng)
    refractive = Observation(
        key=refractive.key,
        omega=refractive.omega,
        alpha=refractive.alpha + chi0,
        sigma=refractive.sigma,
        systematic=None,
        options=refractive.options,
    )
    s2 = TRUE_LINE.f_prior / TRUE_LINE.energy**2 + sum(
        f / e**2 for f, e in TRUE_EFFECTIVE
    )
    moments = Observation(
        key="moments",
        omega=np.zeros(1),
        alpha=np.array([s2]),
        sigma=np.array([1e-6]),
        systematic=None,
        options=DatasetOptions(key="moments"),
        orders=np.array([-2.0]),
    )
    result = fit(_config(chi0=chi0), [refractive, moments], starts=3)
    # The refractive data carry chi0 and the moment does not: both are met only if
    # the model separates them, so each lands at its noise level.
    assert result.dataset_chi2["moments"][1] < 9
    assert result.model.cauchy_moment(-2) == pytest.approx(s2, abs=5e-6)
    assert float(result.model.polarizability_au(np.array(0.0))) == pytest.approx(
        s2 + chi0, abs=5e-6
    )


def test_ratios_in_either_orientation_are_met_by_the_true_model():
    rng = np.random.default_rng(6)
    theory = _observation("theory", np.linspace(0.0, 0.7, 25), 1e-5, rng)
    omega = np.array([0.0345, 0.0428])  # 1319 and 1064 nm
    delta = _true_alpha(2 * omega) - _true_alpha(omega)
    partner = 1.7 * delta  # any fixed partner dispersion
    ratios = []
    for key, numerator in (("partner_over_gas", False), ("gas_over_partner", True)):
        exact = delta / partner if numerator else partner / delta
        sigma = 1e-5 * exact
        ratios.append(
            Observation(
                key=key,
                omega=omega,
                alpha=exact + sigma * rng.standard_normal(2),
                sigma=sigma,
                systematic=None,
                options=DatasetOptions(key=key),
                harmonic=2,
                partner_delta=partner,
                model_in_numerator=numerator,
            )
        )
    result = fit(_config(), [theory, *ratios], starts=3)
    # Each ratio is met at its noise level only if its orientation is right: the
    # wrong one misses by (1.7 - 1/1.7), 10^5 sigma.
    assert result.dataset_chi2["partner_over_gas"][1] < 12
    assert result.dataset_chi2["gas_over_partner"][1] < 12


def test_fit_configuration_rejects_unknown_model_keys(tmp_path):
    text = (FIT_DIRECTORY / "He.toml").read_text()
    path = tmp_path / "bad.toml"
    path.write_text(text.replace("[model]\n", "[model]\nstart = 3\n", 1))
    with pytest.raises(DataError, match=r"unknown \[model\] keys \['start'\]"):
        load_config(path)


def test_fit_configuration_reads_its_number_of_starts():
    assert load_config("D2").starts == 20
    assert load_config("He").starts == 6


def test_warm_start_returns_to_the_same_minimum():
    rng = np.random.default_rng(7)
    data = [_observation("theory", np.linspace(0.0, 0.72, 40), 1e-5, rng)]
    first = fit(_config(), data, starts=3)
    again = fit(_config(), data, starts=0, start=first)
    # Started at the solution, the solver stops at once: the same chi-squared to its
    # convergence tolerance.
    assert again.chi2 == pytest.approx(first.chi2, rel=1e-8)


def test_fit_without_a_starting_point_or_with_another_layout_raises():
    rng = np.random.default_rng(8)
    data = [_observation("theory", np.linspace(0.0, 0.72, 20), 1e-5, rng)]
    with pytest.raises(ValueError, match="at least one"):
        fit(_config(), data, starts=0)
    three = fit(_config(effective=3), data, starts=1)
    with pytest.raises(ValueError, match="layout"):
        fit(_config(), data, starts=0, start=three)


def test_fixed_infrared_terms_enter_measured_values_but_not_electronic_ones():
    rng = np.random.default_rng(9)
    band = ("infrared test", 2e-4, 0.0043)  # f, E (about 945 cm^-1)
    omega = np.linspace(0.02, 0.6, 30)
    electronic = _observation("electronic", omega, 1e-6, rng, infrared=False)
    measured = _observation("measured", omega, 1e-6, rng)
    ir = band[1] / (band[2] ** 2 - omega**2)
    measured = Observation(
        key=measured.key,
        omega=omega,
        alpha=measured.alpha + ir,
        sigma=measured.sigma,
        systematic=None,
        options=measured.options,
    )
    result = fit(_config(infrared=(band,)), [electronic, measured], starts=3)
    # Both datasets are met at their noise level only if the band is added to the
    # measured one and not to the electronic one.
    assert result.dataset_chi2["electronic"][1] < 60
    assert result.dataset_chi2["measured"][1] < 60
    assert result.model.labels[-1] == "infrared test"
    assert result.model.strengths[-1] == band[1]


def _rayleigh_dataset(key):
    from optmat.datasets import datasets  # noqa: PLC0415

    (dataset,) = [
        d for d in datasets(quantity="rayleigh_cross_section") if d.key == key
    ]
    return dataset


def test_rayleigh_cross_section_converts_to_the_polarizability():
    from optmat.constants import BOHR_RADIUS  # noqa: PLC0415
    from optmat.gas.fitting import _rayleigh  # noqa: PLC0415

    dataset = _rayleigh_dataset("sneep_2005/rayleigh_N2")
    obs = _rayleigh(dataset, DatasetOptions(key=dataset.key))
    wavelength = dataset.vacuum_wavelength()[0]
    cross_section = dataset.scaled_values()[0] * 1e-4
    king = dataset.parameters["king_factor"]
    assert isinstance(king, float)
    k = 2 * math.pi / wavelength
    expected = (
        math.sqrt(cross_section / (8 * math.pi / 3 * k**4 * king)) / BOHR_RADIUS**3
    )
    assert obs.alpha[0] == pytest.approx(expected, rel=1e-12)
    # Half the cross section's relative uncertainty; electric (no magnetic term).
    relative = dataset.standard_uncertainty()[0] * 1e-4 / cross_section
    assert obs.sigma[0] / obs.alpha[0] == pytest.approx(relative / 2, rel=1e-12)
    assert obs.options.electric


def test_rayleigh_without_king_factor_or_beyond_resonance_is_refused():
    from optmat.gas.fitting import _rayleigh  # noqa: PLC0415

    for key, message in (
        ("chopra_1974/rayleigh_N2", "king_factor"),
        ("chopra_1974/rayleigh_Kr", "beyond a resonance"),
    ):
        dataset = _rayleigh_dataset(key)
        with pytest.raises(DataError, match=message):
            _rayleigh(dataset, DatasetOptions(key=key))


def test_argon_model_reproduces_sneeps_index_based_cross_section():
    from optmat.constants import BOHR_RADIUS  # noqa: PLC0415
    from optmat.gas.models.fitted import fitted_models  # noqa: PLC0415
    from optmat.gas.models.oscillator import angular_frequency_au  # noqa: PLC0415

    # Sneep & Ubachs's cross section computed from argon's refractive index, against
    # the same from fit/Ar: the conversion's constants and the King factor (1 for an
    # atom) are then the only difference; they agree to 0.2 %.
    dataset = _rayleigh_dataset("sneep_2005/n_based_Ar")
    wavelength = dataset.vacuum_wavelength()
    # Both come from the refractive index, magnetic term included.
    alpha = fitted_models()["fit/Ar"].polarizability_au(
        angular_frequency_au(wavelength)
    )
    k = 2 * math.pi / wavelength[0]
    cross_section = 8 * math.pi / 3 * k**4 * (alpha[0] * BOHR_RADIUS**3) ** 2 * 1e4
    assert dataset.scaled_values()[0] == pytest.approx(cross_section, rel=5e-3)


def test_the_retardation_term_enters_refractive_data_but_not_electric_ones():
    rng = np.random.default_rng(7)
    omega = np.linspace(0.0, 0.7, 25)
    # Helium's size (TC-705): -24 ppm of S(-4), a hundred times the data's noise below.
    retardation = -3.7e-5
    refractive = _observation("refractive", omega, 1e-7, rng)
    refractive = Observation(
        key=refractive.key,
        omega=refractive.omega,
        alpha=refractive.alpha + retardation * omega**2,
        sigma=refractive.sigma,
        systematic=None,
        options=refractive.options,
    )
    electric = _observation("electric", omega, 1e-7, rng, electric=True)
    result = fit(_config(retardation=retardation), [refractive, electric], starts=3)
    # Both datasets are met only if the model adds the term to the refractive one alone.
    assert result.dataset_chi2["refractive"][1] < 2 * len(omega)
    assert result.dataset_chi2["electric"][1] < 2 * len(omega)
    model = result.model
    assert model.retardation_au == retardation
    # The model's polarizability is the refractive one; its moments are electric.
    difference = model.polarizability_au(omega) - _true_alpha(omega)
    assert np.max(np.abs(difference - retardation * omega**2)) < 1e-6


def test_the_retardation_term_enters_the_model_and_its_derivatives():
    base = {
        "gas": Gas.He,
        "name": "demo",
        "citation": "",
        "kind": "fitted",
        "wavelength_range": (60e-9, math.inf),
        "strengths": (0.5,),
        "energies": (0.8,),
        "labels": ("a",),
    }
    plain = OscillatorModel(**base)
    retarded = OscillatorModel(**base, retardation_au=-3e-5)
    omega = np.array([0.1, 0.3])
    added = retarded.polarizability_au(omega) - plain.polarizability_au(omega)
    # Differences of values near 1: round-off of a few 1e-16, so absolute tolerances.
    np.testing.assert_allclose(added, -3e-5 * omega**2, rtol=0, atol=1e-15)
    # d^k (c omega^3) / d omega^k, analytically, for the orders the pages use.
    expected = {
        0: omega**3,
        1: 3 * omega**2,
        2: 6 * omega,
        3: 6 + 0 * omega,
        4: 0 * omega,
    }
    for order, term in expected.items():
        np.testing.assert_allclose(
            retarded.dispersion_au(omega, order) - plain.dispersion_au(omega, order),
            -3e-5 * term,
            atol=1e-15,
        )
    derivative = retarded.dpolarizability_domega_au(omega)
    np.testing.assert_allclose(
        derivative - plain.dpolarizability_domega_au(omega),
        -6e-5 * omega,
        rtol=0,
        atol=1e-15,
    )
    assert retarded.cauchy_moment(-4) == plain.cauchy_moment(-4)


def test_a_temperature_option_brings_values_to_the_reference_temperature():
    from optmat.datasets import datasets  # noqa: PLC0415
    from optmat.gas.fitting import _observation  # noqa: PLC0415

    dataset = next(d for d in datasets() if d.key == "egan_2026/A_R_CO2")
    plain = _observation(dataset, DatasetOptions(key=dataset.key))
    # CO2 at 303 K brought to 293.15 K with 4.22e-6 per K (egan_2026's A_theta).
    options = DatasetOptions(key=dataset.key, temperature=(303.0, 4.22e-6, 0.34e-6))
    moved = _observation(dataset, options, reference_temperature=293.15)
    gap = 293.15 - 303.0
    np.testing.assert_allclose(moved.alpha / plain.alpha, 1 + 4.22e-6 * gap, rtol=1e-12)
    # The coefficient's uncertainty over the gap adds in quadrature (3.3e-6 relative).
    added = np.sqrt(moved.sigma**2 - (plain.sigma * (1 + 4.22e-6 * gap)) ** 2)
    np.testing.assert_allclose(
        added / np.abs(moved.alpha), 0.34e-6 * abs(gap), rtol=1e-6
    )
    with pytest.raises(DataError, match="reference_temperature_k"):
        _observation(dataset, options)


def test_fit_configuration_reads_a_temperature_option(tmp_path):
    path = tmp_path / "demo.toml"
    path.write_text(
        'gas = "CO2"\nrange_nm = [250.0, 1400.0]\n'
        "[model]\neffective = 1\neffective_min_energy_hartree = 0.4\n"
        "reference_temperature_k = 293.15\n"
        '[[fitted]]\nkey = "egan_2026/A_R_CO2"\n'
        "temperature = { measured_k = 303.0, coefficient = 4.22e-6, sigma = 0.34e-6 }\n"
    )
    config = load_config(path)
    assert config.reference_temperature == 293.15
    assert config.fitted[0].temperature == (303.0, 4.22e-6, 0.34e-6)
    path.write_text(path.read_text().replace(", sigma = 0.34e-6", ""))
    with pytest.raises(DataError, match="temperature needs"):
        load_config(path)


def test_a_correction_is_added_with_or_without_a_temperature_option():
    from optmat.datasets import datasets  # noqa: PLC0415
    from optmat.gas.fitting import _observation  # noqa: PLC0415

    dataset = next(d for d in datasets() if d.key == "egan_2026/A_R_CO2")
    plain = _observation(dataset, DatasetOptions(key=dataset.key))
    corrected = _observation(
        dataset, DatasetOptions(key=dataset.key, correction_au=(0.01, 0.002))
    )
    np.testing.assert_allclose(corrected.alpha - plain.alpha, 0.01, rtol=1e-12)
    both = _observation(
        dataset,
        DatasetOptions(
            key=dataset.key,
            correction_au=(0.01, 0.002),
            temperature=(303.0, 4.22e-6, 0.34e-6),
        ),
        reference_temperature=293.15,
    )
    # The temperature factor scales the measured values; the correction is added after.
    expected = plain.alpha * (1 + 4.22e-6 * (293.15 - 303.0)) + 0.01
    np.testing.assert_allclose(both.alpha, expected, rtol=1e-12)

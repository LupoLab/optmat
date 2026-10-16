# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The ionisation rates: their limits, their pieces against quadrature, the models."""

import dataclasses
import math

import numpy as np
import pytest
from scipy import integrate

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD, E_CHARGE
from optmat.datasets import datasets
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas.ionisation import hartree_coefficient_squared, ionisation_parameters
from optmat.gas.ionisation_rates import (
    channel_factor,
    channel_function,
    coulomb_factor,
    cycle_average_factor,
    ionisation_rate_model,
    ionisation_rate_models,
    keldysh_parameter,
    ppt_exponent_function,
    short_range_rate_au,
    tong_lin_factor,
    tunnelling_rate_au,
    zhang_factor,
)
from optmat.gas.molecular_ionisation_rates import molecular_rate_model

ATOMS = ("He", "Ne", "Ar", "Kr", "Xe")
OMEGA_800 = 2 * math.pi * 299792458 / 800e-9 / AU_ANGULAR_FREQUENCY


def _dataset(key: str):
    (dataset,) = [d for d in datasets() if d.key == key]
    return dataset


# The functions of the Keldysh parameter.


def test_the_exponent_function_is_continuous_across_its_series():
    # Series and closed form are both accurate to about 1e-12 at the switch (the
    # series' first omitted term is of relative order gamma^6 = 1e-12; the closed
    # form loses eps / gamma^2 = 2e-12), so they agree to 1e-10 on either side of it.
    below, above = ppt_exponent_function(np.array([0.00999999, 0.01000001]))
    assert below == pytest.approx(above, rel=1e-10, abs=0.0)


def test_the_exponent_function_against_its_integral():
    # Perelomov et al. (1966), eqs. (32)-(33): g = (3 / 2 gamma) int_0^1 ... is the
    # tunnelling action; here checked through the identity (2 F0 / 3F) g = 2 nu
    # [arsinh gamma - gamma sqrt(1 + gamma^2) / (1 + 2 gamma^2)] that Bergé et al.'s
    # (B14) uses, to rounding.
    gamma = np.array([0.1, 0.5, 1.0, 3.0, 10.0])
    bracket = np.arcsinh(gamma) - gamma * np.sqrt(1 + gamma**2) / (1 + 2 * gamma**2)
    expected = 3 * (1 + 2 * gamma**2) * bracket / (4 * gamma**3)
    assert ppt_exponent_function(gamma) == pytest.approx(expected, rel=1e-12, abs=0.0)


@pytest.mark.parametrize("m", [0, 1, 2])
def test_the_channel_function_is_its_integral(m):
    # Eq. (56) by quadrature; quad's own tolerance is 1.5e-8.
    for x in (0.3, 1.0, 2.5, 6.0):
        value, _ = integrate.quad(
            lambda y, x=x: np.exp(y**2 - x**2) * (x**2 - y**2) ** m, 0, x
        )
        assert channel_function(m, np.array(x)) == pytest.approx(
            value, rel=1e-7, abs=0.0
        )


def test_the_channel_function_does_not_overflow():
    # e^{y^2} overflows a float at y = 27; the Kummer form does not contain it. The
    # asymptote is |m|! / (2x), reached to 1 / x^2.
    assert channel_function(1, np.array(100.0)) == pytest.approx(
        1 / 200, rel=1e-3, abs=0.0
    )


@pytest.mark.parametrize("m", [0, 1])
@pytest.mark.parametrize("gamma", [0.3, 1.0, 3.0])
def test_the_integral_over_channels_is_its_closed_form(m, gamma):
    alpha = 2 * (math.asinh(gamma) - gamma / math.sqrt(1 + gamma**2))
    beta = 2 * gamma / math.sqrt(1 + gamma**2)
    value, _ = integrate.quad(
        lambda x: (
            math.exp(-alpha * x)
            * float(channel_function(m, np.asarray(math.sqrt(beta * x))))
        ),
        0,
        np.inf,
        limit=400,
    )
    expected = (
        4
        / math.sqrt(3 * math.pi)
        / math.factorial(m)
        * gamma**2
        / (1 + gamma**2)
        * value
    )
    # quad's tolerance on a slowly decaying integrand.
    assert channel_factor(m, np.array(gamma)) == pytest.approx(
        expected, rel=1e-6, abs=0.0
    )


def test_the_channel_factor_tends_to_one():
    # Perelomov et al.: "For gamma << 1, A_m -> 1"; the correction is of order
    # gamma^2.
    assert channel_factor(0, np.array(1e-2)) == pytest.approx(1.0, abs=1e-3)
    assert channel_factor(1, np.array(1e-2)) == pytest.approx(1.0, abs=1e-3)


def test_the_sum_over_channels_tends_to_its_integral():
    # At gamma = 0.1 the summand varies by alpha = 6.6e-4 per term, so the sum is a
    # Riemann sum of the integral to that order, whatever the threshold's fraction.
    gamma = np.array([0.1])
    for threshold in (20.0, 20.3, 20.9):
        summed = channel_factor(0, gamma, np.array([threshold]))
        assert summed == pytest.approx(channel_factor(0, gamma), rel=0.05, abs=0.0)


def test_the_sum_over_channels_steps_at_a_channel_closing():
    # Just below an integer threshold the lowest channel is about to close and
    # carries a vanishing weight w_0(0) = 0; just above, the next channel is a whole
    # photon away. The two sums differ, which is the kink the integral removes.
    gamma = np.array([1.0])
    below = channel_factor(0, gamma, np.array([9.999]))
    above = channel_factor(0, gamma, np.array([10.001]))
    assert abs(float(below[0] / above[0]) - 1) > 0.03


# The rates of one electron.


def test_hydrogen_static_rate():
    # The exact leading order of hydrogen's ground state, (4 / F) exp(-2 / 3F).
    field = np.array([0.03, 0.05, 0.08])
    expected = 4 / field * np.exp(-2 / (3 * field))
    assert tunnelling_rate_au(field, 0.5, 4.0, 0) == pytest.approx(
        expected, rel=1e-14, abs=0.0
    )


def test_the_tunnelling_rate_is_tong_zhao_and_lins():
    # Tong et al. (2002), eq. (7), m = 0: C_l^2 (2l + 1) / 2 kappa^(1 - 2n*)
    # (2 kappa^3 / F)^(2n* - 1) exp(-2 kappa^3 / 3F), with C_l^2 = C^2 kappa^(2n*+1).
    energy, c_l, field = 0.57916, 2.44, 0.06
    kappa = math.sqrt(2 * energy)
    n_star = 1 / kappa
    expected = (
        c_l**2
        * 1.5
        * kappa ** (1 - 2 * n_star)
        * (2 * kappa**3 / field) ** (2 * n_star - 1)
        * math.exp(-2 * kappa**3 / (3 * field))
    )
    squared = c_l**2 / kappa ** (2 * n_star + 1)
    rate = tunnelling_rate_au(np.array(field), energy, squared, 1)
    assert rate == pytest.approx(expected, rel=1e-13, abs=0.0)


def test_the_short_range_rate_in_the_tunnelling_limit():
    # Eq. (54) for gamma -> 0 is eq. (59): the cycle average of the static rate
    # without its Coulomb factor, with exp[-(2 F0 / 3F)(1 - gamma^2 / 10)]. At
    # gamma = 1e-3 the exponent's correction is 2 F0 gamma^2 / 30 F = 1e-6.
    field, energy = np.array(0.05), 0.5
    omega = 1e-3 * 0.05
    static = tunnelling_rate_au(field, energy, 4.0, 0) / (2 / 0.05) ** 2
    expected = static * cycle_average_factor(field, energy)
    rate = short_range_rate_au(field, omega, energy, 4.0, 0)
    assert rate == pytest.approx(expected, rel=1e-5, abs=0.0)


def test_the_multiphoton_limit_of_the_short_range_rate():
    # Eq. (60): for gamma >> 1 the rate goes as F^(2K), K the number of photons of
    # the lowest open channel; with the sum taken term by term and the threshold
    # just below an integer the slope d ln w / d ln F tends to 2K from above.
    energy, omega = 0.5, 0.0934  # nu0 = 5.35: six photons
    field = np.array([1e-3, 1.1e-3])
    rate = short_range_rate_au(field, omega, energy, 4.0, 0, channel_sum="sum")
    slope = np.log(rate[1] / rate[0]) / np.log(1.1)
    assert slope == pytest.approx(12, abs=0.2)


def test_the_coulomb_factors():
    field, omega, energy = np.array(0.05), 0.05, 0.5  # gamma = 1, n* = 1
    static = coulomb_factor(field, omega, energy)
    assert static == pytest.approx((2 / 0.05) ** 2)
    ilkov = coulomb_factor(field, omega, energy, "ilkov")
    assert ilkov / static == pytest.approx(0.5)
    popruzhenko = coulomb_factor(field, omega, energy, "popruzhenko")
    assert popruzhenko / static == pytest.approx((1 + 2 / math.e) ** -2)


def test_the_static_corrections_agree_for_argon():
    # Both were fitted to exact one-electron rates of argon relative to the same
    # tunnelling rate (TC-509, TC-655); between 0.64 (the lowest field of Zhang et
    # al.'s Table IX) and 1.7 of the barrier-suppression field they agree within 10 %.
    p = ionisation_parameters("Ar")
    barrier = p.ionization_energy_au**2 / 4
    field = np.linspace(0.64, 1.7, 12) * barrier
    assert p.barrier_suppression_alpha is not None
    assert p.barrier_suppression_zhang is not None
    tong = tong_lin_factor(field, p.ionization_energy_au, p.barrier_suppression_alpha)
    zhang = zhang_factor(field, p.ionization_energy_au, p.barrier_suppression_zhang)
    assert zhang == pytest.approx(tong, rel=0.10, abs=0.0)


def test_zhangs_factor_is_held_at_one():
    p = ionisation_parameters("Ar")
    assert p.barrier_suppression_zhang is not None
    factor = zhang_factor(
        np.array([0.0, 1e-3]), p.ionization_energy_au, p.barrier_suppression_zhang
    )
    assert np.all(factor == 1.0)


# The models.


def test_the_default_covers_the_five_atoms():
    models = ionisation_rate_models()["default"]
    assert {gas.value for gas in models} == set(ATOMS)
    for model in models.values():
        assert model.coulomb == "frequency"
        assert model.static_correction == "tong-lin"
        assert model.structure == "tong_2002"
        assert math.isfinite(model.uncertainty_factor)


def test_the_channels_of_the_default():
    (helium,) = ionisation_rate_model("He", 800e-9).channels
    assert helium.weight == 1.0
    # The ions' 2P1/2 levels (NIST ASD), in eV to the digits that matter here.
    splitting = {"Ne": 0.0968, "Ar": 0.1775, "Kr": 0.6658, "Xe": 1.3064}
    for gas, expected in splitting.items():
        lower, upper = ionisation_rate_model(gas, 800e-9).channels
        assert (lower.weight, upper.weight) == pytest.approx((2 / 3, 1 / 3))
        difference = (upper.ionization_energy - lower.ionization_energy) / E_CHARGE
        assert difference == pytest.approx(expected, abs=5e-5)


@pytest.mark.parametrize("name", ["default", "ppt_ilkov", "ppt_popruzhenko"])
@pytest.mark.parametrize("gas", ATOMS)
def test_a_long_wavelength_gives_the_static_rate(name, gas):
    # At 1 cm the Keldysh parameter is below 1e-4 for the fields here. Popruzhenko
    # et al.'s Coulomb factor, the slowest to converge, differs from its limit by
    # 4 n* gamma / e, below 2e-4; the others by gamma^2.
    static = ionisation_rate_model(gas, math.inf, name)
    slow = ionisation_rate_model(gas, 1e-2, name)
    field = np.array([0.5, 1.0]) * static.barrier_suppression_field
    for method in ("cycle_averaged_rate", "instantaneous_rate"):
        expected = getattr(static, method)(field)
        assert getattr(slow, method)(field) == pytest.approx(
            expected, rel=2e-3, abs=0.0
        )


def test_the_literature_forms_are_the_tunnelling_rate_of_adk():
    # ppt_ilkov and ppt_larochelle at a static field are Ammosov et al.'s rate for
    # m = 0: eq. (1) without its cycle average.
    p = ionisation_parameters("Ar")
    squared = hartree_coefficient_squared(p.n_star, p.l_star)
    field = np.array([0.05, 0.08])
    expected = tunnelling_rate_au(field, p.ionization_energy_au, squared, 1)
    for name in ("ppt_ilkov", "ppt_larochelle"):
        model = ionisation_rate_model("Ar", math.inf, name)
        rate = model.instantaneous_rate(field * AU_FIELD) / AU_ANGULAR_FREQUENCY
        assert rate == pytest.approx(expected, rel=1e-12, abs=0.0)


def _cycle_average(model, amplitude):
    phase = np.linspace(0, math.pi / 2, 4001)
    rate = model.instantaneous_rate(amplitude * np.cos(phase))
    return float(np.trapezoid(rate, phase)) * 2 / math.pi


@pytest.mark.parametrize(
    ("gas", "wavelength"),
    [("Ar", 800e-9), ("Ar", 400e-9), ("He", 800e-9), ("Xe", 1030e-9), ("He", math.inf)],
)
def test_the_instantaneous_rate_averages_to_the_cycle_averaged_rate(gas, wavelength):
    # The relation is exact for a power law and, to first order in 1 / s, for any
    # rate of logarithmic slope s. From a tenth of the barrier-suppression
    # intensity, where s is above 12, to that intensity, where it has fallen to 6,
    # the two agree to 3 % (Keldysh parameters 0.3 to 5).
    model = ionisation_rate_model(gas, wavelength)
    for fraction in (0.32, 0.5, 0.7, 1.0):
        amplitude = fraction * model.barrier_suppression_field
        expected = float(model.cycle_averaged_rate(amplitude))
        assert _cycle_average(model, amplitude) == pytest.approx(
            expected, rel=0.03, abs=0.0
        )


def test_the_tunnelling_average_is_the_literatures():
    # averaging = "tunnelling": the instantaneous rate is the cycle-averaged rate
    # over (3F / pi F0)^(1/2), exactly, and its average over a cycle is then too
    # high: by 11 % for argon at 800 nm and 1e14 W/cm^2 (gamma = 1.15).
    model = ionisation_rate_model("Ar", 800e-9, "ppt_ilkov")
    field = np.array([1.5e10, 2.744e10])
    (channel,) = model.channels
    factor = cycle_average_factor(field / AU_FIELD, channel.ionization_energy_au)
    expected = model.cycle_averaged_rate(field) / factor
    assert model.instantaneous_rate(field) == pytest.approx(
        expected, rel=1e-12, abs=0.0
    )
    ratio = _cycle_average(model, 2.744e10) / float(model.cycle_averaged_rate(2.744e10))
    assert ratio == pytest.approx(1.11, abs=0.01)


def test_ppt_ilkov_is_berge_b14():
    # Bergé et al. (2007), eq. (B14), written out independently, for m = 0:
    # W = (4 sqrt(2) / pi) |C|^2 (2 E0 / (E sqrt(1 + g^2)))^(2n* - 3/2) f(l, 0)
    #     exp(-2 nu [asinh g - g sqrt(1 + g^2) / (1 + 2 g^2)]) Ui g^2 / (1 + g^2)
    #     sum_k exp(-alpha (k - nu)) Phi_0(sqrt(beta (k - nu))),
    # with the sum here replaced by its integral, as the model takes it.
    p = ionisation_parameters("Ar")
    energy, n_star = p.ionization_energy_au, p.n_star
    field = np.array([0.03, 0.05338, 0.08])
    gamma = OMEGA_800 * math.sqrt(2 * energy) / field
    root = np.sqrt(1 + gamma**2)
    nu = (energy + field**2 / (4 * OMEGA_800**2)) / OMEGA_800
    alpha = 2 * (np.arcsinh(gamma) - gamma / root)
    beta = 2 * gamma / root
    integral = 0.5 * math.sqrt(math.pi) / (alpha + beta) * np.sqrt(beta / alpha)
    expected = (
        4
        * math.sqrt(2)
        / math.pi
        * hartree_coefficient_squared(n_star, n_star - 1)
        * (2 * (2 * energy) ** 1.5 / (field * root)) ** (2 * n_star - 1.5)
        * 3
        * np.exp(-2 * nu * (np.arcsinh(gamma) - gamma * root / (1 + 2 * gamma**2)))
        * energy
        * gamma**2
        / (1 + gamma**2)
        * integral
    )
    model = ionisation_rate_model("Ar", 800e-9, "ppt_ilkov")
    rate = model.cycle_averaged_rate(field * AU_FIELD) / AU_ANGULAR_FREQUENCY
    assert rate == pytest.approx(expected, rel=1e-10, abs=0.0)


def test_the_cycle_average_of_the_static_rate():
    # The cycle-averaged rate of a static model is the average of its rate over
    # |cos|, through the local slope: to 2 % at half the barrier-suppression field.
    model = ionisation_rate_model("He", math.inf, "one_electron_static")
    amplitude = 0.1 * AU_FIELD
    expected = float(model.cycle_averaged_rate(amplitude))
    assert _cycle_average(model, amplitude) == pytest.approx(
        expected, rel=0.02, abs=0.0
    )


def test_shape_sign_and_zero():
    model = ionisation_rate_model("Ar", 800e-9)
    field = np.array([[0.0, 2e10], [-2e10, 3e10]])
    rate = model.instantaneous_rate(field)
    assert rate.shape == (2, 2)
    assert rate[0, 0] == 0.0
    assert rate[0, 1] == rate[1, 0] > 0
    assert model.cycle_averaged_rate(2e10).shape == ()


def test_the_rate_rises_with_the_field_and_the_frequency():
    field = np.geomspace(5e9, 4e10, 40)
    previous = None
    for wavelength in (math.inf, 1800e-9, 800e-9, 400e-9):
        rate = ionisation_rate_model("Ar", wavelength).cycle_averaged_rate(field)
        assert np.all(np.diff(rate) > 0)
        if previous is not None:
            assert np.all(rate > previous)
        previous = rate


def test_the_effective_electrons_multiply_the_rate():
    model = ionisation_rate_model("He", 800e-9)
    one = dataclasses.replace(model, electrons=1.0)
    rate = model.cycle_averaged_rate(6e10) / one.cycle_averaged_rate(6e10)
    assert rate == pytest.approx(model.electrons)


def test_the_fine_structure_lowers_the_rate():
    # The 2P1/2 channel has the higher ionisation energy, so the rate is below that
    # of three thirds of the 2P3/2 channel, and no lower than two thirds of it.
    for gas in ("Ne", "Ar", "Kr", "Xe"):
        model = ionisation_rate_model(gas, 800e-9)
        lower = model.channels[0]
        single = dataclasses.replace(
            model, channels=(dataclasses.replace(lower, weight=1.0),)
        )
        ratio = model.cycle_averaged_rate(2e10) / single.cycle_averaged_rate(2e10)
        assert 2 / 3 < ratio < 1


def test_the_keldysh_parameter():
    model = ionisation_rate_model("Ar", 800e-9)
    energy = model.channels[0].ionization_energy_au
    expected = keldysh_parameter(np.array(2e10 / AU_FIELD), OMEGA_800, energy)
    assert model.keldysh(2e10) == pytest.approx(expected, rel=1e-6, abs=0.0)
    assert model.keldysh(0.0) == math.inf
    assert ionisation_rate_model("Ar", math.inf).keldysh(2e10) == 0.0


def test_warnings_outside_the_ranges():
    model = ionisation_rate_model("Ar", 800e-9)
    with pytest.warns(ModelRangeWarning, match="static correction"):
        model.instantaneous_rate(1.01 * model.maximum_field)
    with pytest.warns(ModelRangeWarning, match="shortest wavelength"):
        ionisation_rate_model("Ar", 140e-9)


def test_errors():
    with pytest.raises(DataError, match="no ionisation-rate model"):
        ionisation_rate_model("Ar", 800e-9, "other")
    with pytest.raises(DataError, match="does not cover"):
        ionisation_rate_model("N2", 800e-9)
    model = ionisation_rate_model("Ar", 800e-9)
    with pytest.raises(ValueError, match="coulomb"):
        dataclasses.replace(model, coulomb="other")
    with pytest.raises(ValueError, match="electrons"):
        dataclasses.replace(model, electrons=0.0)
    with pytest.raises(ValueError, match="weights"):
        dataclasses.replace(model, channels=model.channels[:1])
    with pytest.raises(ValueError, match="constants"):
        dataclasses.replace(model, static_constants=())


def test_the_citation_names_every_factor():
    citation = ionisation_rate_model("Ar", 800e-9).citation
    for name in ("Perelomov", "Larochelle", "Tong & Lin", "Tong et al."):
        assert name in citation


# Against the references.


def test_argon_one_electron_static_rate_against_zhangs_exact_rates():
    # Zhang et al. (2014), Table IX, Tong and Lin's potential: the exact rates the
    # static correction was fitted to reproduce. Within 8 % up to 1.72 times the
    # barrier-suppression field (0.144 au); the other two potentials, whose rates
    # differ from these by up to 25 %, within 21 %.
    model = ionisation_rate_model("Ar", math.inf, "one_electron_static")
    for potential, tolerance in (("tonglin", 0.08), ("gsz", 0.21), ("muller", 0.21)):
        dataset = _dataset(f"zhang_2014/tableIX_{potential}_Ar")
        field = dataset.abscissa_values()
        keep = field <= 0.145 * AU_FIELD
        rate = model.instantaneous_rate(field[keep]) / AU_ANGULAR_FREQUENCY
        assert rate == pytest.approx(
            dataset.scaled_values()[keep], rel=tolerance, abs=0.0
        )


def test_helium_static_rate_against_the_two_electron_calculations():
    # Parker et al. (2009), Table 1, with which four other all-electron calculations
    # agree to 4 % above 0.09 au. The default's effective number of electrons is the
    # geometric mean of the ratio over 0.07 to 0.25 au, where the ratio itself runs
    # from 0.86 to 1.08 of the mean (TC-658): within 15 %.
    model = ionisation_rate_model("He", math.inf)
    dataset = _dataset("parker_2009/table1_present")
    field = dataset.abscissa_values()
    keep = (field >= 0.0699 * AU_FIELD) & (field <= 0.2501 * AU_FIELD)
    rate = model.instantaneous_rate(field[keep]) / AU_ANGULAR_FREQUENCY
    ratio = dataset.scaled_values()[keep] / rate
    assert ratio == pytest.approx(1.0, abs=0.15)
    assert float(np.exp(np.mean(np.log(ratio)))) == pytest.approx(1.0, abs=0.005)


def _dense_average(log_terms, field_au):
    # The average over a cycle by 2000-point Gauss-Legendre in the phase itself.
    nodes, weights = np.polynomial.legendre.leggauss(2000)
    phase = (nodes + 1) * math.pi / 4
    averages = []
    for value in field_au:
        peak = log_terms(np.array([value]))
        cycle = log_terms(value * np.cos(phase))
        averages.append(
            [
                np.sum(weights / 2 * np.exp(c - p[0]))
                for c, p in zip(cycle, peak, strict=True)
            ]
        )
    return np.array(averages).T


@pytest.mark.parametrize("gas", ["He", "Ar", "Xe"])
def test_the_cycle_average_by_quadrature(gas):
    # The average of each rate over a cycle within 1e-6 of a dense quadrature up to
    # the maximum field and 1e-4 at three times it, at 400 nm and in a static field,
    # where the expansion in slope and curvature was 24 % low (argon, at the maximum
    # field) (TC-828).
    from optmat.gas.ionisation_rates import (  # noqa: PLC0415 - private
        _log_cycle_averages,  # pyright: ignore[reportPrivateUsage]
    )

    for wavelength in (400e-9, math.inf):
        model = ionisation_rate_model(gas, wavelength)
        x = np.array([0.01, 0.1, 0.5, 1.0, 3.0])
        field = x * model.maximum_field / AU_FIELD
        computed = np.exp(np.array(_log_cycle_averages(model._log_terms, field)))
        reference = _dense_average(model._log_terms, field)
        error = np.abs(computed / reference - 1)
        assert np.all(error[:, x <= 1] < 1e-6)
        assert np.all(error < 1e-4)


def test_the_log_sum_exp_is_scipys_to_the_last_bit():
    # The private reduction repeats scipy.special.logsumexp(axis=0)'s arithmetic, so
    # that replacing the call (for its dispatch, 70 us) changes no rate: ties, one
    # row, infinities and NaN included.
    from scipy import special  # noqa: PLC0415 - the reference only

    from optmat.gas.ionisation_rates import (  # noqa: PLC0415 - private
        _log_sum_exp,  # pyright: ignore[reportPrivateUsage]
    )

    rng = np.random.default_rng(1)
    for rows in (1, 2, 3, 9):
        terms = rng.normal(scale=200.0, size=(rows, 500))
        terms[:, :20] = terms[0, :20]
        terms[-1, 20:40] = -np.inf
        terms[:, 40] = -np.inf
        terms[0, 41] = np.inf
        terms[-1, 42] = np.nan
        terms[:, 43] = 800.0
        with np.errstate(all="ignore"):
            expected = special.logsumexp(terms, axis=0)
        np.testing.assert_array_equal(_log_sum_exp(terms), expected)


@pytest.mark.parametrize("wavelength", [400e-9, 800e-9, math.inf])
def test_each_rate_alone_is_the_pair(wavelength):
    # The average over a cycle evaluates each rate alone at its own nodes, the
    # expensive one only at its own: the same numbers, to the last bit, as the pair.
    from optmat.gas.ionisation_rates import (  # noqa: PLC0415 - private
        _log_cycle_averages,  # pyright: ignore[reportPrivateUsage]
    )

    for model in (
        ionisation_rate_model("Ar", wavelength),
        molecular_rate_model("N2", wavelength),
        molecular_rate_model("H2", wavelength),
    ):
        field = np.geomspace(1e-3, 2.0, 50) * model.maximum_field / AU_FIELD
        pair = model._log_terms(field)
        for k in (0, 1):
            np.testing.assert_array_equal(model._log_term(k, field), pair[k])
        averages = _log_cycle_averages(model._log_terms, field)
        given = _log_cycle_averages(
            model._log_terms, field, centre=pair, log_term=model._log_term
        )
        for a, b in zip(averages, given, strict=True):
            np.testing.assert_array_equal(a, b)


def test_the_closed_forms_are_the_functions_of_gamma():
    # The rates' shared evaluation of the closed forms is channel_factor's integral
    # and ppt_exponent_function to the last bit, from the series' end up.
    from optmat.gas.ionisation_rates import (  # noqa: PLC0415 - private
        _closed_forms,  # pyright: ignore[reportPrivateUsage]
    )

    gamma = np.geomspace(1e-2, 1e3, 2001)
    for m in (0, 1, 2):
        factor, exponent = _closed_forms(m, gamma, gamma**2)
        np.testing.assert_array_equal(factor, channel_factor(m, gamma))
        np.testing.assert_array_equal(exponent, ppt_exponent_function(gamma))


@pytest.mark.parametrize("gas", ["He", "Ne", "Ar", "Kr", "Xe"])
def test_the_rates_rise_and_are_held_above_the_range(gas):
    # The two rates rise with the field to the maximum field and are held there
    # above it; they turned over at 1.01 to 1.19 times it before (TC-816, TC-828).
    model = ionisation_rate_model(gas, 800e-9)
    field = np.linspace(0.3, 3.0, 541) * model.maximum_field
    with pytest.warns(ModelRangeWarning, match="held at its value there"):
        averaged = model.cycle_averaged_rate(field)
        instantaneous = model.instantaneous_rate(field)
    for rate in (averaged, instantaneous):
        assert np.all(np.diff(rate) >= 0)
        above = field > model.maximum_field
        assert np.all(rate[above] == rate[above][0])
    at = model.instantaneous_rate(model.maximum_field)
    assert instantaneous[above][0] == pytest.approx(at, rel=1e-12, abs=0.0)
    # A model without a correction keeps its formula above its range.
    literature = ionisation_rate_model(gas, 800e-9, "ppt_ilkov")
    with pytest.warns(ModelRangeWarning, match="static correction"):
        rate = literature.instantaneous_rate(field)
    assert np.all(np.diff(rate) > 0)


def test_the_rate_is_smooth_near_the_maximum_field():
    # The rounding of a second difference left noise of 3.6e-9 in helium's
    # instantaneous rate at 0.9 of the maximum field (TC-828).
    model = ionisation_rate_model("He", 800e-9)
    field = np.linspace(0.9, 0.9009, 201) * model.maximum_field
    t = np.log(field) - np.log(field).mean()
    for rate in (model.instantaneous_rate, model.cycle_averaged_rate):
        log_rate = np.log(rate(field))
        residual = log_rate - np.polyval(np.polyfit(t, log_rate, 4), t)
        assert np.max(np.abs(residual)) < 1e-12

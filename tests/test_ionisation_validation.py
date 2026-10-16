# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The ionisation rates against their references: the bounds the manual states."""

import dataclasses
import math

import numpy as np
import pytest

from optmat.constants import AU_ANGULAR_FREQUENCY
from optmat.datasets import datasets
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import ionisation_validation as validation
from optmat.gas.ionisation_rates import ionisation_rate_model

# The pulses and the averages.


def test_the_amplitude_of_an_intensity():
    # 1e14 W/cm^2 is 0.05338 au of field (3.509e16 W/cm^2 per squared atomic unit).
    field = validation.amplitude(1e18) / 5.14220674763e11
    assert field == pytest.approx(math.sqrt(1e14 / 3.50944e16), rel=1e-5, abs=0.0)


def test_the_gaussian_pulse():
    times, envelope, field = validation.gaussian_pulse(1e18, 6e-15, 800e-9)
    peak = float(validation.amplitude(1e18))
    assert np.max(envelope) == pytest.approx(peak, rel=1e-12, abs=0.0)
    # The full width at half maximum of the intensity is the one asked for, to the
    # grid's spacing of a hundredth of a cycle.
    above = times[envelope**2 >= 0.5 * peak**2]
    assert above[-1] - above[0] == pytest.approx(6e-15, abs=6e-17)
    # A field that is the derivative of a vector potential has no static component.
    assert abs(np.trapezoid(field, times)) < 1e-6 * np.trapezoid(np.abs(field), times)
    # With a phase of 0 the crest of the field is at the peak of the envelope.
    assert times[np.argmax(np.abs(field))] == pytest.approx(0.0, abs=1e-18)


def test_the_ionised_fraction():
    times = np.linspace(0, 1e-14, 101)
    assert validation.ionised_fraction(np.zeros(101), times) == 0.0
    fraction = validation.ionised_fraction(np.full(101, 1e14), times)
    assert fraction == pytest.approx(1 - math.exp(-1), rel=1e-12, abs=0.0)


def test_the_focal_average_of_a_power_law():
    # int_0^I0 I^m dI / I = I0^m / m, exactly; the trapezoid rule in the logarithm
    # of the intensity, with steps h = 0.017, is good to (m h)^2 / 12 = 9e-4 for m = 6.
    intensity = np.geomspace(1e-3, 1.0, 400)
    for power in (2, 6):
        average = validation.focal_average(intensity, intensity**power)
        assert average[-1] == pytest.approx(1 / power, rel=2e-3, abs=0.0)


@pytest.mark.parametrize(("gas", "wavelength"), [("Ar", 800e-9), ("He", 800e-9)])
def test_the_two_rates_give_the_same_yield(gas, wavelength):
    # The instantaneous rate on the carrier and the cycle-averaged rate on the
    # envelope, for a 42 fs pulse, whose envelope is constant over a cycle to 1 %:
    # the yields agree as the rates' cycle averages do (3 %).
    model = ionisation_rate_model(gas, wavelength)
    peak = (
        0.5
        * 8.8541878128e-12
        * 299792458
        * (0.6 * model.barrier_suppression_field) ** 2
    )
    on_envelope = validation.pulse_yield(model, peak, 42e-15)
    on_carrier = validation.pulse_yield(model, peak, 42e-15, carrier=True)
    assert on_carrier == pytest.approx(on_envelope, rel=0.03, abs=0.0)


# The references. The bounds are those the manual page states.


def test_helium_static_rates():
    parker, scrinzi, majety, jagau = validation.static_rates("He")
    field = parker.x / 5.14220674763e11
    core = (field >= 0.0699) & (field <= 0.1501)
    assert parker.ratio[core] == pytest.approx(1.0, abs=0.08)
    wide = (field >= 0.0699) & (field <= 0.3001)
    assert parker.ratio[wide] == pytest.approx(1.0, abs=0.18)
    # Scrinzi et al.'s lowest field, 0.08 au, is 25 % below Parker et al.'s and
    # Majety and Scrinzi's values there; above it the three agree to 3 %.
    assert scrinzi.ratio[1:] == pytest.approx(1.0, abs=0.20)
    assert majety.ratio == pytest.approx(1.0, abs=0.11)
    # Jagau's full configuration interaction in a Gaussian basis, 0.10 to 0.40 au: its
    # own rate at 0.10 au is 19 % below Scrinzi et al.'s (its Table S2).
    assert jagau.ratio == pytest.approx(1.0, abs=0.21)


def test_neon_static_rates():
    jagau, _bahl = validation.static_rates("Ne")
    assert jagau.model == "default"
    field = jagau.x / 5.14220674763e11
    # From 0.14 au, where the reference is smooth in the field (TC-669).
    assert jagau.ratio[field >= 0.1399] == pytest.approx(1.0, abs=0.19)
    # Below, the reference rises above the model, by up to 2.7 at 0.08 au.
    assert np.all(jagau.ratio[field < 0.1399] < 0.85)


def test_argon_static_rates():
    jagau, tong_lin, gsz, muller, trecx, _bahl = validation.static_rates("Ar")
    assert jagau.model == "default"
    field = jagau.x / 5.14220674763e11
    # The atom, from 0.07 au to twice the barrier-suppression field (0.16 au).
    assert jagau.ratio[field >= 0.0699] == pytest.approx(1.0, abs=0.23)
    # Below, the atom's rate rises above the model's, toward the weak-field limit of
    # two electrons (and beyond it at 0.05 au, where the reference is in doubt).
    assert np.all(jagau.ratio[field < 0.0699] < 0.65)
    assert tong_lin.model == "one_electron"
    assert tong_lin.ratio == pytest.approx(1.0, abs=0.07)
    assert gsz.ratio == pytest.approx(1.0, abs=0.21)
    assert muller.ratio == pytest.approx(1.0, abs=0.21)
    # tRecX's rates of the same electron, the excluded field of Zhang et al. included
    # (TC-813).
    assert trecx.model == "one_electron"
    assert trecx.ratio == pytest.approx(1.0, abs=0.07)


def test_the_tdse_static_rates_reproduce_zhang_et_al():
    # 0.04 to 2.4 % at the six fields Zhang et al.'s Table IX keeps; at the seventh,
    # Gamma/2 is the printed Muller cell (TC-813).
    by_key = {d.key: d for d in datasets()}
    trecx = by_key["travers_2026_trecx/static_Ar_m0"]
    zhang = by_key["zhang_2014/tableIX_tonglin_Ar"]
    calculated = dict(zip(trecx.abscissa_values(), trecx.scaled_values(), strict=True))
    printed = dict(zip(zhang.abscissa_values(), zhang.scaled_values(), strict=True))
    common = [f for f in printed if any(np.isclose(f, c) for c in calculated)]
    assert len(common) == 6
    for f in common:
        (match,) = [c for c in calculated if np.isclose(f, c)]
        assert calculated[match] / printed[f] == pytest.approx(1.0, abs=0.025)
    muller = dict(
        zip(
            by_key["zhang_2014/tableIX_muller_Ar"].abscissae,
            by_key["zhang_2014/tableIX_muller_Ar"].values,
            strict=True,
        )
    )
    assert "0.12177" not in muller  # excluded, as the Tong-Lin row
    (row,) = [
        v for f, v in calculated.items() if np.isclose(f, 0.12177 * 5.14220674763e11)
    ]
    assert row / 2 == pytest.approx(0.013229, rel=5e-4, abs=0.0)


def test_parkers_fit_reproduces_the_papers_values():
    # Parker et al. (2007), p. 13, as the source record checks them: the fits give
    # 0.0107 au at 30e14 W/cm^2 (390 nm) and 0.00142 au at 14.5e14 W/cm^2 (780 nm).
    assert validation.parker_rate_au(30e18, 390) == pytest.approx(0.0107, abs=5e-5)
    assert validation.parker_rate_au(14.5e18, 780) == pytest.approx(0.00142, abs=5e-6)
    with pytest.raises(DataError, match="390 and 780"):
        validation.parker_rate_au(1e18, 800)


def test_helium_in_a_laser_field():
    # With the static factor (TC-660, TC-828):
    infrared, blue = validation.helium_laser_rates("default_static")
    assert infrared.ratio == pytest.approx(1.07, abs=0.13)  # 0.95 to 1.20
    assert blue.ratio == pytest.approx(1.36, abs=0.10)  # 1.27 to 1.45
    # The default (TC-848): 0.85 to 1.13 at 780 nm, 0.85 to 1.28 at 390 nm.
    infrared, blue = validation.helium_laser_rates()
    assert np.all((infrared.ratio > 0.84) & (infrared.ratio < 1.14))
    assert np.all((blue.ratio > 0.84) & (blue.ratio < 1.29))


def test_the_factor_of_ilkov_is_too_low_for_helium():
    # With (1 + gamma^2)^(-n*) the rate at 780 nm falls below Parker et al.'s by
    # more than 2 at 1e14 W/cm^2 (gamma = 1.47), and by 4 at 390 nm.
    infrared, blue = validation.helium_laser_rates("ppt_ilkov")
    assert infrared.ratio[0] < 0.45
    assert blue.ratio[0] < 0.25


@pytest.mark.parametrize("gas", ["Ar", "He"])
def test_the_one_electron_yields(gas):
    ionised = validation.tdse_yields(gas, "one_electron_static")[0]
    # Against the ionised fraction the model is 1.02 to 1.67 times (argon), the most at
    # the largest Keldysh parameter, and 1.16 to 1.23 times (helium) (TC-812, TC-828).
    assert np.all(ionised.ratio > 1.0)
    assert np.all(ionised.ratio < (1.7 if gas == "Ar" else 1.25))


def test_the_other_coulomb_factors_against_the_one_electron_yields():
    for name in ("ppt_ilkov", "ppt_popruzhenko"):
        beyond, *_ = validation.tdse_yields("Ar", name, carrier=False)
        gamma_above_two = beyond.x < 3e17  # below 3e13 W/cm^2
        assert np.all(beyond.ratio[gamma_above_two] < 0.4)


@pytest.mark.parametrize(
    ("gas", "low", "high"), [("Ar", 1.6, 4.0), ("Kr", 1.6, 4.0), ("Xe", 4.0, 8.0)]
)
def test_wahlstrands_yields_are_below_the_models(gas, low, high):
    # The open disagreement of TC-660, recorded as it stands: a change of the model
    # or of the reference that moves it fails here and is looked at.
    comparison = validation.measured_yields(gas, "default_static")
    assert np.all(comparison.ratio > low)
    assert np.all(comparison.ratio < high)


@pytest.mark.parametrize("gas", ["Ar", "Kr", "Xe"])
def test_wallaces_yields_are_matched_within_their_calibration(gas):
    test = validation.intensity_scale_test(gas, "default_static")
    # Two standard uncertainties of the calibration, and the step of the search
    # for the scale (0.2 %).
    assert test.scale == pytest.approx(1.0, abs=2 * test.calibration + 0.002)
    assert test.residual < 0.02


def test_kryptons_measured_points():
    test = validation.measured_points_test("Kr", "default_static")
    assert len(test.measured) == 37
    # The calibration of the intensities is 1.5 %.
    assert test.scale == pytest.approx(1.0, abs=0.015)
    # The authors' fit function has 369 (wallace_2016, the formula's check); the
    # printed errors are smaller than the points' scatter, so no model has 37.
    assert test.chi_square < 300
    # 4 % root mean square, of which the point at 3.49e14 W/cm^2 is 15 %.
    assert test.spread < 0.04
    assert np.max(np.abs(test.predicted / test.measured - 1)) < 0.16
    # One electron without the fine structure of the ion needs 4 % in intensity.
    assert validation.measured_points_test("Kr", "one_electron_static").scale < 0.97
    # The literature's forms have twice the chi-square.
    assert validation.measured_points_test("Kr", "ppt_ilkov").chi_square > 550


def test_the_factor_of_ilkov_fails_wallaces_test():
    assert validation.intensity_scale_test("Ar", "ppt_ilkov").scale > 1.1


def test_errors():
    with pytest.raises(DataError, match="no static"):
        validation.static_rates("N2")
    with pytest.raises(DataError, match="no TDSE"):
        validation.tdse_yields("Ne")
    with pytest.raises(DataError, match="no yield"):
        validation.measured_yields("He")
    with pytest.raises(DataError, match="no yield"):
        validation.wallace_yield("He", 1e18)
    with pytest.raises(DataError, match="no measured points"):
        validation.measured_points_test("Ar")


# The one-electron calculation of argon in the 42 fs pulse (TC-702).


def test_argon_in_the_42_fs_pulse_against_the_models_and_the_measurement():
    # The one_electron model is 0.99 to 1.8 times the calculated ionised fraction,
    # 1.24 on geometric average (the calculation's resonance structure); the
    # calculation, one electron, is 1.3 to 4.1 times Wahlstrand et al.'s measured yield
    # of the atom.
    one = validation.long_pulse_yield("one_electron_static")
    assert len(one.x) == 17
    assert np.all((one.ratio > 0.95) & (one.ratio < 1.85))
    assert float(np.exp(np.mean(np.log(one.ratio)))) == pytest.approx(1.24, abs=0.02)
    measured = validation.measured_yields(
        "Ar", "default_static", intensities=list(one.x)
    )
    over = one.expected / measured.expected
    assert np.all((over > 1.25) & (over < 4.1))


# The all-electron R-matrix rates at 390 nm (TC-700, TC-701).


def _geometric_mean(values):
    return float(np.exp(np.mean(np.log(values))))


def test_the_floquet_curve_read_from_an_image_matches_its_vector_copy():
    # van der Hart 2006's Ar curve, read from the image of its Fig. 6, against the same
    # curve as van der Hart 2007's Fig. 4 draws it in vectors: within 4 % where the
    # vector curve has vertices within 0.004e14 W/cm2 (data rule 8, TC-701).
    by_key = {c.reference: c for c in validation.rmatrix_rates("Ar")}
    raster = by_key["vanderhart_2006/fig6_floquet_Ar"]
    vector = next(d for d in datasets() if d.key == "vanderhart_2007/fig4_floquet_Ar")
    x, rate = vector.abscissa_values(), vector.scaled_values() / AU_ANGULAR_FREQUENCY
    compared = 0
    for intensity, expected in zip(raster.x, raster.expected, strict=True):
        near = np.abs(x - intensity) < 0.004e18
        if near.any():
            copy = float(np.exp(np.median(np.log(rate[near]))))
            assert expected / copy == pytest.approx(1.0, abs=0.04)
            compared += 1
    assert compared >= 20


def test_helium_at_390_nm_against_the_floquet_rates():
    # Parker et al.'s two-electron rate over van der Hart's Floquet rate: 0.8 on
    # geometric average, 0.96 to 1.25 away from the resonance at 1.5e14 W/cm2; the
    # default, with its 1.47 electrons, 1.13.
    (floquet,) = validation.rmatrix_rates("He", "default_static")
    parker = validation.parker_rate_au(floquet.x, 390) / floquet.expected
    assert _geometric_mean(parker) == pytest.approx(0.8, abs=0.05)
    assert _geometric_mean(floquet.ratio) == pytest.approx(1.13, abs=0.05)


def test_argon_at_390_nm_against_the_rmatrix_rates():
    # The default is 1.6 to 2.9 times the all-electron rates from 1.25e14 W/cm2 (a
    # Keldysh parameter of 2.1 to 1.4), and 7 to 12 times the time-dependent ones at
    # 0.25 to 0.75e14 W/cm2 (TC-701, open).
    for comparison in validation.rmatrix_rates("Ar", "default_static"):
        high = comparison.x >= 1.24e18
        assert np.all((comparison.ratio[high] > 1.5) & (comparison.ratio[high] < 3.0))
        low = comparison.x <= 0.76e18
        if comparison.reference != "vanderhart_2006/fig6_floquet_Ar":
            assert np.all(comparison.ratio[low] > 7.0)
    with pytest.raises(DataError, match="no R-matrix rates"):
        validation.rmatrix_rates("Kr")


# Plasma interferometry: Nambu et al. 2026 (TC-817).


def test_the_phase_of_a_fully_ionised_column():
    beams = {
        "wavelength_nm": 800,
        "probe_wavelength_nm": 400,
        "m2": 1.25,
        "pump_waist_um": 60.0,
        "probe_waist_um": 42.0,
        "density_m3": 2.7e25,
    }
    # Nie et al.'s helium at 5e14 W/cm2 (Fig. 3b, c): a peak ionised fraction of
    # about 5.5e-3, the yield's slope about 7, gives their averaged 0.27 rad to the
    # precision those are read; spot sizes taken as diameters would give a quarter.
    grid = np.geomspace(1e16, 1e20, 400)
    power_law = np.minimum(5.5e-3 * (grid / 5e18) ** 7, 1.0)
    phase = float(validation.plasma_phase_shift(power_law, grid, 5e18, beams))
    assert phase == pytest.approx(0.27, rel=0.15, abs=0.0)


def test_the_envelope_yields_are_the_pulse_yields():
    model = ionisation_rate_model("Ar", 800e-9)
    peaks = np.array([1e18, 2e18, 4e18])
    envelope = validation._envelope_yields(model, peaks, 55e-15)
    direct = [validation.pulse_yield(model, peak, 55e-15) for peak in peaks]
    assert envelope == pytest.approx(direct, rel=1e-6, abs=0.0)
    static = ionisation_rate_model("Ar", math.inf)
    # The turning field's magnitude is the amplitude over the square root of 2.
    times = np.linspace(-3, 3, 20001) * 55e-15
    magnitude = (
        validation.amplitude(4e18)
        / math.sqrt(2)
        * np.exp(-2 * math.log(2) * (times / 55e-15) ** 2)
    )
    expected = -math.expm1(
        -float(np.trapezoid(static.instantaneous_rate(magnitude), times))
    )
    assert validation.circular_yield(static, 4e18, 55e-15) == pytest.approx(
        expected, rel=1e-5, abs=0.0
    )


@pytest.mark.parametrize(
    ("gas", "low", "high"), [("Ar", 0.97, 1.05), ("Kr", 0.92, 1.0), ("Xe", 0.85, 0.93)]
)
def test_the_default_matches_nambus_absolute_yields_at_800_nm(gas, low, high):
    # The scale on the measurement's intensities at which the default's phase fits:
    # 1.01 (Ar), 0.96 (Kr), 0.89 (Xe), against Wahlstrand et al.'s 1.7 to 3.9 times
    # in the yield (TC-660, TC-817).
    test = validation.phase_shift_test(f"nambu_2026/fig6_{gas}_LP", "default_static")
    assert low < test.scale < high


@pytest.mark.parametrize("gas", ["Ar", "Kr", "Xe"])
def test_the_default_is_too_high_at_400_nm(gas):
    # The static Coulomb factor at Keldysh parameters 2 to 6: 2.9 to 5.2 times the
    # measured phase on geometric average, scales 0.50 to 0.69 (TC-817).
    test = validation.phase_shift_test(f"nambu_2026/fig7_{gas}_LP", "default_static")
    assert test.scale < 0.72
    assert float(np.exp(np.mean(np.log(test.ratio)))) > 2.5
    # Popruzhenko et al.'s factor is within 20 % in the intensity.
    other = validation.phase_shift_test(f"nambu_2026/fig7_{gas}_LP", "ppt_popruzhenko")
    assert 1.0 < other.scale < 1.2


def test_a_phase_shift_test_needs_a_phase_shift():
    with pytest.raises(DataError, match="not an ionization_phase_shift"):
        validation.phase_shift_test("wahlstrand_2018/yield_fit_c1_Ar")
    # Nie et al. state the beam quality of the 800 nm laser only.
    with pytest.raises(DataError, match="does not state m2"):
        validation.phase_shift_test("nie_2022/fig6_H2_400nm")


def _intensity_at_equal_phase(reproduced, original):
    by_key = {d.key: d for d in datasets()}
    a, b = by_key[original], by_key[reproduced]
    order = np.argsort(a.scaled_values())
    phase, intensity = b.scaled_values(), b.abscissa_values()
    keep = (phase >= a.scaled_values().min()) & (phase <= a.scaled_values().max())
    at = np.exp(
        np.interp(
            np.log(phase[keep]),
            np.log(a.scaled_values()[order]),
            np.log(a.abscissa_values()[order]),
        )
    )
    return float(np.median(intensity[keep] / at))


@pytest.mark.parametrize(
    ("gas", "low", "high"), [("Ar", 0.86, 0.91), ("Kr", 0.86, 0.91), ("Xe", 0.80, 0.85)]
)
@pytest.mark.parametrize("polarisation", ["LP", "CP"])
def test_nambus_reproduction_of_nies_points_lowers_the_intensities(
    gas, low, high, polarisation
):
    # At equal phase, 0.88 to 0.89 (Ar, Kr) and 0.82 (Xe) times Nie et al.'s
    # intensities (TC-818).
    ratio = _intensity_at_equal_phase(
        f"nambu_2026/fig6_{gas}_{polarisation}", f"nie_2022/fig4_{gas}_{polarisation}"
    )
    assert low < ratio < high


@pytest.mark.parametrize("source", ["nambu_2026/fig6", "nie_2022/fig4"])
def test_the_absolute_measurements_put_xenon_high_against_argon(source):
    # One intensity calibration serves the three gases of a measurement: xenon's
    # scale is 0.82 to 0.88 of argon's, krypton's 0.94 to 0.98 (TC-819).
    scale = {
        gas: validation.phase_shift_test(f"{source}_{gas}_LP", "default_static").scale
        for gas in ("Ar", "Kr", "Xe")
    }
    assert 0.8 < scale["Xe"] / scale["Ar"] < 0.9
    assert 0.92 < scale["Kr"] / scale["Ar"] < 1.0


@pytest.mark.parametrize(
    ("gas", "low", "high"),
    [("Ne", 1.05, 1.36), ("Ar", 0.9, 1.13), ("Kr", 1.08, 1.23), ("Xe", 0.86, 1.06)],
)
def test_the_one_electron_rate_against_bahls_static_rates(gas, low, high):
    # One electron, m = 0, one method for every gas (TC-821): the model's static
    # limit is within 0.86 to 1.36 of it, krypton and xenon included; the excess of
    # xenon against argon in the measurements is not in the static limit.
    (bahl,) = [c for c in validation.static_rates(gas) if "bahl_2016" in c.reference]
    assert np.all((bahl.ratio > low) & (bahl.ratio < high))


def test_krypton_against_the_spin_orbit_resolved_static_rate():
    # Santra et al.'s 4p total at 0.10 au: the default is 0.77 of it (TC-821).
    (santra,) = [
        c for c in validation.static_rates("Kr") if "santra_2006" in c.reference
    ]
    # 0.77 with argon's 1 electron (TC-821); with krypton's 1.35 (TC-848), 1.04.
    assert float(santra.ratio[0]) == pytest.approx(1.04, abs=0.02)


def test_the_one_electron_yields_at_400_nm():
    # 10 fs at 400 nm (TC-824): the default's static Coulomb factor is 1.4 to 4.0 times
    # the calculation at Keldysh parameters 1.4 to 4.2, the measured excess of Nambu et
    # al. at 400 nm (TC-817) for one electron; the literature's form with Popruzhenko
    # et al.'s factor is 0.59 to 0.90 of it above 2, on the carrier.
    one = validation.tdse_yields("Ar", "one_electron_static")[1]
    assert one.wavelength == pytest.approx(400e-9)
    gamma = ionisation_rate_model("Ar", 400e-9, "one_electron_static").keldysh(
        validation.amplitude(one.x)
    )
    large = gamma > 1.3
    assert np.all((one.ratio[large] > 1.25) & (one.ratio[large] < 4.1))
    other = validation.tdse_yields("Ar", "ppt_popruzhenko")[1]
    assert np.all(other.ratio[gamma > 2] < 0.95)


def test_the_one_electron_yields_at_1030_nm():
    # 10 fs at 1030 nm (TC-826): the default's static Coulomb factor is 1.04 to 1.35
    # times the calculation at Keldysh parameters 0.5 to 1.8, much less than at 400 nm
    # with the same pulse at the same Keldysh parameter; the literature's form with
    # Popruzhenko et al.'s factor is 0.40 to 0.71 of it above 1.
    one = validation.tdse_yields("Ar", "one_electron_static")[2]
    assert one.wavelength == pytest.approx(1030e-9)
    assert one.x.size == 8  # rows 0 and 1 are below the floor
    gamma = ionisation_rate_model("Ar", 1030e-9, "one_electron_static").keldysh(
        validation.amplitude(one.x)
    )
    unsaturated = one.expected < 0.5
    assert np.all((one.ratio[unsaturated] > 1.0) & (one.ratio[unsaturated] < 1.4))
    short = validation.tdse_yields("Ar", "one_electron_static")[1]
    short_gamma = ionisation_rate_model("Ar", 400e-9, "one_electron_static").keldysh(
        validation.amplitude(short.x)
    )
    at = 1.8
    assert np.interp(at, short_gamma[::-1], short.ratio[::-1]) > 1.5 * np.interp(
        at, gamma[::-1], one.ratio[::-1]
    )
    other = validation.tdse_yields("Ar", "ppt_popruzhenko")[2]
    assert np.all(other.ratio[gamma > 1] < 0.75)


def test_the_one_electron_yields_with_one_pulse_at_three_wavelengths():
    # The 10 fs pulse at 400, 800 and 1030 nm (TC-824, TC-830, TC-826): at gamma = 1.8
    # the default's static factor is 1.9, 1.7 and 1.1 times the calculation, so the
    # excess grows with the frequency at one pulse length.
    by_wavelength = {
        round(c.wavelength * 1e9): c
        for c in validation.tdse_yields("Ar", "one_electron_static")
        if "gauss10" in c.reference
    }
    assert sorted(by_wavelength) == [400, 800, 1030]
    excess = {}
    for nm, comparison in by_wavelength.items():
        gamma = ionisation_rate_model("Ar", nm * 1e-9, "one_electron_static").keldysh(
            validation.amplitude(comparison.x)
        )
        excess[nm] = float(
            np.exp(np.interp(1.8, gamma[::-1], np.log(comparison.ratio[::-1])))
        )
    assert excess[400] > excess[800] > excess[1030]
    assert excess[800] == pytest.approx(1.7, abs=0.1)
    assert excess[1030] == pytest.approx(1.1, abs=0.1)


def test_the_frequency_candidate():
    # The frequency-dependent Coulomb factor (TC-831), fitted to the 10 fs sets: within
    # 10 % of every one-electron argon calculation on geometric average (the static
    # factor: 1.17 to 2.2), including the 6 fs set it was not fitted to; Wallace et
    # al.'s intensity scale moves outside its calibration, 1.03 to 1.13.
    for comparison in validation.tdse_yields("Ar", "one_electron"):
        if "gauss30" in comparison.reference:
            continue  # not fitted: test_the_candidate_on_the_30_fs_pulse
        below = comparison.expected < 0.75
        mean = float(np.exp(np.mean(np.log(comparison.ratio[below]))))
        assert mean == pytest.approx(1.0, abs=0.1), comparison.reference
    scales = [
        validation.intensity_scale_test(gas, "default").scale
        for gas in ("Ar", "Kr", "Xe")
    ]
    # With krypton's 1.35 electrons (TC-846, TC-848) and xenon's 1.5 (TC-850).
    assert scales == pytest.approx([1.033, 1.023, 1.026], abs=0.01)
    candidate = ionisation_rate_model("Ar", 800e-9)
    assert candidate.identity.status == "validated"
    assert candidate.coulomb == "frequency"
    # In a static field it is the default.
    static = ionisation_rate_model("Ar", math.inf, "default_static")
    default = ionisation_rate_model("Ar", math.inf)
    assert static.instantaneous_rate(5e10) == pytest.approx(
        default.instantaneous_rate(5e10), rel=1e-12, abs=0.0
    )


def test_hydrogen():
    # The one-electron calculations of hydrogen (TC-833 to TC-835): below gamma = 1
    # every set agrees with the static factor within 1.05 to 1.32 and with the
    # candidate frequency factor within 0.92 to 1.03 (geometric means, yields below
    # 0.75); Popruzhenko et al.'s factor is low by 2.
    model = validation.hydrogen_model(800e-9)
    assert model.static_constants == (6.0,)
    by_form = {
        form: validation.hydrogen_yields(form) + validation.hydrogen_spectra(form)
        for form in ("static", "frequency", "popruzhenko")
    }
    assert len(by_form["static"]) == len(validation.HYDROGEN_YIELDS) + len(
        validation.HYDROGEN_SPECTRA
    )
    for k, comparison in enumerate(by_form["static"]):
        if math.isnan(comparison.wavelength):
            continue
        gamma = validation.hydrogen_model(comparison.wavelength).keldysh(
            validation.amplitude(comparison.x)
        )
        tunnelling = (gamma < 1) & (comparison.expected < 0.75)
        if not tunnelling.any():
            continue
        means = {
            form: float(np.exp(np.mean(np.log(found[k].ratio[tunnelling]))))
            for form, found in by_form.items()
        }
        assert 1.0 < means["static"] < 1.35, comparison.reference
        assert 0.9 < means["frequency"] < 1.05, comparison.reference
        assert means["popruzhenko"] < 0.6, comparison.reference
    spectrum = by_form["static"][-1]
    assert spectrum.reference == "bauer_2018/fig8_tdse_H_spectrum_F0063"
    assert np.all(np.diff(spectrum.x) > 0)
    assert spectrum.x.max() == pytest.approx(0.25, abs=0.003)


def test_the_pulse_in_cycles():
    times, envelope, field = validation.cycle_pulse(1e18, 800e-9, 2, 10)
    omega = 2 * math.pi * 299792458.0 / 800e-9
    potential = -(envelope / omega) * np.sin(omega * times)
    # The field is -dA/dt (to the grid's accuracy) and has no static part.
    assert np.max(np.abs(-np.gradient(potential, times) - field)) < 1e-3 * np.max(
        envelope
    )
    assert abs(np.trapezoid(field, times)) < 1e-12 * np.max(envelope) * times[-1]
    sine = validation.cycle_pulse(1e18, 800e-9, 3)[1]
    assert sine[0] == pytest.approx(0.0, abs=1e-12 * sine.max())
    assert sine[len(sine) // 2] == pytest.approx(sine.max())


def test_the_one_electron_calculation_in_wallaces_pulse():
    # Wallace et al.'s 5.5 fs pulse (TC-836): the static factor is 1.06 to 1.45 times
    # the calculation at both phases, the candidate frequency factor 0.86 to 1.04; put
    # through Wallace et al.'s test, the calculation itself matches their measured
    # yield at an intensity scale of 1.002 (calibration 1.3 %), and an atom's factor
    # of 0.8 on its rate moves it only to 1.016.
    for name, low, high in (
        ("one_electron_static", 1.05, 1.46),
        ("one_electron", 0.85, 1.05),
    ):
        wallace = [
            c for c in validation.tdse_yields("Ar", name) if "wallace" in c.reference
        ]
        assert len(wallace) == 2
        for comparison in wallace:
            assert np.all((comparison.ratio > low) & (comparison.ratio < high))
    exact = validation.tdse_intensity_scale_test("Ar")
    assert exact.scale == pytest.approx(1.002, abs=0.003)
    assert exact.scale - 1 < exact.calibration
    assert validation.tdse_intensity_scale_test("Ar", 0.8).scale == pytest.approx(
        1.016, abs=0.003
    )
    with pytest.raises(DataError, match="no one-electron calculation"):
        validation.tdse_intensity_scale_test("Ne")


def test_argon_at_390_nm_with_the_candidate_factor():
    # Van der Hart's all-electron rates of argon at 390 nm (TC-842): the static factor
    # puts one electron 1.7 to 3.1 times above them below gamma = 2.2, the candidate
    # frequency factor 0.85 to 1.37, about 1.15 on geometric average: the atom at
    # about 0.87 of one electron.
    logs = []
    for name, low, high in (
        ("one_electron_static", 1.6, 3.1),
        ("one_electron", 0.5, 1.4),
    ):
        for comparison in validation.rmatrix_rates("Ar", name):
            gamma = ionisation_rate_model("Ar", 390e-9).keldysh(
                validation.amplitude(comparison.x)
            )
            tunnelling = gamma < 2.2
            assert np.all(comparison.ratio[tunnelling] < high)
            if name == "one_electron_static":
                assert np.all(comparison.ratio[tunnelling & (gamma < 2.0)] > low)
            else:
                logs += list(np.log(comparison.ratio[tunnelling]))
    assert float(np.exp(np.mean(logs))) == pytest.approx(1.155, abs=0.02)


def test_hydrogen_from_two_codes():
    # tRecX's hydrogen and Zhao et al.'s read from their figure, the same 15 fs pulse at
    # 800 nm (TC-843): 0.96 on geometric average, 0.72 to 1.17 point by point through
    # the resonances between tRecX's coarser intensities.
    by_key = {d.key: d for d in datasets()}
    ours = by_key["travers_2026_trecx/ionised_gauss15_800nm_H"]
    theirs = by_key["zhao_2016/fig2b_tdse_H_800nm"]
    intensity, values = ours.abscissa_values(), ours.scaled_values()
    x, y = theirs.abscissa_values(), theirs.scaled_values()
    inside = (x >= intensity.min()) & (x <= intensity.max())
    interpolated = np.exp(
        np.interp(np.log(x[inside]), np.log(intensity), np.log(values))
    )
    ratio = y[inside] / interpolated
    assert float(np.exp(np.mean(np.log(ratio)))) == pytest.approx(0.96, abs=0.03)
    assert np.all((ratio > 0.7) & (ratio < 1.2))


def test_the_candidate_on_the_30_fs_pulse():
    # The 30 fs calculations at 400 and 1030 nm (TC-844), not fitted: the candidate
    # frequency factor is 0.87 and 0.95 times them on geometric average (yields below
    # 0.75), the static factor 2.0 and 1.18.
    expected = {
        ("one_electron_static", 400): 2.01,
        ("one_electron_static", 1030): 1.175,
        ("one_electron", 400): 0.87,
        ("one_electron", 1030): 0.95,
    }
    for name in ("one_electron_static", "one_electron"):
        for comparison in validation.tdse_yields("Ar", name):
            if "gauss30" not in comparison.reference:
                continue
            below = comparison.expected < 0.75
            mean = float(np.exp(np.mean(np.log(comparison.ratio[below]))))
            key = (name, round(comparison.wavelength * 1e9))
            assert mean == pytest.approx(expected[key], abs=0.03), key


def test_krypton_and_xenon_one_electron():
    # Kr 4p and Xe 5p in A2's 6 fs pulse (TC-845), yields below 0.75: the static factor
    # is 1.74 and 1.61 times the calculation on geometric average, the candidate
    # frequency factor 1.16 and 0.94. Krypton's 1.16 is its static normalisation: the
    # one-electron model is 1.08 to 1.17 times tRecX's own static rates of the same
    # potential, xenon's 0.84 to 1.04.
    expected = {
        ("one_electron_static", "Kr"): 1.74,
        ("one_electron_static", "Xe"): 1.61,
        ("one_electron", "Kr"): 1.16,
        ("one_electron", "Xe"): 0.94,
    }
    for (name, gas), value in expected.items():
        (comparison,) = [
            c for c in validation.tdse_yields(gas, name) if "gauss6" in c.reference
        ]
        below = comparison.expected < 0.75
        mean = float(np.exp(np.mean(np.log(comparison.ratio[below]))))
        assert mean == pytest.approx(value, abs=0.03), (name, gas)
    for gas, low, high in (("Kr", 1.07, 1.18), ("Xe", 0.83, 1.05)):
        (own,) = [c for c in validation.static_rates(gas) if "travers" in c.reference]
        assert np.all((own.ratio > low) & (own.ratio < high))


def test_krypton_in_wallaces_pulse():
    # Kr 4p, m = 0, in Wallace et al.'s pulse (TC-846): the exact one-electron
    # calculation through their test needs the intensities scaled by 1.045, beyond the
    # calibration of 1.5 %, where argon's matched at 1.002: the krypton atom ionises
    # more than one m = 0 electron. A factor of 1.47 on the calculation's rate brings
    # it to 0.989. (Santra et al.'s rate, 1.34 times tRecX's one-electron rate at
    # 0.10 au, is an independent-electron sum, not an all-electron one; TC-851.)
    exact = validation.tdse_intensity_scale_test("Kr")
    assert exact.scale == pytest.approx(1.045, abs=0.003)
    assert exact.scale - 1 > exact.calibration
    assert validation.tdse_intensity_scale_test("Kr", 1.47).scale == pytest.approx(
        0.989, abs=0.003
    )
    for name, low, high in (
        ("one_electron_static", 1.25, 2.35),
        ("one_electron", 1.05, 1.6),
    ):
        wallace = [
            c for c in validation.tdse_yields("Kr", name) if "wallace" in c.reference
        ]
        assert len(wallace) == 2
        for comparison in wallace:
            below = comparison.expected < 0.75
            assert np.all(
                (comparison.ratio[below] > low) & (comparison.ratio[below] < high)
            )


def test_xenon_in_wallaces_pulse():
    # Xe 5p, m = 0, in Wallace et al.'s pulse (TC-849): the exact one-electron
    # calculation through their test needs the intensities scaled by 1.019, within the
    # calibration of 2.5 %; a factor of 1.09 on its rate brings it to 1. The default,
    # with xenon's provisional one electron, is at 1.133: its fine structure lowers the
    # rate by about 0.74, and it would need about 1.5 electrons for 1.026.
    exact = validation.tdse_intensity_scale_test("Xe")
    assert exact.scale == pytest.approx(1.019, abs=0.003)
    assert exact.scale - 1 < exact.calibration
    assert validation.tdse_intensity_scale_test("Xe", 1.09).scale == pytest.approx(
        1.0, abs=0.003
    )
    # The frequency-dependent factor holds in a second pulse for xenon: 0.80 to 1.18
    # (geometric means 0.96, 1.00), the static 1.24 to 1.73.
    for name, low, high in (
        ("one_electron_static", 1.2, 1.8),
        ("one_electron", 0.78, 1.2),
    ):
        wallace = [
            c for c in validation.tdse_yields("Xe", name) if "wallace" in c.reference
        ]
        assert len(wallace) == 2
        for comparison in wallace:
            below = comparison.expected < 0.75
            assert np.all(
                (comparison.ratio[below] > low) & (comparison.ratio[below] < high)
            )


def test_the_argon_atom_against_one_electron():
    # haCC's all-electron argon (A2's 6 fs pulse, norm-loss at +1000 au) over one m = 0
    # electron at the same ionisation energy (TC-847): 0.76, 0.80 and 0.87 at 4.0e13,
    # 9.6e13 and 2.3e14 W/cm^2 with 4 channels, 0.81 with 3; the argon atom ionises
    # like about 0.8 of one electron, not like two or six.
    # With mMax 1, 0.86 (TC-861); in a 20 fs pulse, 0.84 against the same set-up's 0.81
    # at 6 fs (the lMax 30 pair: 0.81 against 0.75).
    c4, c3, c3m1, c3_20fs = validation.hacc_ratios()
    assert c4.reference.endswith("_c4")
    assert c4.ratio == pytest.approx([0.760, 0.797, 0.871], abs=0.005)
    assert c3.ratio == pytest.approx([0.809], abs=0.005)
    assert c3m1.ratio == pytest.approx([0.857], abs=0.005)
    assert c3_20fs.reference.endswith("gauss20_Ar_c3")
    assert c3_20fs.ratio == pytest.approx([0.836], abs=0.005)
    assert c4.ionization_energy == pytest.approx(0.5631718)
    assert c4.one_electron[1] < c4.atom[1] / 0.75


def test_the_weak_field_limit_on_the_formulas_scale():
    # Tolstikhina et al. 2014, Table I (TC-852): the atom's weak-field rate over the
    # model's one electron, (g/g_1)^2, g_1 = sqrt(2(2l+1)/kappa) 2^(-1/kappa) C_l.
    # Hydrogen fixes the normalisation: C_l = 2, kappa = 1 gives their g = sqrt 2.
    assert math.sqrt(2 * 1 / 1) * 2 ** (-1 / 1) * 2 == pytest.approx(math.sqrt(2))
    factors = {f.gas.value: f.factor for f in validation.weak_field_factors()}
    assert factors == pytest.approx(
        {"He": 1.665, "Ne": 0.525, "Ar": 1.723, "Kr": 1.784, "Xe": 1.824}, abs=0.005
    )
    # The factor 2 of the doubly occupied orbital, less the excess of Tong, Zhao &
    # Lin's C_l^2 over Hartree-Fock's, for all but neon.
    assert all(1.6 < factors[gas] < 2 for gas in ("He", "Ar", "Kr", "Xe"))


def test_the_weak_field_limit_needs_a_structure_coefficient(monkeypatch):
    original = validation.ionisation_rate_model

    def without(gas, wavelength, model="default"):
        return dataclasses.replace(
            original(gas, wavelength, model),
            structure="hartree",
            structure_coefficient=None,
        )

    monkeypatch.setattr(validation, "ionisation_rate_model", without)
    with pytest.raises(DataError, match="no structure coefficient"):
        validation.weak_field_factors()


def test_uiterwaals_saturation_intensities():
    # Uiterwaal et al. 1998, 248.6 nm, 500 fs (TC-853): the rate's integral over the
    # pulse at the measured saturation intensity. The frequency-dependent factor
    # removes most of the static factor's excess; helium and argon stay 3 to 4 high.
    exponents = {
        (gas, name): validation.saturation_test(gas, name).exponent
        for gas in ("He", "Ar", "Kr")
        for name in ("default", "default_static")
    }
    assert exponents == pytest.approx(
        {
            ("He", "default"): 3.38,
            ("Ar", "default"): 4.29,
            ("Kr", "default"): 0.55,
            ("He", "default_static"): 5.65,
            ("Ar", "default_static"): 31.6,
            ("Kr", "default_static"): 13.0,
        },
        rel=0.02,
        abs=0.0,
    )
    argon = validation.saturation_test("Ar")
    assert argon.keldysh == pytest.approx(3.70, abs=0.01)
    assert argon.scale == pytest.approx(0.578, abs=0.005)
    assert argon.sigma / argon.measured == pytest.approx(0.15)
    with pytest.raises(DataError, match="no measured saturation intensity"):
        validation.saturation_test("Xe")


@pytest.mark.filterwarnings("ignore::optmat.exceptions.ModelRangeWarning")
def test_argon_three_photon_cross_section_at_193_nm():
    # McKenna & van der Hart's R-matrix Floquet sigma_3 of argon at 6.4241 eV (193 nm),
    # read from their Fig. 8 (TC-856), against the default's effective sigma_3 = W /
    # Phi^3 at their 1e13 W/cm2 (gamma 15): the channel sum gives 0.74 of it, the
    # default's integral over the channels 1.5, and 4.3 and 12 times at 1e12 and 1e11.
    (dataset,) = [d for d in datasets() if d.key == "mckenna_2004/fig8_sigma3_Ar"]
    energy, sigma = dataset.abscissa_values(), dataset.scaled_values()
    order = np.argsort(energy)
    reference = float(np.exp(np.interp(6.4241, energy[order], np.log(sigma[order]))))
    assert reference == pytest.approx(6.43e-83, rel=0.01, abs=0.0)
    photon = 6.4241 * 1.602176634e-19  # J
    ratios = {}
    with pytest.warns(ModelRangeWarning, match="below the shortest wavelength"):
        argon = ionisation_rate_model("Ar", 193e-9)
    for channel_sum in ("integral", "sum"):
        model = dataclasses.replace(argon, channel_sum=channel_sum)
        for intensity in (1e11, 1e13):  # W/cm2
            rate = float(
                model.cycle_averaged_rate(validation.amplitude(intensity * 1e4))
            )
            ratios[channel_sum, intensity] = (
                rate / (intensity / photon) ** 3 / reference
            )
    assert ratios == pytest.approx(
        {
            ("integral", 1e11): 11.6,
            ("integral", 1e13): 1.50,
            ("sum", 1e11): 1.10,
            ("sum", 1e13): 0.742,
        },
        rel=0.02,
        abs=0.0,
    )


def test_published_many_electron_calculations_over_one_electron():
    # Published time-dependent CASSCF (and Pathak et al.'s own single-electron curve)
    # over one m = 0 electron in the same sin^2 pulse, both outside 20 au at the end
    # of the pulse (TC-858 to TC-860).
    ratios = {r.reference: r.ratio for r in validation.published_atom_ratios()}
    expected = {
        "pathak_2020/fig2_tdcasscf_Ar": [1.102, 0.947, 0.868],
        "pathak_2020/fig2_sae_Ar": [0.892, 0.869, 0.884],
        "pathak_2022/fig2_tdcasscf_Kr": [0.905],
        "sato_2018/fig3a_casscf_Ar_800nm": [0.797],
        "sato_2016/fig6_casscf_He_800nm": [1.707, 1.604, 1.531],
    }
    assert set(ratios) == set(expected)
    for key, values in expected.items():
        assert ratios[key] == pytest.approx(values, abs=0.005), key


def test_helium_with_two_electrons_at_390_nm():
    # Part K (TC-1050): the two-electron helium yield at 390 nm, 6 fs, over one electron
    # at the same ionisation energy (norm-loss 1000 au after the pulse, BOX 40): 1.71
    # and 1.79; with |m| <= 2, 1.72. haCC's helium over one electron 1.30 and 1.41,
    # 0.76 and 0.79 of the two-electron yield. Parker et al.'s rate integrated over the
    # pulse: 7.00e-3 and 8.49e-4, within 1 to 5 % of the two-electron yield.
    ratios = {
        r.reference.split("/")[1]: r for r in validation.helium_two_electron_ratios()
    }
    two = ratios["twoelectron_normloss_gauss6_390nm_He"]
    assert two.ratio == pytest.approx([1.709, 1.794], abs=0.005)
    assert ratios["twoelectron_normloss_gauss6_390nm_He_m5"].ratio == pytest.approx(
        [1.722], abs=0.005
    )
    hacc = ratios["hacc_normloss_gauss6_390nm_He"]
    assert hacc.ratio == pytest.approx([1.305, 1.410], abs=0.005)
    assert hacc.atom / two.atom == pytest.approx([0.763, 0.786], abs=0.005)
    parker = [validation.parker_pulse_yield(i, 6e-15, 390) for i in two.intensity]
    assert parker == pytest.approx([7.000e-3, 8.491e-4], rel=0.002, abs=0.0)
    assert two.atom / np.array(parker) == pytest.approx([1.047, 1.016], abs=0.005)

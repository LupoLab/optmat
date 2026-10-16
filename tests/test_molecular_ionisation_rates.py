# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The molecular rates' literature forms: their formulas, and parity with solprop."""

import json
import math
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from optmat.constants import AU_ANGULAR_FREQUENCY, AU_FIELD, E_CHARGE
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.gas import _tunnelling
from optmat.gas.ionisation import ionisation_parameters
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.ionisation_validation import measured_yields, pulse_yield
from optmat.gas.molecular_ionisation_rates import (
    isotropic_structure_factors,
    moadk_rate_au,
    molecular_rate_model,
    molecular_rate_models,
)
from optmat.gas.molecular_ionisation_validation import (
    fitted_constants,
    orientation_average,
    orientation_averaged_rates,
    vibrational_factors,
)

FIXTURE = json.loads(
    (Path(__file__).parent / "data" / "solprop_molecular_rates.json").read_text(
        encoding="utf-8"
    )
)
FIELDS = np.array(FIXTURE["fields_v_per_m"])
# solprop's names of the tabulations, and the models they are here.
MOADK = {"zhao2010": "moadk_zhao_2010", "tong2002": "moadk_tong_2002"}


def _with_solprop_energy(model, gas):
    """Evaluate a model with solprop's ionisation energy: the formulas are compared."""
    energy = FIXTURE["molecules"][gas]["potential_j"]
    return model.with_ionization_energy(energy, "solprop's")


def _rates(method, fields):
    with warnings.catch_warnings():
        # The fixture's fields run to 0.117 au, above the barrier-suppression fields.
        warnings.simplefilter("ignore", ModelRangeWarning)
        return method(fields)


def _nonzero(expected):
    # solprop returns exactly 0 below its threshold field and clips at its ceiling of
    # one per atomic unit of time; neither is part of the formula.
    ceiling = 0.999 * AU_ANGULAR_FREQUENCY
    return (np.asarray(expected) > 0) & (np.asarray(expected) < ceiling)


# The formulas.


def test_the_isotropic_average_does_not_depend_on_the_orbitals_m():
    factors = isotropic_structure_factors({2: 0.52, 4: 0.03})
    assert sorted(factors) == list(range(-4, 5))
    assert factors[1] == factors[-1]
    # <B^2(m')> = sum_l C_l^2 (l + |m'|)! / (2 (l - |m'|)!): 1/2 for m' = 0, and for
    # m' = 2, 4!/2 = 12 (l = 2) and 6!/(2 2!) = 180 (l = 4).
    assert factors[0] == pytest.approx((0.52**2 + 0.03**2) / 2, rel=1e-12, abs=0.0)
    assert factors[2] == pytest.approx(0.52**2 * 12 + 0.03**2 * 180, rel=1e-12, abs=0.0)


def test_the_molecular_rate_is_the_private_validators():
    factors = isotropic_structure_factors({0: 2.68, 2: 1.10, 4: 0.06})
    for field in (0.03, 0.08):
        expected = _tunnelling.static_rate_au(field, 0.5726, factors)
        assert float(moadk_rate_au(np.array(field), 0.5726, factors)) == pytest.approx(
            expected, rel=1e-12, abs=0.0
        )


def test_an_atoms_single_partial_wave_gives_the_atomic_tunnelling_rate():
    # One partial wave l with C_l, averaged over orientations, is the atomic rate
    # averaged over the orbital's m: for argon's 3p, (w_0 + 2 w_1) / 3.
    argon = ionisation_rate_model("Ar", math.inf, "one_electron_static")
    energy = argon.channels[0].ionization_energy / 4.3597447222071e-18
    field = np.array([0.04, 0.07])
    coefficient = argon.structure_coefficient
    assert coefficient is not None
    factors = isotropic_structure_factors({1: coefficient})
    kappa = math.sqrt(2 * energy)
    # m = 0 alone, without the static correction: B^2(0) = 3 C^2 / 2 aligned.
    aligned = moadk_rate_au(field, energy, {0: 1.5 * coefficient**2})
    averaged = moadk_rate_au(field, energy, factors)
    ratio = averaged / aligned
    # w_1 / w_0 = F / (2 kappa^3) for a p orbital (f(1, 1) = f(1, 0) = 3).
    expected = (1 + 2 * field / (2 * kappa**3)) / 3
    assert ratio == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_the_cycle_average_of_the_static_forms_is_the_tunnelling_factor():
    for name in ("moadk_zhao_2010", "adk_hydrogenic"):
        model = molecular_rate_model("N2", 800e-9, name)
        field = np.array([2e10, 3e10])
        f0 = (2 * model.ionization_energy_au) ** 1.5
        factor = np.sqrt(3 * field / AU_FIELD / (math.pi * f0))
        assert model.cycle_averaged_rate(field) == pytest.approx(
            factor * model.instantaneous_rate(field), rel=1e-12, abs=0.0
        )


def test_the_effective_charge_form_tends_to_its_tunnelling_limit():
    # At 30 um the Keldysh parameter is 0.03: the rate is the tunnelling rate of an
    # s electron with the effective charge, to gamma^2.
    model = molecular_rate_model("N2", 30e-6, "ppt_talebpour_ilkov")
    field = np.array([3e10])
    energy, charge = model.ionization_energy_au, model.charge
    kappa = math.sqrt(2 * energy)
    f_au = field / AU_FIELD
    static = (
        energy
        * model.coefficient_squared
        * (2 * kappa**3 / f_au) ** (2 * charge / kappa - 1)
        * np.exp(-2 * kappa**3 / (3 * f_au))
    )
    assert float(model.keldysh(field)[0]) < 0.04
    assert model.instantaneous_rate(field) / AU_ANGULAR_FREQUENCY == pytest.approx(
        static, rel=3e-3, abs=0.0
    )


def test_zero_field_and_the_range():
    model = molecular_rate_model("O2", 800e-9, "moadk_zhao_2010")
    assert float(model.instantaneous_rate(0.0)) == 0.0
    assert model.instantaneous_rate(np.array([-2e10, 2e10]))[0] == pytest.approx(
        float(model.instantaneous_rate(2e10))
    )
    with pytest.warns(ModelRangeWarning, match="barrier-suppression"):
        model.instantaneous_rate(1.1 * model.barrier_suppression_field)


# Parity with solprop (the fixture of tools/solprop_molecular_rates.py).


@pytest.mark.parametrize("gas", ["H2", "D2", "N2", "O2"])
@pytest.mark.parametrize("tabulation", sorted(MOADK))
def test_parity_of_the_molecular_tunnelling_rate(gas, tabulation):
    entry = FIXTURE["molecules"][gas]["moadk"][tabulation]
    assert entry["occupancy"] == 1
    model = _with_solprop_energy(
        molecular_rate_model(gas, 800e-9, MOADK[tabulation]), gas
    )
    expected = np.array(entry["instantaneous"])
    keep = _nonzero(expected)
    # The same formula evaluated twice: rounding only.
    assert _rates(model.instantaneous_rate, FIELDS)[keep] == pytest.approx(
        expected[keep], rel=1e-10, abs=0.0
    )


@pytest.mark.parametrize("gas", ["H2", "D2", "CH4", "SF6", "N2O"])
def test_parity_of_the_hydrogenic_rate(gas):
    entry = FIXTURE["molecules"][gas]["adk"]
    assert entry["potential_only"]
    model = replace(
        _with_solprop_energy(molecular_rate_model(gas, 800e-9, "adk_hydrogenic"), gas),
        electrons=float(entry["occupancy"]),
    )
    for name, method in (
        ("instantaneous", model.instantaneous_rate),
        ("cycle_averaged", model.cycle_averaged_rate),
    ):
        expected = np.array(entry[name])
        keep = _nonzero(expected)
        assert _rates(method, FIELDS)[keep] == pytest.approx(
            expected[keep], rel=1e-10, abs=0.0
        )


@pytest.mark.parametrize("gas", ["N2", "O2"])
@pytest.mark.parametrize("wavelength", ["800", "400"])
def test_parity_of_the_effective_charge_rate(gas, wavelength):
    entry = FIXTURE["molecules"][gas]["ppt"][wavelength]
    settings = dict(entry["settings"])
    # solprop's defaults for a molecule: no barrier-suppression factor, and no
    # polarizability for its Stark shift or dipole factor to act with.
    assert settings["bsi"] == "none"
    assert settings["sum_integral"]
    model = replace(
        _with_solprop_energy(
            molecular_rate_model(gas, float(wavelength) * 1e-9, "ppt_talebpour_ilkov"),
            gas,
        ),
        electrons=float(entry["occupancy"]),
    )
    for name, method in (
        ("instantaneous", model.instantaneous_rate),
        ("cycle_averaged", model.cycle_averaged_rate),
    ):
        expected = np.array(entry[name])
        keep = _nonzero(expected)
        # The functions of the Keldysh parameter are evaluated by different
        # expressions in the two codes: 1e-9 of rounding in an exponent of 30.
        assert _rates(method, FIELDS)[keep] == pytest.approx(
            expected[keep], rel=1e-8, abs=0.0
        )


@pytest.mark.parametrize(
    ("gas", "difference", "change"),
    [
        # optmat's energy less solprop's, meV, and the rate with optmat's over that
        # with solprop's at 2e10 V/m (1e14 W/cm^2 is 2.7e10 V/m).
        ("H2", -0.3, 1.0007),
        ("D2", -0.5, 1.0011),
        ("N2", -0.5, 1.0011),
        ("O2", -0.1, 1.0001),
        ("CH4", -5.2, 1.0102),
        ("SF6", 30.0, 0.9382),
        ("N2O", -12.2, 1.0243),
    ],
)
def test_what_optmats_energies_change(gas, difference, change):
    # solprop's energies are rounded hartree values (SF6's 15.32 eV); optmat's are
    # the primary determinations. With them the rates move by 0.1 % (H2, N2, O2), 1
    # to 2 % (CH4, N2O) and 6 % (SF6), which is what parity cannot mean.
    model = molecular_rate_model(gas, 800e-9, "adk_hydrogenic")
    other = _with_solprop_energy(model, gas)
    shift = (model.ionization_energy - other.ionization_energy) / E_CHARGE * 1e3
    assert shift == pytest.approx(difference, abs=0.05)
    field = np.array(2e10)
    ratio = float(model.instantaneous_rate(field) / other.instantaneous_rate(field))
    assert ratio == pytest.approx(change, abs=1e-4)


# The models.


def test_the_models_and_their_molecules():
    models = molecular_rate_models()
    assert sorted(g.value for g in models["moadk_zhao_2010"]) == [
        "CO",
        "CO2",
        "D2",
        "H2",
        "N2",
        "O2",
    ]
    assert sorted(g.value for g in models["ppt_talebpour"]) == ["D2", "N2", "O2"]
    assert sorted(g.value for g in models["ppt_talebpour_ilkov"]) == ["N2", "O2"]
    assert "N2O" in {g.value for g in models["moadk_zhao_2017"]}
    for name, by_gas in models.items():
        for model in by_gas.values():
            # The forms of the literature are of one electron; the default counts the
            # shell's.
            assert (model.electrons == 1.0) == (name != "default")
            # The forms take the recommended (adiabatic) energy, the default the
            # vertical one (TC-679). For CO no adiabatic energy is curated, and the
            # vertical one stands for both (TC-683).
            expected = (
                {"vertical"}
                if name == "default" or str(model.gas) == "CO"
                else {"adiabatic", "not stated"}
            )
            assert model.ionization_energy_kind in expected
            assert model.citation
    nitrogen = models["ppt_talebpour"][
        next(g for g in models["ppt_talebpour"] if g.value == "N2")
    ]
    assert nitrogen.charge == 0.9
    assert nitrogen.sources["charge"] == "talebpour_1999/zeff_N2"


def test_the_charges_were_fitted_within_the_static_coulomb_factor():
    # The rate of the source over solprop's is (1 + gamma^2)^(n*): at a Keldysh
    # parameter of 1, 2^(n*), 1.80 for N2 with Z_eff = 0.9 (n* = 0.841).
    source = molecular_rate_model("N2", 800e-9, "ppt_talebpour")
    solprop = molecular_rate_model("N2", 800e-9, "ppt_talebpour_ilkov")
    amplitude = np.array([3.0e10])
    gamma = float(source.keldysh(amplitude)[0])
    n_star = source.charge / math.sqrt(2 * source.ionization_energy_au)
    ratio = source.cycle_averaged_rate(amplitude) / solprop.cycle_averaged_rate(
        amplitude
    )
    assert ratio == pytest.approx((1 + gamma**2) ** n_star, rel=1e-12, abs=0.0)
    assert n_star == pytest.approx(0.841, abs=1e-3)


def test_zhao_2017_against_zhao_2010():
    # The same group's later coefficients: +60 % in the rate of N2, +52 % for O2,
    # -7 % for H2, -22 % for CO2 (the ratio does not depend on the field to 1 %).
    field = np.array([2.5e10])
    for gas, change in (("N2", 1.60), ("O2", 1.52), ("H2", 0.93), ("CO2", 0.78)):
        new = molecular_rate_model(gas, 800e-9, "moadk_zhao_2017")
        old = molecular_rate_model(gas, 800e-9, "moadk_zhao_2010")
        ratio = float(
            (new.instantaneous_rate(field) / old.instantaneous_rate(field))[0]
        )
        assert ratio == pytest.approx(change, abs=0.02)


def test_errors():
    with pytest.raises(DataError, match="unknown molecular rate model"):
        molecular_rate_model("N2", 800e-9, "mo_ppt")
    with pytest.raises(DataError, match="not given for CH4"):
        molecular_rate_model("CH4", 800e-9, "moadk_zhao_2010")
    with pytest.raises(ValueError, match="finite wavelength"):
        molecular_rate_model("N2", math.inf, "ppt_talebpour").instantaneous_rate(2e10)


# The static forms against the orientation average of angle-resolved rates.


def test_the_orientation_average():
    angles = np.arange(0.0, 91.0, 15.0)
    cosine = np.cos(np.radians(angles))
    # 1 + cos^2 averages to 4/3 over the sphere; the interpolation of seven angles is
    # good to 1e-3.
    assert float(orientation_average(angles, 1 + cosine**2)) == pytest.approx(
        4 / 3, rel=1e-3, abs=0.0
    )
    # Two angles: (parallel + 2 perpendicular) / 3.
    two = orientation_average(np.array([0.0, 90.0]), np.array([[3.0], [1.5]]))
    assert two == pytest.approx([2.0], rel=1e-12, abs=0.0)
    # A polar molecule, 0 to 180 degrees: 2 + cos averages to 2, and 1 + cos^2 to
    # 4/3, both through the interpolation of thirteen angles in cos.
    polar = np.arange(0.0, 181.0, 15.0)
    cosine = np.cos(np.radians(polar))
    assert float(orientation_average(polar, 2 + cosine)) == pytest.approx(
        2, rel=1e-4, abs=0.0
    )
    assert float(orientation_average(polar, 1 + cosine**2)) == pytest.approx(
        4 / 3, rel=1e-3, abs=0.0
    )
    with pytest.raises(ValueError, match="span 0 to 90 or 0 to 180"):
        orientation_average(np.array([0.0, 45.0]), np.array([1.0, 2.0]))


def test_nitrogen_against_its_orientation_averaged_static_rates():
    zhao = orientation_averaged_rates("N2", "moadk_zhao_2010")
    field = zhao.field / AU_FIELD
    assert list(zhao.angles) == [0, 15, 30, 45, 60, 75, 90]
    # The rate along the axis is 2.6 times the average at 0.03 au and 1.2 times at
    # 0.14 au: the anisotropy falls as the barrier is suppressed.
    assert zhao.anisotropy[0] == pytest.approx(2.57, abs=0.01)
    assert zhao.anisotropy[-1] == pytest.approx(1.21, abs=0.01)
    # Without a correction at finite field the form is within 10 % at the lowest
    # fields and 2.8 times too high at the barrier-suppression field (0.08 au).
    assert zhao.ratio[field < 0.045] == pytest.approx(1.0, abs=0.11)
    assert zhao.ratio[np.isclose(field, 0.08)] == pytest.approx(2.78, abs=0.01)
    assert np.all(np.diff(zhao.ratio) > 0)
    tong = orientation_averaged_rates("N2", "moadk_tong_2002")
    # Tong et al.'s coefficients give 0.56 of Zhao et al.'s rate.
    assert tong.ratio / zhao.ratio == pytest.approx(0.556, abs=0.003)


def test_hydrogen_against_its_orientation_averaged_static_rates():
    zhao = orientation_averaged_rates("H2", "moadk_zhao_2010")
    assert list(zhao.angles) == [0, 90]
    # One electron of the model against both of the molecule, fixed nuclei: the form
    # is 1.9 to 2.2 times too high from 0.06 to 0.12 au.
    assert zhao.ratio == pytest.approx(2.0, abs=0.2)
    with pytest.raises(DataError, match="no angle-resolved"):
        orientation_averaged_rates("O2", "moadk_zhao_2010")


# optmat's rate of molecules (the model "default").


def test_the_default_models_constants():
    expected = {
        # gas: electrons, alpha, energy of the rate (eV)
        "H2": (2.0, 6.0, 16.06),
        "D2": (2.0, 6.0, 16.06),
        "N2": (2.0, 14.0, 15.58),
        "O2": (2.0, 6.1, 12.307),
        "CH4": (6.0, 18.0, 14.4),
        "CO2": (4.0, 14.0, 13.79),
        "N2O": (4.0, 14.0, 12.89),
        "CO": (2.0, 14.0, 14.01),
        "SF6": (6.0, 14.0, 15.69),
    }
    for gas, (electrons, alpha, energy) in expected.items():
        model = molecular_rate_model(gas, 800e-9)
        assert model.form == "molecular"
        assert model.coulomb == "static"
        assert model.averaging == "consistent"
        assert model.electrons == electrons
        assert model.tong_lin_alpha == alpha
        assert model.ionization_energy / E_CHARGE == pytest.approx(energy, abs=1e-9)
        assert model.ionization_energy_kind == "vertical"
        # The adiabatic energy is the lower, by up to 0.64 eV (H2); N2's vertical
        # energy is printed to 0.01 eV and is 0.7 meV below the adiabatic one.
        assert model.adiabatic_ionization_energy <= model.ionization_energy + 1e-3 * (
            E_CHARGE
        )
        assert model.uncertainty_factor >= 1.5
    assert molecular_rate_model("H2", 800e-9).vibrational_constant == 0.00165
    assert molecular_rate_model("N2", 800e-9).vibrational_constant == 0.0
    assert molecular_rate_model("N2O", 800e-9).sources["structure"] == (
        "zhao_2017/table2_N2O"
    )


def test_one_partial_wave_is_the_atoms_rate():
    # The molecular form with one s wave is the noble gases' rate: helium's default
    # (C_0 = 3.13, alpha = 7, N_eff = 1.47, no fine structure) at every wavelength,
    # both rates.
    # The molecular default keeps the static Coulomb factor (TC-848).
    helium = ionisation_rate_model("He", 800e-9, "default_static")
    energy = helium.channels[0].ionization_energy
    coefficient = helium.structure_coefficient
    assert coefficient is not None
    like = replace(
        molecular_rate_model("H2", 800e-9),
        ionization_energy=energy,
        structure_factors=isotropic_structure_factors({0: coefficient}),
        tong_lin_alpha=helium.static_constants[0],
        electrons=helium.electrons,
        vibrational_constant=0.0,
    )
    # Within both ranges: the molecule's ends at 1.5 times the barrier-suppression
    # field, the atom's at 2 (1.6e11 and 2.1e11 V/m), above which each is held.
    field = np.array([6e10, 1e11, 1.5e11])
    for wavelength in (800e-9, 400e-9, math.inf):
        atom = helium.at_wavelength(wavelength)
        molecule = like.at_wavelength(wavelength)
        # The same formula evaluated through two paths: rounding only.
        assert _rates(molecule.instantaneous_rate, field) == pytest.approx(
            _rates(atom.instantaneous_rate, field), rel=1e-9, abs=0.0
        )
        if not math.isinf(wavelength):
            assert _rates(molecule.cycle_averaged_rate, field) == pytest.approx(
                _rates(atom.cycle_averaged_rate, field), rel=1e-9, abs=0.0
            )


def test_the_default_in_a_static_field():
    # N_eff V S times the orientation-averaged tunnelling rate of one electron.
    model = molecular_rate_model("H2", math.inf)
    field = np.array([0.05, 0.08])
    energy = model.ionization_energy_au
    bare = moadk_rate_au(field, energy, model.structure_factors)
    correction = np.exp(
        -model.tong_lin_alpha * field / (energy * (2 * energy) ** 1.5)
        + model.vibrational_constant / field**2
    )
    expected = model.electrons * correction * bare * AU_ANGULAR_FREQUENCY
    assert model.instantaneous_rate(field * AU_FIELD) == pytest.approx(
        expected, rel=1e-12, abs=0.0
    )


@pytest.mark.parametrize("gas", ["N2", "O2", "H2"])
def test_the_two_rates_of_the_default_give_the_same_yield(gas):
    # The instantaneous rate on the carrier and the cycle-averaged rate on the
    # envelope of a 42 fs pulse, whose envelope is constant over a cycle to 1 %: the
    # yields agree as the rates' cycle averages do (3 %).
    model = molecular_rate_model(gas, 800e-9)
    peak = (
        0.5
        * 8.8541878128e-12
        * 299792458
        * (0.6 * model.barrier_suppression_field) ** 2
    )
    on_envelope = pulse_yield(model, peak, 42e-15)
    on_carrier = pulse_yield(model, peak, 42e-15, carrier=True)
    assert on_carrier == pytest.approx(on_envelope, rel=0.03, abs=0.0)


def test_the_default_tends_to_its_static_rate():
    # At 30 um the Keldysh parameter is below 0.04: the cycle-averaged rate is the
    # average over a cycle of the static rate (a trapezoid over 4000 phases), to the
    # frequency's own correction, of order gamma^2 (0.2 %).
    static = molecular_rate_model("N2", math.inf)
    slow = molecular_rate_model("N2", 30e-6)
    amplitude = 3e10
    phase = np.linspace(0, math.pi / 2, 4001)
    average = np.trapezoid(
        static.instantaneous_rate(amplitude * np.cos(phase)), phase
    ) / (math.pi / 2)
    assert float(slow.cycle_averaged_rate(amplitude)) == pytest.approx(
        average, rel=0.005, abs=0.0
    )


def test_the_fitted_constants_are_the_curated_ones():
    fitted = fitted_constants()
    # The curated alphas are the least-squares values rounded to 0.1.
    assert fitted.alpha["N2"] == pytest.approx(14.0, abs=0.05)
    assert fitted.alpha["H2"] == pytest.approx(6.0, abs=0.05)
    # With the count free the references ask for 1.76 (N2) and 2.12 (H2) electrons:
    # the weak-field limit of 2, within 12 %.
    assert fitted.free_electrons["N2"] == pytest.approx(2.0, rel=0.13, abs=0.0)
    assert fitted.free_electrons["H2"] == pytest.approx(2.0, rel=0.07, abs=0.0)
    assert fitted.vibrational_constant == pytest.approx(0.00165, abs=5e-6)
    # CO, which takes N2's 14 unfitted, would ask for more: the offset of 1.3 is
    # read as a slope when the count is held at 2 (TC-683).
    assert fitted.alpha["CO"] > fitted.alpha["N2"]
    # O2 and CH4 against Jagau's CCSD rates (TC-698).
    assert fitted.alpha["O2"] == pytest.approx(6.1, abs=0.05)
    assert fitted.alpha["CH4"] == pytest.approx(18.0, abs=0.05)


def test_the_default_against_the_static_references():
    nitrogen = orientation_averaged_rates("N2")
    field = nitrogen.field / AU_FIELD
    # 0.81 at 0.14 au (1.7 times the barrier-suppression field) to 1.13 at 0.08 au:
    # the references end inside the correction's range (TC-828).
    assert np.all(nitrogen.field <= molecular_rate_model("N2", math.inf).maximum_field)
    assert nitrogen.ratio == pytest.approx(1.0, abs=0.2)
    assert nitrogen.ratio[field <= 0.1201] == pytest.approx(1.0, abs=0.14)
    hydrogen = orientation_averaged_rates("H2")
    assert hydrogen.ratio == pytest.approx(1.0, abs=0.035)
    vibration = vibrational_factors()
    assert vibration.ratio == pytest.approx(1.0, abs=0.035)


def test_carbon_monoxide_tests_nitrogens_alpha():
    # CO takes N2's alpha without a fit: against the average over its thirteen angles
    # (0 to 180 degrees; majety_2015 Table A2) the default is 1.26 to 1.32 times the
    # reference from 0.03 to 0.09 au (0.45 to 1.36 times its barrier-suppression
    # field), a constant offset in the magnitude. At 0.025 au it is 2.25, where the
    # reference's anisotropy jumps from 1.4 to 2.7 (TC-683).
    carbon = orientation_averaged_rates("CO")
    field = carbon.field / AU_FIELD
    assert len(carbon.angles) == 13
    assert carbon.angles[-1] == 180
    inside = field >= 0.0299
    assert carbon.ratio[inside] == pytest.approx(1.29, abs=0.035)
    assert carbon.ratio[~inside] == pytest.approx([2.25], abs=0.01)
    # The forms without a correction rise against it as N2's do.
    zhao = orientation_averaged_rates("CO", "moadk_zhao_2010")
    assert np.all(np.diff(zhao.ratio[inside]) > 0)
    assert zhao.ratio[-1] > 6


def test_jagaus_references():
    # Jagau 2018's correlated rates, averaged over orientations (TC-698). H2 (full CI)
    # and CO (CCSD) cross-check Majety and Scrinzi's: the default is within 4 % of H2's
    # from 0.08 au (1.30 at 0.06 au, where the rate is 2e-5, near the method's floor)
    # and 1.14 to 1.54 times CO's, Majety's being 7 to 17 % higher at 0.06 au. O2 and
    # CH4 are the references their alphas were fitted to.
    hydrogen = orientation_averaged_rates("H2", reference="jagau_2018")
    assert hydrogen.ratio[1:] == pytest.approx(1.0, abs=0.045)
    carbon = orientation_averaged_rates("CO", reference="jagau_2018")
    assert list(carbon.field / AU_FIELD) == pytest.approx([0.06, 0.09])
    assert carbon.ratio == pytest.approx([1.54, 1.14], abs=0.01)
    majety = orientation_averaged_rates("CO")
    at = int(np.argmin(abs(majety.field / AU_FIELD - 0.06)))
    assert majety.expected[at] / carbon.expected[0] == pytest.approx(1.18, abs=0.02)
    oxygen = orientation_averaged_rates("O2", reference="jagau_2018")
    assert len(oxygen.angles) == 9
    assert list(oxygen.field / AU_FIELD) == pytest.approx([0.06, 0.10])
    assert oxygen.ratio == pytest.approx(1.0, abs=0.1)
    methane = orientation_averaged_rates("CH4", reference="jagau_2018")
    assert list(methane.angles) == pytest.approx([0.0, 54.7356, 180.0], abs=1e-4)
    # 0.04 au is the method's low-field limit, excluded in the record.
    assert list(methane.field / AU_FIELD) == pytest.approx([0.06, 0.08, 0.10, 0.12])
    assert methane.ratio == pytest.approx(1.0, abs=0.28)
    # Lebedev's weights, 0.3 A + 0.4 C + 0.3 B, of the three orientations.
    assert methane.expected == pytest.approx(
        0.3 * methane.resolved[0]
        + 0.4 * methane.resolved[1]
        + 0.3 * methane.resolved[2]
    )
    with pytest.raises(DataError, match="no angle-resolved static rates of jagau_2018"):
        orientation_averaged_rates("N2", reference="jagau_2018")


def test_the_default_against_wahlstrands_yields():
    # The measurement of TC-660: the model is above it for the molecules as for the
    # atoms, N2 by 5 to 9 and O2 by 5.0 to 6.5 (Ar: 1.7 to 3.8), O2's since its alpha
    # was fitted to Jagau's CCSD rates (TC-698).
    nitrogen = measured_yields("N2")
    assert np.all((nitrogen.ratio > 5.0) & (nitrogen.ratio < 9.0))
    oxygen = measured_yields("O2")
    assert np.all((oxygen.ratio > 5.0) & (oxygen.ratio < 6.6))


def test_the_default_against_guos_molecule_to_atom_yields():
    # Guo et al. 1998 (800 nm, 30 fs; the text, not the figures): N2 behaves like Ar,
    # and 'the O2+ LP ion yield is roughly ten times below Xe+ in the unsaturated
    # range'. The default's N2/Ar is 0.84 to 1.4 from 2e13 to 1.5e14 W/cm2; its
    # O2/Xe is 0.035 to 0.055 from 2e13 to 7e13 W/cm2 (Xe's yield below 0.3), 2 to 3
    # times lower than stated; with N2's alpha it was 4.5 times (TC-683, TC-698).
    bounds = (("N2", "Ar", 0.8, 1.45), ("O2", "Xe", 0.035, 0.056))
    for molecule, atom, low, high in bounds:
        for intensity in (2e17, 3e17, 5e17, 7e17):
            ratio = pulse_yield(
                molecular_rate_model(molecule, 800e-9), intensity, 30e-15
            ) / pulse_yield(
                ionisation_rate_model(atom, 800e-9, "default_static"), intensity, 30e-15
            )
            assert low < ratio < high


def test_methanes_three_orbitals_are_combined():
    # Zhao et al. 2011: three degenerate 1t2 orbitals, complex C_lm at several m.
    structure = ionisation_parameters("CH4").structure["zhao_2011"]
    orbitals = {key[0] for key in structure.projections}
    assert len(orbitals) == 3
    # HOMO2's C_{1,+1} is printed 2.23i: the imaginary column is read.
    assert any(
        value == complex(0.0, 2.23)
        for (_, wave, m), value in structure.projections.items()
        if (wave, m) == (1, 1)
    )
    squares = structure.partial_wave_squares
    # Per orbital: 9.9458, 9.9458 and 9.9856 at l = 1 (the table's rounding).
    assert squares[1] == pytest.approx((2 * 9.9458 + 9.9856) / 3, rel=1e-9, abs=0.0)
    assert squares[2] == pytest.approx(1.28, rel=1e-9, abs=0.0)
    model = molecular_rate_model("CH4", 800e-9)
    assert model.electrons == 6.0
    assert sum(model.structure_factors.values()) > 0


def test_the_rate_falls_to_zero_with_the_field():
    # The vibrational factor's Gaussian is not followed beyond two standard
    # deviations, so the rate of H2 falls monotonically to zero with the field,
    # static and at 800 nm, both rates.
    field = np.geomspace(1e-4, 0.1, 200) * AU_FIELD
    for wavelength in (math.inf, 800e-9):
        model = molecular_rate_model("H2", wavelength)
        rates = [_rates(model.instantaneous_rate, field)]
        if not math.isinf(wavelength):
            rates.append(_rates(model.cycle_averaged_rate, field))
        for rate in rates:
            assert rate[0] < 1e-10
            assert np.all(np.diff(rate) >= 0)


def test_sulphur_hexafluorides_three_orbitals_are_combined():
    # Zhao et al. 2011, read on the page (TC-684): three degenerate 1t1g orbitals,
    # l = 4 to 10. Their sums of |C_lm|^2 over m agree per l to the table's rounding.
    structure = ionisation_parameters("SF6").structure["zhao_2011"]
    assert len({key[0] for key in structure.projections}) == 3
    squares = structure.partial_wave_squares
    assert sorted(squares) == [4, 6, 8, 10]
    assert squares[4] == pytest.approx(145.49, abs=0.01)
    assert squares[6] == pytest.approx(3.25, abs=0.01)
    # The l = 4 wave dominates: the default is 80 times N2's at 2e14 W/cm2, for
    # nearly the same ionisation energy, and 100 times the hydrogen-like form.
    field = math.sqrt(2 * 2e18 / (8.8541878128e-12 * 299792458))
    rate = float(molecular_rate_model("SF6", 800e-9).cycle_averaged_rate(field))
    nitrogen = molecular_rate_model("N2", 800e-9)
    hydrogenic = molecular_rate_model("SF6", 800e-9, "adk_hydrogenic")
    assert rate / float(nitrogen.cycle_averaged_rate(field)) == pytest.approx(
        83, rel=0.05, abs=0.0
    )
    assert rate / float(hydrogenic.cycle_averaged_rate(field)) == pytest.approx(
        104, rel=0.05, abs=0.0
    )


@pytest.mark.parametrize("gas", ["H2", "D2", "N2", "O2", "CO2", "N2O", "CH4", "SF6"])
def test_the_default_rises_to_its_maximum_field(gas):
    # With N2's alpha of 14 the rate stops rising at 1.7 times the
    # barrier-suppression field; the range the default is given in ends at 1.5
    # (TC-683). Both rates, static, at 800 nm and at 1.8 um.
    for wavelength in (math.inf, 800e-9, 1800e-9):
        model = molecular_rate_model(gas, wavelength)
        field = np.geomspace(0.3, 1.0, 200) * model.maximum_field
        for rate in (
            model.instantaneous_rate(field),
            model.cycle_averaged_rate(field),
        ):
            assert np.all(np.diff(rate) > 0)


def test_the_vibrational_factor_is_continuous():
    # At two standard deviations, x = sqrt(2c)/F = 2, the Gaussian's x^2/2 and the
    # tail's 2x - 2 meet with equal slopes.
    model = molecular_rate_model("H2", math.inf)
    switch = math.sqrt(2 * model.vibrational_constant) / 2
    fixed = replace(model, vibrational_constant=0.0)
    field = np.array([switch * (1 - 1e-6), switch * (1 + 1e-6)]) * AU_FIELD
    ratio = model.instantaneous_rate(field) / fixed.instantaneous_rate(field)
    assert ratio == pytest.approx(math.exp(2.0), rel=1e-4, abs=0.0)


@pytest.mark.parametrize(
    "gas", ["H2", "D2", "N2", "O2", "CO", "CO2", "N2O", "CH4", "SF6"]
)
def test_the_correction_range_lies_below_the_turn_over(gas):
    # Each molecule's range ends at the largest field its alpha was fitted to or
    # checked against, below the field where the rate would turn over (TC-828).
    model = molecular_rate_model(gas, 800e-9)
    assert 1.6 < model.correction_range < 2.0
    field = np.linspace(0.3, 1.0, 141) * model.maximum_field
    for rate in (model.cycle_averaged_rate, model.instantaneous_rate):
        assert np.all(np.diff(rate(field)) > 0)
    with pytest.warns(ModelRangeWarning, match="held at its value there"):
        above = model.instantaneous_rate(2 * model.maximum_field)
    assert above == pytest.approx(model.instantaneous_rate(model.maximum_field))
    with pytest.raises(ValueError, match="range"):
        replace(model, correction_range=0.0)

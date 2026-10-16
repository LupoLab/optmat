# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Named formulas and parameters reproduce independently recorded Luna values.

``tools/luna_reference.jl`` records ``tests/data/luna_master.json`` in Luna's
own environment. These reference checks do not require the private comparison pages.
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import EPS0, C
from optmat.exceptions import ModelRangeWarning
from optmat.gas.eos import number_density
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.models.registry import model
from optmat.gas.molecular_ionisation_rates import molecular_rate_model

ROOT = Path(__file__).parents[1]
LUNA = json.loads((ROOT / "tests" / "data" / "luna_master.json").read_text())
STATE = State(LUNA["pressure_bar"] * 1e5, LUNA["temperature_k"])
# Luna's dispersion formulas whose reference state optmat's copy shares (Börzsönyi's
# are referred to the ideal-gas density here and to the real-gas one in Luna).
FORMULAS = {
    "He": "luna_physdata/jct_He",
    "Ar": "luna_physdata/kingston_fit_Ar",
    "H2": "peck_1977/eq7",
    "D2": "peck_1977/eq7",
}


@pytest.mark.parametrize("gas", sorted(FORMULAS))
def test_the_named_formula_gives_lunas_refractivity(gas):
    entry = LUNA["gases"][gas]
    wavelength = np.array(entry["wavelength"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        ours = model(FORMULAS[gas]).refractivity(wavelength, STATE)
    np.testing.assert_allclose(ours, entry["refractivity"], rtol=2e-5)


@pytest.mark.parametrize("gas", sorted(LUNA["gases"]))
def test_lunas_density_is_the_real_gas_density(gas):
    ours = number_density(gas, STATE)
    assert LUNA["gases"][gas]["density"] == pytest.approx(ours, rel=1e-5)


@pytest.mark.parametrize("gas", ["He", "Ne", "Ar", "Kr", "Xe"])
def test_lunas_bare_ppt_rate_is_ppt_ilkov_in_the_tunnelling_regime(gas):
    record = LUNA["gases"][gas]["rate"]
    intensity = np.array(record["intensity_w_cm2"])
    field = np.sqrt(2 * intensity * 1e4 / (C * EPS0))
    form = ionisation_rate_model(gas, record["wavelength"], "ppt_ilkov")
    tunnelling = form.keldysh(field) < 0.5
    # The comparison is of the formula, so it runs above the barrier-suppression field.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        ours = form.instantaneous_rate(field[tunnelling])
    theirs = np.array(record["peak_bare"])[tunnelling]
    np.testing.assert_allclose(theirs, ours, rtol=0.05)


@pytest.mark.parametrize("gas", ["H2", "D2"])
def test_lunas_adk_rate_is_adk_hydrogenic(gas):
    record = LUNA["gases"][gas]["rate"]
    intensity = np.array(record["intensity_w_cm2"])
    field = np.sqrt(2 * intensity * 1e4 / (C * EPS0))
    form = molecular_rate_model(gas, record["wavelength"], "adk_hydrogenic")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ModelRangeWarning)
        ours = form.instantaneous_rate(field)
    np.testing.assert_allclose(record["peak_bare"], ours, rtol=5e-3)


def test_lunas_third_order_coefficients_are_solprops():
    solprop = json.loads((ROOT / "tests" / "data" / "solprop_gamma3.json").read_text())
    for gas, entry in LUNA["gases"].items():
        theirs = solprop["gamma3"][gas]["gamma3"]
        assert entry["gamma3"] == pytest.approx(theirs, rel=1e-7), gas


@pytest.mark.parametrize("gas", ["H2", "D2"])
def test_lunas_raman_parameters_are_solprops(gas):
    luna = LUNA["gases"][gas]["raman"]
    solprop = json.loads((ROOT / "tests" / "data" / "solprop_raman.json").read_text())
    rotational = solprop["raman"][gas]["rotational"]
    vibrational = solprop["raman"][gas]["vibrational"]
    pairs = [
        (luna["B"], rotational["b"]),
        (luna["D"], rotational["d"]),
        (luna["qJodd"], rotational["degeneracy_odd"]),
        (luna["qJeven"], rotational["degeneracy_even"]),
        (luna["delta_alpha"], rotational["polarisability_anisotropy"]),
        (luna["Arhor"], rotational["linewidth"]["inverse_density"]),
        (luna["omegav"], vibrational["frequency"]),
        (luna["dalpha_dq"], vibrational["polarisability_derivative"]),
        (luna["mu"], vibrational["reduced_mass"]),
        (luna["Brhov"], vibrational["linewidth"]["density"]),
        (luna["Arhov"], vibrational["linewidth"]["inverse_density"]),
    ]
    for theirs, ours in pairs:
        assert theirs == pytest.approx(ours, rel=1e-6)
    # The reference records differ in D2's rotational width per amagat.
    ratio = luna["Brhor"] / rotational["linewidth"]["density"]
    assert ratio == pytest.approx(1.0 if gas == "H2" else 1.20e8 / 1.32e8, rel=5e-3)

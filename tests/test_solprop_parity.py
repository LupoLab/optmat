# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The formulas here are the ones solprop uses today, coefficient for coefficient.

The fixture ``data/solprop_sellmeier.json`` holds solprop's susceptibility at each
model's reference state (``tools/solprop_parity.py``, run in solprop's environment).
At its own reference state a formula's refractivity is independent of how the density
is computed (solprop scales linearly; this package applies Lorentz-Lorenz and the
virial term, and inverts them at the same state), so the comparison tests the
coefficients alone.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from optmat.gas.models.formulas import published_models

FIXTURE = json.loads(
    (Path(__file__).parent / "data" / "solprop_sellmeier.json").read_text()
)

# solprop (gas, source) -> this package's model, and the relative agreement expected.
EXACT = 1e-11  # the same coefficients: round-off in the two evaluations
# Luna (and so solprop) multiplies Bideau-Mehu's printed 0.012055 into the printed B
# coefficients and rounds the products to six figures; this package keeps both as
# printed. The rounding is at most 3.3e-6 of the refractivity (krypton).
LUNA_ROUNDED = 5e-6
MODELS = {
    ("He", "ermolov"): ("luna_physdata/jct_He", EXACT),
    ("He", "borzsonyi"): ("borzsonyi_2008/sellmeier_He", EXACT),
    ("Ne", "borzsonyi"): ("borzsonyi_2008/sellmeier_Ne", EXACT),
    ("Ne", "bideau-mehu"): ("bideau_mehu_1981/sellmeier_Ne", LUNA_ROUNDED),
    ("Ar", "kingston"): ("luna_physdata/kingston_fit_Ar", EXACT),
    ("Ar", "borzsonyi"): ("borzsonyi_2008/sellmeier_Ar", EXACT),
    ("Ar", "bideau-mehu"): ("bideau_mehu_1981/sellmeier_Ar", LUNA_ROUNDED),
    ("Kr", "kingston"): ("luna_physdata/kingston_fit_Kr_1atm", EXACT),
    ("Kr", "borzsonyi"): ("borzsonyi_2008/sellmeier_Kr", EXACT),
    ("Kr", "bideau-mehu"): ("bideau_mehu_1981/sellmeier_Kr", LUNA_ROUNDED),
    ("Xe", "borzsonyi"): ("borzsonyi_2008/sellmeier_Xe", EXACT),
    ("N2", "borzsonyi"): ("borzsonyi_2008/sellmeier_N2", EXACT),
    ("H2", "peck"): ("peck_1977/eq7", EXACT),
    # solprop applies the hydrogen formula to deuterium (Orr 1936).
    ("D2", "peck"): ("peck_1977/eq7", EXACT),
    ("O2", "zhang"): ("kren_2011/kren_O2", EXACT),
    ("CH4", "quanfu-he"): ("he_2021/CH4", EXACT),
    ("SF6", "quanfu-he"): ("he_2021/SF6", EXACT),
    ("N2O", "quanfu-he"): ("he_2021/N2O", EXACT),
}
# Reference states that differ on purpose: solprop reads Ermolov's helium formula at
# 1 atm (TC-5), where Luna's legacy copy keeps Luna's 1 bar.
STATE_DIFFERS = {("He", "ermolov")}
ENTRIES = [e for e in FIXTURE["models"] if e["gas"] != "Air"]  # no air here


def test_every_solprop_gas_model_is_mapped():
    assert {(e["gas"], e["source"]) for e in ENTRIES} == set(MODELS)


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda e: f"{e['gas']}-{e['source']}")
def test_formula_matches_solprop_at_its_reference_state(entry):
    name, tolerance = MODELS[entry["gas"], entry["source"]]
    formula = published_models()[name]
    assert formula.state is not None
    wavelength = np.array(entry["wavelength_m"])
    expected = np.sqrt(1 + np.array(entry["chi"])) - 1
    refractivity = formula.refractivity(
        wavelength, formula.state, ideal=formula.ideal, warn=False
    )
    np.testing.assert_allclose(refractivity, expected, rtol=tolerance, atol=0)
    same_state = formula.state.pressure == pytest.approx(
        entry["pressure_bar"] * 1e5
    ) and formula.state.temperature == pytest.approx(entry["temperature"])
    assert same_state is ((entry["gas"], entry["source"]) not in STATE_DIFFERS)

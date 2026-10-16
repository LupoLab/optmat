# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The window crystals' recommended models (TC-545, TC-546)."""

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from optmat.datasets import datasets
from optmat.exceptions import OptmatError
from optmat.solid.models.formulas import FormulaModel
from optmat.solid.models.registry import default_model, default_reason, model

SNAPSHOT = Path(__file__).parents[1] / "data" / "solprop_solids.json"

# Each formula against its own paper's table, relative to the paper's medium (the model
# with its reference gas removed). Five-decimal tables are reproduced to half a unit of
# the fifth decimal (5.1e-6 with the formula's own rounding); the rows listed are the
# papers' own unexplained departures, recorded in the register (TC-546).
_TABLES = [
    ("daimon_2002/sellmeier_20C", "daimon_2002/n_20C", 2e-6, ()),
    ("daimon_2002/sellmeier_25C", "daimon_2002/n_25C", 2e-6, ()),
    ("malitson_1963/sellmeier", "malitson_1963/computed_24C", 5.1e-6, (253.7, 365.0)),
    ("dodge_1984/sellmeier_o", "dodge_1984/n_o_19C", 5.1e-6, ()),
    ("dodge_1984/sellmeier_e", "dodge_1984/n_e_19C", 5.1e-6, ()),
    # 0.334148 um sits on a rounding edge (5.4e-6); 0.486133 um is printed 1.4e-5 high.
    ("malitson_1964/sellmeier", "malitson_1964/computed_25C", 5.5e-6, (486.1,)),
    (
        "malitson_1962/sellmeier_o",
        "malitson_1962/computed_o_24C",
        5.1e-6,
        (3243.9, 3266.8, 3302.6, 3330.3, 3422.0, 3507.0),
    ),
]


def _dataset(key):
    return next(d for d in datasets() if d.key == key)


@pytest.mark.parametrize(("formula", "table", "tolerance", "outliers"), _TABLES)
def test_each_formula_reproduces_its_paper(formula, table, tolerance, outliers):
    chosen = model(formula)
    assert isinstance(chosen, FormulaModel)
    relative = dataclasses.replace(chosen, air=None)
    data = _dataset(table)
    wavelength = data.vacuum_wavelength()
    residual = np.real(relative.refractive_index(wavelength, warn=False))
    residual = residual - data.scaled_values()
    # The outliers are named by their printed (air) wavelength, within 0.5 nm.
    printed = data.abscissa_values()
    printed_nm = printed * (1e3 if data.abscissa_unit == "um" else 1.0)
    skip = np.zeros(wavelength.shape, dtype=bool)
    for nm in outliers:
        skip |= np.abs(printed_nm - nm) < 0.5
    assert skip.sum() == len(outliers)
    assert np.max(np.abs(residual[~skip])) < tolerance


def test_the_recommendations_are_absolute():
    # Each recommended formula carries its reference gas, which raises n by the gas's
    # refractivity times n: 3.7-6.5e-4 across the crystals' ranges (TC-546).
    for material, axis in [("CaF2", None), ("MgF2", "o"), ("BaF2", None)]:
        chosen = default_model(material, axis)
        assert isinstance(chosen, FormulaModel)
        assert chosen.air is not None
    calcium = default_model("CaF2")
    assert isinstance(calcium, FormulaModel)
    assert calcium.reference_gas == "nitrogen"
    assert "Daimon" in default_reason("CaF2")


def test_an_anisotropic_crystal_needs_its_axis():
    with pytest.raises(OptmatError, match="give the axis"):
        default_model("MgF2")
    assert default_model("MgF2", "e").axis == "e"
    # Sapphire has both rays since TC-651.
    assert default_model("sapphire", "e").axis == "e"


def test_solprop_uses_the_coefficients_as_absolute_indices():
    # solprop's CaF2, MgF2 and BaF2 are these papers' coefficients evaluated as
    # absolute indices at vacuum wavelengths; optmat's recommended models are higher by
    # the reference gas's refractivity times n (3.7-6.5e-4; TC-545, TC-546).
    # Sapphire is a different formula (Dodge's handbook formula at 20 C), so its
    # difference also holds that of the formulas (TC-546).
    snapshot = json.loads(SNAPSHOT.read_text(encoding="utf-8"))["models"]
    pairs = {
        ("CaF2", "daimon"): (None, 3.8e-4, 6.6e-4),
        ("MgF2", "dodge-o"): ("o", 3.6e-4, 4.6e-4),
        ("BaF2", "malitson"): (None, 3.8e-4, 4.5e-4),
        ("Al2O3", "malitson-o"): ("o", 4.4e-4, 6.5e-4),
    }
    for entry in snapshot:
        key = (entry["material"], entry["source"])
        if key not in pairs:
            continue
        axis, low, high = pairs[key]
        chosen = default_model(entry["material"], axis)
        wavelength = np.asarray(entry["wavelength_m"])
        start, stop = chosen.wavelength_range
        inside = (wavelength >= start) & (wavelength <= stop)
        ours = np.real(chosen.refractive_index(wavelength[inside]))
        difference = ours - np.asarray(entry["n"])[inside]
        assert low < difference.min() and difference.max() < high, key

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import re
from pathlib import Path

import numpy as np
import pytest

from optmat.solid.kerr import kerr_models
from optmat.solid.models.formulas import FormulaModel
from optmat.solid.models.registry import default_model, models

PAGE = (Path(__file__).parents[2] / "docs" / "silica.md").read_text()


def _table(name):
    block = PAGE.split(f"<!-- {name}:start -->")[1].split(f"<!-- {name}:end -->")[0]
    lines = [line for line in block.strip().splitlines()[2:] if line.strip()]
    return [[c.strip() for c in line.strip("|").split("|")] for line in lines]


@pytest.mark.parametrize("row", _table("vuv"), ids=lambda r: r[0])
def test_the_vuv_table_is_the_recommended_models(row):
    wavelength = float(row[0]) * 1e-9
    index = default_model("SiO2").refractive_index(wavelength)
    # Printed to 4 decimals in n and 3 significant digits in k.
    assert float(np.real(index)) == pytest.approx(float(row[1]), abs=5e-5)
    assert float(np.imag(index)) == pytest.approx(float(row[2]), rel=5e-3)


def test_the_window_table_is_the_five_models():
    names = [
        "leviton_2008/sellmeier",
        "malitson_1965/sellmeier",
        "arosa_2020/sellmeier_fit3",
        "franta_2016/udm_lithosil_q2",
        "fit/SiO2",
    ]
    rows = {row[0]: row[1:] for row in _table("window")}
    for column, name in enumerate(names):
        model = models("SiO2")[name]
        if name == "leviton_2008/sellmeier":
            assert isinstance(model, FormulaModel)
            model = model.at_temperature(293.15)  # the column is at 20 degrees C
        assert np.real(model.refractive_index(800e-9)) == pytest.approx(
            float(rows["n"][column]), abs=5e-6
        )
        assert model.group_index(800e-9) == pytest.approx(
            float(rows["group index"][column]), abs=5e-6
        )
        assert model.gvd(800e-9) * 1e27 == pytest.approx(
            float(rows["GVD (fs^2/mm)"][column]), abs=5e-3
        )


@pytest.mark.parametrize("row", _table("n2"), ids=lambda r: r[0])
def test_the_n2_table_is_the_stored_fit(row):
    model = kerr_models()["kerr/SiO2"]
    wavelength = float(row[0]) * 1e-9
    # Printed to 0.01e-20 m^2/W.
    assert model.n2(wavelength) * 1e20 == pytest.approx(float(row[1]), abs=5e-3)
    assert model.n2_sigma(wavelength) * 1e20 == pytest.approx(float(row[2]), abs=5e-3)
    assert model.n2_electronic(wavelength) * 1e20 == pytest.approx(
        float(row[3]), abs=5e-3
    )


def test_the_solprop_values_quoted_are_malitsons_pole():
    # "Re n = 0.000 at 110 nm and 2.94 at 120 nm": Malitson's printed formula evaluated
    # at vacuum wavelengths far outside its range, with a complex square root, as
    # solprop evaluates it (optmat's model converts to air wavelengths and gives 2.957
    # at 120 nm; the pole makes the difference large).
    b = np.array([0.6961663, 0.4079426, 0.8974794])
    poles = np.array([0.0684043, 0.1162414, 9.896161])
    lam = np.array([0.110, 0.120])[:, None]
    n = np.sqrt(1 + np.sum(b * lam**2 / (lam**2 - poles**2), axis=1) + 0j)
    assert re.search(r"0\.000 at 110 nm and 2\.94 at 120 nm", PAGE)
    assert np.real(n[0]) == pytest.approx(0.0, abs=5e-4)
    assert np.real(n[1]) == pytest.approx(2.94, abs=5e-3)

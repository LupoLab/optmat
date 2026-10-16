# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

from optmat.conditions import State
from optmat.datasets import datasets
from optmat.exceptions import ModelRangeWarning
from optmat.gas.models.formulas import published_models


def _fit_column(key):
    (dataset,) = [d for d in datasets() if d.key == key]
    return dataset


@pytest.mark.parametrize(
    ("model", "column", "rtol"),
    [
        # Rounding of the printed column (4 significant figures): 1e-4.
        ("bideau_mehu_1981/sellmeier_Ne", "bideau_mehu_1981/fit_Ne", 1.5e-4),
        ("bideau_mehu_1981/sellmeier_Ar", "bideau_mehu_1981/fit_Ar", 3e-4),
        # Next to the 146.9 nm resonance the column is sensitive to the rounding of the
        # printed pole (see the formula's check): 1 %.
        ("bideau_mehu_1981/sellmeier_Xe", "bideau_mehu_1981/fit_Xe", 1.1e-2),
    ],
)
def test_published_formula_reproduces_its_papers_column(model, column, rtol):
    formula = published_models()[model]
    fit = _fit_column(column)
    assert fit.state is not None
    computed = formula.refractivity(fit.vacuum_wavelength(), fit.state, warn=False)
    np.testing.assert_allclose(computed, fit.scaled_values(), rtol=rtol)


def test_bideau_mehu_krypton_formula_is_known_not_to_reproduce_its_column():
    formula = published_models()["bideau_mehu_1981/sellmeier_Kr"]
    fit = _fit_column("bideau_mehu_1981/fit_Kr")
    assert fit.state is not None
    computed = formula.refractivity(fit.vacuum_wavelength(), fit.state, warn=False)
    worst = np.max(np.abs(computed / fit.scaled_values() - 1))
    assert 5e-3 < worst < 1e-2  # recorded in the formula's check: 0.7 % at 140.4 nm


def test_formula_warns_outside_its_range():
    formula = published_models()["bideau_mehu_1981/sellmeier_Ar"]
    with pytest.warns(ModelRangeWarning, match="outside 140.4-567.7 nm"):
        formula.polarizability(100e-9)


@pytest.mark.parametrize(
    ("model", "wavelength", "state", "expected", "rtol"),
    [
        # Stone & Stejskal: n(633 nm, 101 325 Pa, 20 C) = 1.000 032 426 00(8). With
        # the refractivity virial term (Achtermann's B_R = -0.068 cm^6/mol^2; Stone
        # uses -0.061 at 293 K) the value closes to 0.8 ppm: the tolerance is the
        # difference of the two B_R, 1e-6 of n - 1.
        ("stone_2004/eq8", 632.99e-9, (101325.0, 293.15), 3242.600e-8, 1.5e-6),
        # Mansfield & Peck: (n - 1) = 3.4950e-5 at 5462.258 A, 0 C, 760 Torr; the
        # formula's rms residual is 3.6e-10 (1e-5 relative).
        ("mansfield_peck_1969/eq3", 546.2258e-9, (101325.0, 273.15), 3.4950e-5, 1e-4),
        # Peck & Khanna: 299.14e-6 at 546.2 nm, 0 C (rounding of the printed value).
        ("peck_1966/eq4", 546.2258e-9, (101325.0, 273.15), 299.14e-6, 5e-5),
        # Zhang et al.: argon 2605.925e-7 at 800 nm, 20 C (Table 2): rounding.
        ("zhang_2008/zhang_Ar", 800e-9, (101325.0, 293.15), 2605.925e-7, 1e-6),
        # Peck & Fisher: argon 2677.08e-7 at 5462.258 A, 15 C: rounding.
        ("peck_1964/peck_15C", 546.2258e-9, (101325.0, 288.15), 2677.08e-7, 1e-5),
        # Kren: oxygen n = 1.000251323 at 632.99 nm, 20 C: rounding.
        ("kren_2011/kren_O2", 632.99e-9, (101325.0, 293.15), 251.323e-6, 5e-6),
        # Borzsonyi xenon with C1 corrected: 707.09e-6 at 546.2258 nm and 1013.25 mbar
        # (checked with ideal-gas scaling, as the paper): rounding.
        (
            "borzsonyi_2008/sellmeier_Xe",
            546.2258e-9,
            (101325.0, 273.15),
            707.09e-6,
            2e-4,
        ),
    ],
)
def test_published_formula_reproduces_a_printed_value(
    model, wavelength, state, expected, rtol
):
    formula = published_models()[model]
    dn = float(formula.refractivity(wavelength, State(*state), ideal=formula.ideal))
    assert dn == pytest.approx(expected, rel=rtol, abs=0.0)

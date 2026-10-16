# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import numpy as np
import pytest

from optmat.gas.comparison import database_models, tabulated_entries
from optmat.gas.models.formulas import published_models


@pytest.mark.parametrize(
    ("database", "curated"),
    [
        ("Ar/Bideau-Mehu", "bideau_mehu_1981/sellmeier_Ar"),
        ("Ne/Bideau-Mehu", "bideau_mehu_1981/sellmeier_Ne"),
        ("Xe/Bideau-Mehu", "bideau_mehu_1981/sellmeier_Xe"),
        ("Kr/Bideau-Mehu", "bideau_mehu_1981/sellmeier_Kr"),
    ],
)
def test_database_bideau_mehu_agrees_with_the_curated_formula(database, curated):
    wavelength = np.linspace(150e-9, 500e-9, 20)
    ours = published_models()[curated].refractivity(wavelength, warn=False)
    theirs = database_models()[f"refractiveindex.info/{database}"].refractivity(
        wavelength, warn=False
    )
    # The database rounds B_i f_i products to six digits (2.50141e-3 for
    # 1.2055e-2 x 0.2075): a few 1e-7, 2e-6 near krypton's first pole.
    np.testing.assert_allclose(theirs, ours, rtol=3e-6)


def test_database_reads_ermolov_at_one_atmosphere_not_the_papers_bar():
    wavelength = np.array([200e-9, 800e-9])
    database = database_models()["refractiveindex.info/He/Ermolov"]
    paper = published_models()["ermolov_2015/sellmeier_He"]
    ratio = database.refractivity(wavelength) / paper.refractivity(wavelength)
    # The 1 atm / 1 bar real-gas density ratio, 1.01325: the database references the
    # formula to 1 atm, as solprop does (TC-5).
    np.testing.assert_allclose(ratio, 1 / 1.01325, rtol=1e-5)


def test_formula_types_cover_the_snapshot():
    kinds = {model.formula_type for model in database_models().values()}
    assert kinds == {"formula 2", "formula 5", "formula 6"}


def test_tabulated_entries_are_physical():
    for entry in tabulated_entries():
        assert np.all(entry.refractive_index > 1)
        assert np.all(entry.wavelength > 50e-9)

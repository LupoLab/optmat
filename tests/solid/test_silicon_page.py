# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

from pathlib import Path

import numpy as np
import pytest

from optmat.datasets import datasets
from optmat.solid.models.formulas import FormulaModel
from optmat.solid.models.registry import default_model, models

PAGE = (Path(__file__).parents[2] / "docs" / "silicon.md").read_text(encoding="utf-8")
KEYED = {d.key: d for d in datasets()}
EV_NM = 1239.84198  # hc in eV nm


def _table(name):
    block = PAGE.split(f"<!-- {name}:start -->")[1].split(f"<!-- {name}:end -->")[0]
    lines = [line for line in block.strip().splitlines()[2:] if line.strip()]
    return [[c.strip() for c in line.strip("|").split("|")] for line in lines]


def _k(key, wavelength):
    """Return a dataset's k at a vacuum wavelength, interpolated in log k."""
    dataset = KEYED[key]
    x = dataset.vacuum_wavelength()
    y = dataset.scaled_values()
    if dataset.quantity == "absorption_coefficient":
        y = y * 100 * x / (4 * np.pi)  # base e, cm-1
    order = np.argsort(x)
    return float(np.exp(np.interp(wavelength, x[order], np.log(y[order]))))


def _value(key, wavelength):
    dataset = KEYED[key]
    x = dataset.vacuum_wavelength()
    order = np.argsort(x)
    return float(np.interp(wavelength, x[order], dataset.scaled_values()[order]))


@pytest.mark.parametrize("row", _table("overview"), ids=lambda r: r[0])
def test_the_overview_is_the_recommended_model(row):
    index = default_model("Si").refractive_index(float(row[0]) * 1e-9)
    # Printed to 4 decimals in n and 3 significant digits in k.
    assert float(np.real(index)) == pytest.approx(float(row[1]), abs=5e-5)
    assert float(np.imag(index)) == pytest.approx(float(row[2]), rel=5e-3)


def test_the_window_table_is_the_data_and_models():
    rows = {row[0]: row[1:] for row in _table("window")}
    all_models = models("Si")
    formula = all_models["frey_2006/sellmeier_Si"]
    assert isinstance(formula, FormulaModel)
    frey = formula.at_temperature(293.15)
    newad2 = all_models["franta_2017/udm_c_si"]
    refit = default_model("Si")
    table = KEYED["frey_2006/n_Si_295K"]
    slope = KEYED["frey_2006/dndt_Si_295K"].scaled_values()
    taken = table.scaled_values() - slope * (295 - 293.15)  # to 20 C with their dn/dT
    for label, um in (("n at 1.5 µm", 1.5), ("n at 2.0 µm", 2.0)):
        wavelength = um * 1e-6
        row = [float(v) for v in rows[label]]
        at = np.argmin(np.abs(table.vacuum_wavelength() - wavelength))
        expected = [
            taken[at],
            float(np.real(frey.refractive_index(wavelength))),
            _value("li_1980/n_293K", wavelength),
            _value("franta_2017/rii_20C_n", wavelength),
            float(np.real(newad2.refractive_index(wavelength))),
            float(np.real(refit.refractive_index(wavelength))),
        ]
        # Printed to 5 decimals.
        np.testing.assert_allclose(row, expected, atol=5e-6)
    gvd = [v for v in rows["GVD at 2.0 µm (fs^2/mm)"] if v]
    # Printed to 0.1 fs^2/mm.
    np.testing.assert_allclose(
        [float(v) for v in gvd],
        [m.gvd(2e-6) * 1e27 for m in (frey, newad2, refit)],
        atol=0.05,
    )


@pytest.mark.parametrize("row", _table("edge"), ids=lambda r: r[0])
def test_the_edge_table_is_the_data_and_models(row):
    wavelength = float(row[0]) * 1e-9
    expected = [
        _k("schinke_2015/alpha_295K", wavelength),
        _k("green_2008/alpha_300K", wavelength),
        _k("franta_2017/rii_20C_k", wavelength),
        float(
            np.imag(models("Si")["franta_2017/udm_c_si"].refractive_index(wavelength))
        ),
        float(np.imag(default_model("Si").refractive_index(wavelength))),
    ]
    # Printed to 3 significant digits.
    np.testing.assert_allclose([float(v) for v in row[1:]], expected, rtol=5e-3)


@pytest.mark.parametrize("row", _table("uv"), ids=lambda r: r[0])
def test_the_ultraviolet_table_is_the_data_and_the_model(row):
    wavelength = EV_NM / float(row[0]) * 1e-9
    index = default_model("Si").refractive_index(wavelength)
    expected = [
        _value("jellison_1992/n_100", wavelength),
        _value("jellison_1992/k_100", wavelength),
        _value("aspnes_1983/n_Si", wavelength),
        _value("aspnes_1983/k_Si", wavelength),
        _value("franta_2017/rii_20C_n", wavelength),
        _value("franta_2017/rii_20C_k", wavelength),
        float(np.real(index)),
        float(np.imag(index)),
    ]
    # Printed to 3 decimals.
    np.testing.assert_allclose([float(v) for v in row[1:]], expected, atol=5e-4)

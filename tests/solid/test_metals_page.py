# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The metals page's table is the recommended table models (TC-619, TC-646, TC-647)."""

import json
from pathlib import Path

import numpy as np
import pytest

from optmat.solid.models.registry import default_model

PAGE = (Path(__file__).parents[2] / "docs" / "metals.md").read_text()
_BLOCK = PAGE.split("<!-- metals:start -->")[1].split("<!-- metals:end -->")[0]
_ROWS = [
    [c.strip() for c in line.strip("|").split("|")]
    for line in _BLOCK.strip().splitlines()[2:]
]
_WAVELENGTHS = (100e-9, 400e-9, 800e-9, 1030e-9, 1550e-9, 10e-6)
SNAPSHOT = Path(__file__).parents[1] / "data" / "solprop_solids.json"


@pytest.mark.parametrize("row", _ROWS, ids=lambda r: r[0])
def test_each_row_is_the_recommended_table(row):
    model = default_model(row[0])
    assert model.name == f"table/{row[0]}"
    low, high = (float(v) for v in row[1].split("-"))
    assert model.wavelength_range[0] * 1e9 == pytest.approx(low, abs=0.05)
    assert model.wavelength_range[1] * 1e9 == pytest.approx(high, abs=0.5)
    for cell, wavelength in zip(row[2:], _WAVELENGTHS, strict=True):
        if not cell:
            continue
        real, imaginary = (float(v.rstrip("i")) for v in cell.split(" + "))
        index = complex(model.refractive_index(wavelength))
        # Printed to three decimals.
        assert index.real == pytest.approx(real, abs=5e-4)
        assert index.imag == pytest.approx(imaginary, abs=5e-4)


@pytest.mark.parametrize(
    ("material", "ranges"),
    [
        ("Ag", [(25e-9, 70e-9), (207e-9, 5e-6)]),
        ("Al", [(33e-6, 200e-6)]),
    ],
)
def test_solprop_agrees_where_the_source_is_the_same(material, ranges):
    # solprop's tables are refractiveindex.info's copies of the same papers in these
    # ranges; 1 % in n allows its four-digit rounding of n ~ 0.05 and the two
    # interpolants between the points, 0.2 % in k their difference where k is large.
    entry = next(
        e
        for e in json.loads(SNAPSHOT.read_text())["models"]
        if e["material"] == material
    )
    model = default_model(material)
    wavelength = np.asarray(entry["wavelength_m"])
    for low, high in ranges:
        inside = (wavelength >= low) & (wavelength <= high)
        assert inside.any()
        ours = model.refractive_index(wavelength[inside])
        n = np.asarray(entry["n"])[inside]
        k = np.asarray(entry["k"])[inside]
        np.testing.assert_allclose(np.real(ours), n, rtol=1.1e-2)
        np.testing.assert_allclose(np.imag(ours), k, rtol=2e-3)

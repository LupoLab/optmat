# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The crystals page's table is the recommended models (TC-546, TC-651)."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from optmat.solid.models.registry import default_model

PAGE = (Path(__file__).parents[2] / "docs" / "crystals.md").read_text(encoding="utf-8")
_BLOCK = PAGE.split("<!-- crystals:start -->")[1].split("<!-- crystals:end -->")[0]
_ROWS = [
    [c.strip() for c in line.strip("|").split("|")]
    for line in _BLOCK.strip().splitlines()[2:]
]

_HISTORY = (
    Path(__file__).parents[1] / "data/historical_documentation_tables.md"
).read_text(encoding="utf-8")
_HISTORY_BLOCK = _HISTORY.split("<!-- crystals:start -->")[1].split(
    "<!-- crystals:end -->"
)[0]
_HISTORY_ROWS = {
    (cells[0], cells[1]): cells
    for line in _HISTORY_BLOCK.strip().splitlines()[2:]
    for cells in [[c.strip() for c in line.strip("|").split("|")]]
}


def _solprop(material, source):
    # solprop's formula from the snapshot's poles, as solprop evaluates it: an absolute
    # index at the vacuum wavelength.
    path = Path(__file__).parents[1] / "data" / "solprop_solids.json"
    entry = next(
        e
        for e in json.loads(path.read_text(encoding="utf-8"))["models"]
        if (e["material"], e["source"]) == (material, source)
    )
    square = 1.0  # lambda = 1 um
    return math.sqrt(
        1
        + sum(
            b * square / (square - pole**2)
            for b, pole in entry["sellmeier_b_lambda_um"]
        )
    )


_SOLPROP = {("CaF2", ""): "daimon", ("MgF2", "o"): "dodge-o", ("BaF2", ""): "malitson"}
_SOLPROP[("Al2O3", "o")] = "malitson-o"
_SOLPROP[("N-BK7", "")] = "schott"


@pytest.mark.parametrize("row", _ROWS, ids=lambda r: f"{r[0]}{r[1]}")
def test_each_row_is_the_recommended_model(row):
    material, axis, name = row[0], row[1], row[2].strip("`")
    chosen = default_model(material, axis or None)
    assert chosen.name == name
    low, high = (float(v) for v in row[3].split("-"))
    assert chosen.wavelength_range[0] * 1e9 == pytest.approx(low, abs=0.5)
    assert chosen.wavelength_range[1] * 1e9 == pytest.approx(high, abs=0.5)
    # Printed to 5 decimals.
    for column, wavelength in ((6, 500e-9), (7, 1e-6)):
        index = float(np.real(chosen.refractive_index(wavelength)))
        assert index == pytest.approx(float(row[column]), abs=5e-6)
    if (material, axis) in _SOLPROP:
        difference = _solprop(material, _SOLPROP[(material, axis)])
        difference -= float(np.real(chosen.refractive_index(1e-6)))
        # Printed to two significant digits.
        printed = float(
            _HISTORY_ROWS[(material, axis)][8].replace("\N{MINUS SIGN}", "-")
        )
        assert difference == pytest.approx(printed, abs=5e-6)

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""solprop's recorded values, the baseline of solprop SEP 3 (OEP 1, Phase 3; TC-787).

Every fixture is recorded at one named solprop commit. Its transport values and
densities are compared with optmat's: the molar masses agree exactly, the densities by
the equations' gas constants alone, and the viscosities and heat-capacity ratios by the
amounts in ``tests/data/historical_documentation_tables.md``.
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from optmat.conditions import State
from optmat.constants import K_B, N_A
from optmat.exceptions import StateError
from optmat.gas.eos import air_equation_of_state, equation_of_state
from optmat.gas.transport import transport_model

DATA = Path(__file__).parent / "data"
TRANSPORT = json.loads((DATA / "solprop_transport.json").read_text(encoding="utf-8"))
EOS = json.loads((DATA / "solprop_eos.json").read_text(encoding="utf-8"))
R_CODATA = K_B * N_A


def test_every_fixture_records_one_solprop_commit():
    commits = {
        json.loads(path.read_text(encoding="utf-8"))["solprop_commit"]
        for path in DATA.glob("solprop_*.json")
    }
    assert commits == {"d57741490dbcb05ceffca91a0cdddd069963ec69"}


def _page_rows():
    page = (DATA / "historical_documentation_tables.md").read_text(encoding="utf-8")
    section = page.split("<!-- transport-solprop:start -->")[1].split(
        "<!-- transport-solprop:end -->"
    )[0]
    return {
        cells[0]: cells
        for line in section.splitlines()
        if line.startswith("| ") and not line.startswith("| Gas")
        for cells in [[c.strip() for c in line.strip("|").split("|")]]
    }


@pytest.mark.parametrize("gas", sorted(TRANSPORT["gases"]))
def test_the_transport_against_solprop(gas):
    entry = TRANSPORT["gases"][gas]
    model = transport_model(gas)
    assert model.molar_mass == entry["molar_mass_kg_per_mol"]
    temperatures = np.array(entry["temperatures_k"])
    viscosity = model.viscosity(temperatures) / np.array(entry["viscosity_pa_s"]) - 1
    cp = model.heat_capacity(temperatures)
    ratio = cp / (cp - R_CODATA) / np.array(entry["heat_capacity_ratio"]) - 1
    worst = int(np.argmax(np.abs(viscosity)))
    cells = _page_rows()[gas]
    match = re.match(r"(\S+) at (\d+) K", cells[2])
    assert match is not None
    printed, at = match.groups()
    # The page's figures, to their two digits.
    assert viscosity[worst] == pytest.approx(
        float(printed.replace("\N{MINUS SIGN}", "-")), rel=0.05
    )
    assert temperatures[worst] == float(at)
    assert np.max(np.abs(ratio)) == pytest.approx(float(cells[3]), rel=0.05)


@pytest.mark.parametrize("gas", sorted(EOS["gases"]))
def test_the_densities_differ_by_the_gas_constant_alone(gas):
    eos = air_equation_of_state() if gas == "Air" else equation_of_state(gas)
    # solprop forms the density from each equation's own gas constant, optmat from Z
    # and CODATA's k_B (TC-758): their ratio is R_eos / R_CODATA exactly.
    expected = eos.gas_constant / R_CODATA - 1
    for row in EOS["gases"][gas]:
        state = State(
            pressure=row["pressure_bar"] * 1e5, temperature=row["temperature_k"]
        )
        if gas == "SF6" and row["pressure_bar"] == 30.0:
            # solprop returns a liquid's density (8 times the ideal gas) where SF6 is
            # condensed; optmat refuses the state (OEP 1, section 4.2; TC-787).
            assert row["number_density_m3"] > 7 * state.pressure / (
                K_B * state.temperature
            )
            with pytest.raises(StateError, match="condensed"):
                eos.number_density(state)
            continue
        ratio = eos.number_density(state) / row["number_density_m3"] - 1
        # Beyond the gas constant, the two flashes' iteration tolerances: up to 3e-10,
        # for N2O at 30 bar and 273 K, near its vapour pressure (31 bar).
        assert ratio == pytest.approx(expected, abs=1e-9)
        assert row["round_trip_bar"] == pytest.approx(row["pressure_bar"], rel=1e-14)

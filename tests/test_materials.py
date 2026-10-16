# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import pytest

from optmat.datasets import datasets, load_source
from optmat.exceptions import DataError, OptmatError
from optmat.gas.gases import Gas
from optmat.materials import Material, as_material, material, materials


def test_every_gas_is_a_registered_gas_material():
    for gas in Gas:
        assert materials()[gas].kind == "gas"


def test_aliases_resolve_to_their_id():
    assert as_material("fused silica") == "SiO2"
    assert as_material("SiO2") == "SiO2"
    assert material("quartz").id == "SiO2_quartz"


def test_an_unknown_name_is_refused():
    with pytest.raises(OptmatError, match="unknown material"):
        as_material("unobtainium")


def test_axes_match_the_optical_symmetry():
    for entry in materials().values():
        expected = {"isotropic": 0, "uniaxial": 2, "biaxial": 3}[entry.optics]
        assert len(entry.axes) == expected, entry.id
    with pytest.raises(OptmatError, match="axes"):
        Material("X", "x", "crystal", "X", "uniaxial", ("o",), ())


def test_every_dataset_names_a_registered_id():
    registered = set(materials())
    assert {d.material for d in datasets()} <= registered


def _record(tmp_path, material_line, extra=""):
    record = tmp_path / "demo_2000.toml"
    record.write_text(
        'citation = "demo"\nkind = "measurement"\n[[datasets]]\nid = "d"\n'
        f"{material_line}\n"
        'quantity = "polarizability"\nunit = "au"\nabscissa_unit = "nm"\n'
        'points = [["1", "500", "11.2"]]\n'
        f"{extra}",
        encoding="utf-8",
    )
    return record


def test_a_record_must_use_the_id_not_an_alias(tmp_path):
    with pytest.raises(DataError, match="not a registered id"):
        load_source(_record(tmp_path, 'material = "fused silica"'))


def test_a_solid_has_no_gas_state(tmp_path):
    record = _record(tmp_path, 'material = "SiO2"', 'state = "1atm_0C"\n')
    with pytest.raises(DataError, match="no gas state"):
        load_source(record)


def test_an_axis_must_be_one_of_the_materials(tmp_path):
    good = _record(tmp_path, 'material = "MgF2"', '[datasets.parameters]\naxis = "o"\n')
    assert load_source(good).datasets[0].parameters["axis"] == "o"
    bad = _record(tmp_path, 'material = "CaF2"', '[datasets.parameters]\naxis = "o"\n')
    with pytest.raises(DataError, match="axis"):
        load_source(bad)

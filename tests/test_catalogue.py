# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Single records of the curated catalogue, read without loading all of it."""

import subprocess
import sys

import numpy as np
import pytest

from optmat import _catalogue
from optmat.datasets import DATA_DIRECTORY, datasets, load_source, sources


@pytest.fixture
def unloaded():
    # Read records as a fresh process would, before any record is read.
    _catalogue.source.cache_clear()
    sources.cache_clear()
    yield
    _catalogue.source.cache_clear()
    sources.cache_clear()


def _same(alone, catalogued):
    assert (alone.key, alone.unit, alone.quantity, alone.material) == (
        catalogued.key,
        catalogued.unit,
        catalogued.quantity,
        catalogued.material,
    )
    assert np.array_equal(alone.abscissa_values(), catalogued.abscissa_values())
    assert np.array_equal(alone.scaled_values(), catalogued.scaled_values())


@pytest.mark.parametrize("source_id", ["ogawa_1958", "malitson_1965", "huber_1974"])
def test_a_record_read_alone_is_the_catalogues(unloaded, source_id):
    # The catalogue takes the record already read rather than reading it again, and
    # that record is what the loader reads from the file.
    alone = _catalogue.source(source_id)
    assert sources()[source_id] is alone
    loaded = load_source(DATA_DIRECTORY / "sources" / f"{source_id}.toml")
    assert (alone.citation, alone.formulas) == (loaded.citation, loaded.formulas)
    assert len(alone.datasets) == len(loaded.datasets)
    for one, other in zip(alone.datasets, loaded.datasets, strict=True):
        _same(one, other)


def test_missing_records_and_datasets_raise_key_errors(unloaded):
    with pytest.raises(KeyError, match="no_such_source"):
        _catalogue.source("no_such_source")
    with pytest.raises(KeyError, match="no_such_source/x"):
        _catalogue.dataset("no_such_source/x")
    with pytest.raises(KeyError, match="ogawa_1958/no_such_dataset"):
        _catalogue.dataset("ogawa_1958/no_such_dataset")
    assert _catalogue.records().get("no_such_source") is None


def test_a_loaded_catalogue_serves_its_own_records():
    catalogued = sources()["ogawa_1958"]
    assert _catalogue.source("ogawa_1958") is catalogued
    assert _catalogue.records()["ogawa_1958"] is catalogued


def test_lookups_by_key_filter_as_the_catalogue_does():
    key = "ogawa_1958/absolute"
    _same(_catalogue.dataset(key), {d.key: d for d in datasets()}[key])
    assert key in _catalogue.by_key()
    assert key in _catalogue.by_key("O3", quantity="absorption_cross_section")
    assert key not in _catalogue.by_key("Ar")
    assert key not in _catalogue.by_key(quantity="refractivity")
    assert list(_catalogue.by_key("O3")) == [d.key for d in datasets("O3")]
    assert len(_catalogue.by_key(quantity="refractivity")) == len(
        datasets(quantity="refractivity")
    )
    assert list(_catalogue.records()) == list(sources())
    assert len(_catalogue.records()) == len(sources())


def test_the_formula_records_are_every_record_with_formulas():
    # Only the records with a [[formulas]] table are read; they must be exactly the
    # sources with formulas, in the catalogue's order.
    expected = [s.id for s in sources().values() if s.formulas]
    assert [s.id for s in _catalogue.formula_sources()] == expected


def test_solprops_entry_points_never_load_the_catalogue():
    # The models a propagator asks for read only their own records (TC-810).
    code = """
import numpy as np
import optmat.datasets
from optmat.conditions import ONE_ATM_0C
from optmat.gas.models.registry import default_model as gas_model
from optmat.gas.ionisation import ionisation_parameters
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.raman_response import raman_response
from optmat.gas.transport import transport_model
from optmat.solid.bundle import silica_bundle
from optmat.solid.models.registry import default_model as solid_model
gas_model("N2").refractive_index(np.array([800e-9]), ONE_ATM_0C)
solid_model("CaF2").refractive_index(np.array([800e-9]))
solid_model("Ag").refractive_index(np.array([800e-9]))
silica_bundle()
ionisation_rate_model("Ar", 800e-9)
ionisation_parameters("Ar")
raman_response("N2", 293.15)
transport_model("Ar")
print(optmat.datasets.sources.cache_info().currsize)
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "0"


@pytest.mark.parametrize(
    "quantity", ["mole_fraction", "refractivity", "mixture_compressibility_factor"]
)
def test_the_datasets_of_a_quantity_are_the_catalogues(quantity):
    # Only the records whose text sets the quantity are read; they must hold exactly
    # the catalogue's datasets of it, in its order.
    expected = [d.key for d in datasets(quantity=quantity)]
    assert [d.key for d in _catalogue.quantity_datasets(quantity)] == expected


def test_croaks_entry_points_read_only_their_own_records():
    # croak (pulse retrieval) asks for the recommended solids and gases and for air;
    # each reads its own records, neither the catalogue nor every formula record.
    code = """
import numpy as np
import optmat._catalogue
import optmat.datasets
from optmat.conditions import State
from optmat.gas.air import air_model
from optmat.gas.models.registry import default_model as gas_model
from optmat.solid.models.registry import default_model as solid_model
wavelength = np.array([800e-9])
for name in ("SiO2", "N-BK7", "CaF2", "BaF2", "Si"):
    solid_model(name).refractive_index(wavelength)
solid_model("MgF2", axis="o").refractive_index(wavelength)
gas_model("He").refractive_index(wavelength, State(1e5, 293.15), ideal=True)
air_model().refractive_index(wavelength, State(1e5, 293.15), ideal=True)
print(optmat.datasets.sources.cache_info().currsize,
      optmat._catalogue.formula_sources.cache_info().currsize)
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.split() == ["0", "0"]

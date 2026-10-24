# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The catalogue's index: committed, current, and checked against the records."""

import json
import warnings

import pytest

from optmat import _catalogue, _index
from optmat.datasets import datasets, sources
from optmat.exceptions import CatalogueIndexWarning
from optmat.gas.gases import Gas


def test_the_committed_index_is_the_records():
    # Fails after a record is added or edited until tools/catalogue_index.py is run.
    assert _index.INDEX_FILE.read_text(encoding="utf-8") == _index.dumps(_index.build())


def test_the_committed_index_raises_no_warning():
    _index.index.cache_clear()
    with warnings.catch_warnings():
        warnings.simplefilter("error", CatalogueIndexWarning)
        _index.index()


@pytest.mark.parametrize(
    ("material", "quantity"),
    [("He", None), ("SiO2", None), (None, "mole_fraction"), ("Ar", "refractivity")],
)
def test_a_search_finds_what_the_whole_catalogue_holds(material, quantity):
    expected = [
        d.key
        for record in sources().values()
        for d in record.datasets
        if (material is None or d.material == material)
        and (quantity is None or d.quantity == quantity)
    ]
    assert [d.key for d in datasets(material, quantity=quantity)] == expected


def test_the_formula_sources_of_some_materials_are_the_catalogues():
    gases = frozenset(g.value for g in Gas)
    expected = [
        s.id
        for s in sources().values()
        if any(str(f.get("material")) in gases for f in s.formulas)
    ]
    assert [s.id for s in _catalogue.formula_sources(gases)] == expected


RECORD = 'citation = "synthetic"\nkind = "theory"\n'
DATASET = (
    '[[datasets]]\nid = "{id}"\nmaterial = "{material}"\nquantity = "refractivity"\n'
)


@pytest.fixture
def catalogue(tmp_path, monkeypatch):
    # A catalogue of two records and its index, in place of the curated one.
    directory = tmp_path / "sources"
    directory.mkdir()
    (directory / "a.toml").write_text(
        RECORD + DATASET.format(id="x", material="He"), encoding="utf-8"
    )
    (directory / "b.toml").write_text(
        RECORD + DATASET.format(id="y", material="Ar"), encoding="utf-8"
    )
    monkeypatch.setattr(_index, "_SOURCES", directory)
    monkeypatch.setattr(_index, "INDEX_FILE", tmp_path / "catalogue_index.json")
    _index.INDEX_FILE.write_text(_index.dumps(_index.build()), encoding="utf-8")
    _index.index.cache_clear()
    yield directory
    _index.index.cache_clear()


def _index_quietly():
    with warnings.catch_warnings():
        warnings.simplefilter("error", CatalogueIndexWarning)
        return _index.index()


def test_an_index_in_step_is_used_without_reading_the_records(catalogue, monkeypatch):
    def unread(path):
        raise AssertionError(f"read {path}")

    monkeypatch.setattr(_index, "entry", unread)
    assert _index_quietly()["a"]["datasets"] == [["He", "refractivity"]]
    assert _index.sources_with("Ar") == ("b",)


def test_an_edit_that_keeps_the_entries_is_read_again_quietly(catalogue):
    (catalogue / "a.toml").write_text(
        RECORD + 'notes = "longer"\n' + DATASET.format(id="x", material="He"),
        encoding="utf-8",
    )
    assert _index_quietly()["a"]["size"] == (catalogue / "a.toml").stat().st_size
    assert _index.sources_with("He") == ("a",)


@pytest.mark.parametrize("change", ["edited", "added", "removed", "unindexed"])
def test_a_stale_index_warns_and_the_records_are_used(catalogue, change):
    # An edit is seen by its size, so this one changes it (He to SiO2).
    if change == "edited":
        (catalogue / "a.toml").write_text(
            RECORD + DATASET.format(id="x", material="SiO2"), encoding="utf-8"
        )
    elif change == "added":
        (catalogue / "c.toml").write_text(
            RECORD + DATASET.format(id="z", material="Ne"), encoding="utf-8"
        )
    elif change == "removed":
        (catalogue / "b.toml").unlink()
    else:
        _index.INDEX_FILE.write_text(
            json.dumps({"format": 0, "sources": {}}), encoding="utf-8"
        )
    with pytest.warns(CatalogueIndexWarning, match="out of date"):
        entries = _index.index()
    assert entries == _index.build()["sources"]


def test_a_missing_index_reads_every_record(catalogue):
    _index.INDEX_FILE.unlink()
    stale = r"for 2 source record\(s\) \(a, b\)"
    with pytest.warns(CatalogueIndexWarning, match=stale):
        assert _index.index() == _index.build()["sources"]


def test_a_long_list_of_stale_records_is_shortened(catalogue):
    for n in range(6):
        (catalogue / f"n{n}.toml").write_text(RECORD, encoding="utf-8")
    with pytest.warns(CatalogueIndexWarning, match=r"for 6 .*n4, \.\.\."):
        _index.index()

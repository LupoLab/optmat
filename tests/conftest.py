# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Test session configuration.

Shape and dtype annotations are checked at runtime during tests only: the jaxtyping
import hook instruments every ``optmat`` module imported after it with
beartype, configured with the PEP 484 numeric tower as in solprop.
"""

from collections.abc import Mapping

import pytest
from jaxtyping import install_import_hook

JAXTYPING_HOOK = install_import_hook(
    "optmat",
    "beartype.beartype(conf=beartype.BeartypeConf(is_pep484_tower=True))",
)

import optmat  # noqa: E402, F401 - imported after the hook so it is checked


class _Overlay(Mapping):
    """A mapping that serves some entries itself and reads any other from a fallback."""

    def __init__(self, served, fallback):
        self._served = served
        self._fallback = fallback

    def __getitem__(self, key):
        if key in self._served:
            return self._served[key]
        return self._fallback(key)

    def __iter__(self):
        return iter(self._served)

    def __len__(self):
        return len(self._served)


@pytest.fixture
def curated(monkeypatch):
    """Return a function that adds synthetic datasets (and records) to the catalogue.

    The models look datasets and records up by key (``optmat._catalogue``), so a
    contract test serves synthetic ones there; any other key is read from the curated
    catalogue as usual. Iterating serves only the synthetic entries.
    """
    from optmat import _catalogue  # noqa: PLC0415 - imported after the hook
    from optmat.materials import as_material  # noqa: PLC0415 - as above

    real_dataset, real_source = _catalogue.dataset, _catalogue.source

    def use(*datasets, records=()):
        def by_key(material=None, *, quantity=None):
            wanted = None if material is None else as_material(material)
            served = {
                d.key: d
                for d in datasets
                if (wanted is None or d.material == wanted)
                and (quantity is None or d.quantity == quantity)
            }

            def fallback(key):
                found = real_dataset(key)
                if (wanted is not None and found.material != wanted) or (
                    quantity is not None and found.quantity != quantity
                ):
                    raise KeyError(key)
                return found

            return _Overlay(served, fallback)

        def dataset(key):
            return by_key()[key]

        sources = _Overlay({r.id: r for r in records}, real_source)
        monkeypatch.setattr(_catalogue, "by_key", by_key)
        monkeypatch.setattr(_catalogue, "dataset", dataset)
        monkeypatch.setattr(_catalogue, "records", lambda: sources)
        monkeypatch.setattr(_catalogue, "source", sources.__getitem__)

    return use

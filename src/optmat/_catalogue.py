# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Single records of the curated catalogue, read without loading all of it.

:func:`optmat.datasets.sources` loads every source record and its data files, which
takes seconds; most models need a few records. These read one record with
:func:`optmat.datasets.load_source`, the loader the catalogue uses, so a record read
here is the catalogue's own. Each record is read once per process: the catalogue,
when loaded, takes the records already read here.
"""

import functools
from collections.abc import Iterator, Mapping

from optmat._index import sources_with_formulas
from optmat.datasets import (
    DATA_DIRECTORY,
    Dataset,
    Source,
    datasets,
    load_source,
    sources,
)
from optmat.materials import as_material

_SOURCES = DATA_DIRECTORY / "sources"


@functools.cache
def source(source_id: str) -> Source:
    """Return one curated source, loading only its record.

    Parameters
    ----------
    source_id : str
        The source's id, its record's file name.

    Returns
    -------
    Source
        The source, errata applied; each dataset is built when first asked for.

    Raises
    ------
    KeyError
        If there is no such record.
    """
    path = _SOURCES / f"{source_id}.toml"
    if not path.is_file():
        raise KeyError(source_id)
    return load_source(path, defer=True)


def dataset(key: str) -> Dataset:
    """Return one curated dataset, reading only its source's record and its data file.

    Parameters
    ----------
    key : str
        ``source/dataset``.

    Returns
    -------
    Dataset
        The dataset.

    Raises
    ------
    KeyError
        If there is no such dataset.
    """
    source_id, _, dataset_id = key.partition("/")
    try:
        return source(source_id).dataset(dataset_id)
    except KeyError:
        raise KeyError(key) from None


class _Records(Mapping[str, Source]):
    """The curated sources by id; a lookup reads only that record."""

    def __getitem__(self, source_id: str) -> Source:
        """Return a source by id.

        Parameters
        ----------
        source_id : str
            The source's id.

        Returns
        -------
        Source
            The source.
        """
        return source(source_id)

    def __iter__(self) -> Iterator[str]:
        """Iterate over the ids, which loads the whole catalogue.

        Returns
        -------
        Iterator of str
            The ids, in source order.
        """
        return iter(sources())

    def __len__(self) -> int:
        """Return the number of sources, which loads the whole catalogue.

        Returns
        -------
        int
            The number of curated sources.
        """
        return len(sources())


_RECORDS: Mapping[str, Source] = _Records()


def records() -> Mapping[str, Source]:
    """Return the curated sources by id, read on lookup.

    ``records()[source_id]`` is ``sources()[source_id]`` without loading every record.

    Returns
    -------
    Mapping of str to Source
        The sources; iterating over it loads the whole catalogue.
    """
    return _RECORDS


class _ByKey(Mapping[str, Dataset]):
    """The curated datasets by key, filtered as ``datasets()`` filters them.

    A lookup reads only its source's record.

    Parameters
    ----------
    material : str or None
        Only this material (its id or an alias).
    quantity : str or None
        Only this quantity.
    """

    def __init__(self, material: str | None, quantity: str | None) -> None:
        """Hold the filters."""
        self._material = material
        self._wanted = None if material is None else as_material(material)
        self._quantity = quantity

    def __getitem__(self, key: str) -> Dataset:
        """Return a dataset by key.

        Parameters
        ----------
        key : str
            ``source/dataset``.

        Returns
        -------
        Dataset
            The dataset.

        Raises
        ------
        KeyError
            If there is no such dataset, or it is not of the material or quantity.
        """
        found = dataset(key)
        if (self._wanted is not None and found.material != self._wanted) or (
            self._quantity is not None and found.quantity != self._quantity
        ):
            raise KeyError(key)
        return found

    def __iter__(self) -> Iterator[str]:
        """Iterate over the keys, which loads the whole catalogue.

        Returns
        -------
        Iterator of str
            The keys, in source order.
        """
        return (d.key for d in datasets(self._material, quantity=self._quantity))

    def __len__(self) -> int:
        """Return the number of datasets, which loads the whole catalogue.

        Returns
        -------
        int
            The number of matching datasets.
        """
        return len(datasets(self._material, quantity=self._quantity))


def by_key(
    material: str | None = None, *, quantity: str | None = None
) -> Mapping[str, Dataset]:
    """Return the curated datasets by key, read on lookup.

    ``by_key(material, quantity=q)[key]`` is ``{d.key: d for d in datasets(material,
    quantity=q)}[key]`` without loading every record.

    Parameters
    ----------
    material : str, optional
        Only this material: its id or an alias.
    quantity : str, optional
        Only this quantity.

    Returns
    -------
    Mapping of str to Dataset
        The datasets; iterating over it loads the whole catalogue.
    """
    return _ByKey(material, quantity)


@functools.cache
def formula_sources(materials: frozenset[str] | None = None) -> tuple[Source, ...]:
    """Return the curated sources that record published formulas, in source order.

    Only the records the catalogue's index lists with formulas are read.

    Parameters
    ----------
    materials : frozenset of str, optional
        Only the sources with a formula of one of these material ids.

    Returns
    -------
    tuple of Source
        The sources, in the catalogue's order.
    """
    return tuple(source(s) for s in sources_with_formulas(materials))


@functools.cache
def quantity_datasets(quantity: str) -> tuple[Dataset, ...]:
    """Return the curated datasets of one quantity, reading only the records with it.

    ``quantity_datasets(q)`` is ``datasets(quantity=q)``, which reads only the records
    the catalogue's index lists for the quantity.

    Parameters
    ----------
    quantity : str
        The quantity (e.g. ``"mole_fraction"``).

    Returns
    -------
    tuple of Dataset
        The datasets, in source order.

    Examples
    --------
    >>> {d.quantity for d in quantity_datasets("mole_fraction")}
    {'mole_fraction'}
    """
    return datasets(quantity=quantity)

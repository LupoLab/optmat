# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Single records of the curated catalogue, read without loading all of it.

:func:`optmat.datasets.sources` loads every source record and its data files, which
takes seconds; most models need a few records. These read one record with
:func:`optmat.datasets.load_source`, the loader the catalogue uses, so a record read
here is the catalogue's own. Once the catalogue is loaded they take its records.
"""

import functools
import re
from collections.abc import Iterator, Mapping

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
        The source, with its datasets read and errata applied.

    Raises
    ------
    KeyError
        If there is no such record.
    """
    if sources.cache_info().currsize:
        return sources()[source_id]
    path = _SOURCES / f"{source_id}.toml"
    if not path.is_file():
        raise KeyError(source_id)
    return load_source(path)


def dataset(key: str) -> Dataset:
    """Return one curated dataset, loading only its source's record.

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
    try:
        record = source(key.split("/", 1)[0])
    except KeyError:
        raise KeyError(key) from None
    for found in record.datasets:
        if found.key == key:
            return found
    raise KeyError(key)


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
def formula_sources() -> tuple[Source, ...]:
    """Return every curated source that records published formulas, in source order.

    Only the records whose text has a ``[[formulas]]`` table are read (a note may
    mention one; a record read without formulas is left out).

    Returns
    -------
    tuple of Source
        The sources, in the catalogue's order.
    """
    candidates = (
        source(path.stem)
        for path in sorted(_SOURCES.glob("*.toml"))
        if "[[formulas]]" in path.read_text(encoding="utf-8")
    )
    return tuple(record for record in candidates if record.formulas)


@functools.cache
def quantity_datasets(quantity: str) -> tuple[Dataset, ...]:
    """Return the curated datasets of one quantity, reading only the records with it.

    ``quantity_datasets(q)`` is ``datasets(quantity=q)`` without loading every record:
    only the records whose text sets ``quantity = "<q>"`` are read.

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
    pattern = re.compile(rf'quantity\s*=\s*"{re.escape(quantity)}"')
    return tuple(
        found
        for path in sorted(_SOURCES.glob("*.toml"))
        if pattern.search(path.read_text(encoding="utf-8"))
        for found in source(path.stem).datasets
        if found.quantity == quantity
    )

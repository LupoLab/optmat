# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The diatomic electronic states of ``nist_webbook_2026`` against the stored pages."""

import gzip
import hashlib
import re
from html.parser import HTMLParser

import numpy as np
import pytest

from optmat.constants import E_CHARGE, HBAR, C
from optmat.datasets import DATA_DIRECTORY, Dataset, sources

DIRECTORY = DATA_DIRECTORY / "external" / "nist_webbook"
GASES = ("H2", "D2", "N2", "O2", "CO")
# A printed number, possibly in H&H's parentheses or brackets (doubtful or deduced).
NUMBER = re.compile(r"^[\[\(]*([0-9]+(?:\.[0-9]+)?)[\]\)]*$")
EV_PER_CM1 = 2 * np.pi * HBAR * C * 100 / E_CHARGE


class _Rows(HTMLParser):
    """Collect the cells of every table row of a page."""

    def __init__(self) -> None:
        super().__init__()
        self.rows: list[list[str]] = []
        self._row: list[str] | None = None
        self._cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: object) -> None:
        if tag == "tr":
            self._row = []
        elif tag in {"td", "th"} and self._row is not None:
            self._cell = []
        elif tag == "br" and self._cell is not None:
            self._cell.append(" ")

    def handle_endtag(self, tag: str) -> None:
        if tag in {"td", "th"} and self._cell is not None and self._row is not None:
            self._row.append(re.sub(r"\s+", " ", "".join(self._cell)).strip())
            self._cell = None
        elif tag == "tr" and self._row is not None:
            self.rows.append(self._row)
            self._row = None

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell.append(data)


def _page(gas: str) -> bytes:
    return gzip.decompress(
        (DIRECTORY / f"{gas}_diatomic_constants.html.gz").read_bytes()
    )


def _states(gas: str) -> dict[str, list[str]]:
    parser = _Rows()
    parser.feed(_page(gas).decode("utf-8"))
    # The constants table has 13 columns: State, Te, ..., Trans., nu00.
    return {row[0]: row for row in parser.rows if len(row) >= 13}


def _datasets(gas: str) -> dict[str, Dataset]:
    record = sources()["nist_webbook_2026"]
    return {
        d.id: d
        for d in record.datasets
        if d.material == gas and d.quantity == "electronic_state_energy"
    }


@pytest.mark.parametrize("gas", GASES)
def test_the_pages_are_the_ones_recorded(gas: str):
    provenance = (DIRECTORY / "PROVENANCE.md").read_text(encoding="utf-8")
    digest = hashlib.sha256(_page(gas)).hexdigest()
    assert f"`{gas}_diatomic_constants.html.gz`" in provenance
    assert digest in provenance


@pytest.mark.parametrize("gas", GASES)
def test_the_states_are_the_pages_strings(gas: str):
    rows = _states(gas)
    datasets = _datasets(gas)
    assert datasets
    for dataset in datasets.values():
        parameters = dataset.parameters
        row = rows[str(parameters["state"])]
        (label,) = dataset.abscissae
        assert label == parameters["state"]
        match = NUMBER.match(row[1].split()[0])
        assert match is not None
        assert dataset.values == (match.group(1),)
        assert parameters["kind"] == "Te"
        assert dataset.use == "reproduced"
        if "nu00_cm1" in parameters:
            nu00 = NUMBER.match(row[12].split()[0])
            assert nu00 is not None
            assert parameters["nu00_cm1"] == nu00.group(1)
            assert parameters["transition"] == row[11]
        for rule in ("one_photon", "two_photon"):
            assert parameters[rule] in {"allowed", "forbidden", "spin_forbidden"}


def _origin(gas: str, state: str) -> float:
    (dataset,) = [d for d in _datasets(gas).values() if d.id == f"state_{gas}_{state}"]
    return float(str(dataset.parameters["nu00_cm1"]))


def _lowest(source: str, dataset_id: str) -> float:
    (dataset,) = [d for d in sources()[source].datasets if d.id == dataset_id]
    # The datasets' abscissa is the photon energy in eV, read as a vacuum wavelength.
    return float(np.max(dataset.vacuum_wavelength()))


@pytest.mark.parametrize(
    ("gas", "state", "source", "dataset_id"),
    [
        ("H2", "B", "sakamoto_2010", "lyman_H2"),
        ("CO", "A", "sakamoto_2010", "discrete_CO"),
    ],
)
def test_the_band_origins_are_the_lowest_measured_bands(
    gas: str, state: str, source: str, dataset_id: str
):
    # The 0-0 band against the lowest band of a dipole (e,e) spectrum, printed to
    # 1 meV (8 cm^-1): agreement within 2 meV.
    measured = EV_PER_CM1 / (100 * _lowest(source, dataset_id))
    assert _origin(gas, state) * EV_PER_CM1 == pytest.approx(measured, abs=2e-3)


def test_n2s_b_state_through_the_a_state_is_chan_et_als_first_band():
    # b -> a and a -> a' and a' <-> X: the b(0) - X(0) interval, 12.500 eV in
    # chan_1993_n2 (printed to 1 meV).
    t0 = sum(_origin("N2", state) for state in ("b", "a", "ap"))
    measured = EV_PER_CM1 / (100 * _lowest("chan_1993_n2", "discrete_present"))
    assert t0 * EV_PER_CM1 == pytest.approx(measured, abs=2e-3)

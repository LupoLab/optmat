# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import csv

import numpy as np
import pytest

from optmat.datasets import DATA_DIRECTORY, sources

DIRECTORY = DATA_DIRECTORY / "external" / "nist_asd"
FILE = DIRECTORY / "ionization_energies.csv"
IONS = ("Ne", "Ar", "Kr", "Xe")


def _rows() -> dict[str, list[str]]:
    # NIST's CSV writes each cell as ="..." so that spreadsheets keep it as text.
    with FILE.open(newline="", encoding="utf-8") as handle:
        rows = [
            [cell[2:-1] if cell.startswith('="') else cell for cell in row]
            for row in csv.reader(handle)
            if row and row[0].startswith('="')
        ]
    return {row[1]: row for row in rows}


def _levels(gas: str, spectrum: str = "II") -> list[dict[str, str]]:
    # NIST's tab-delimited output: a header row, then every cell in double quotes.
    path = DIRECTORY / f"{gas}_{spectrum}_levels.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quotechar='"'))


def test_the_inline_ionisation_energies_are_nists_strings():
    rows = _rows()
    record = sources()["nist_asd_2026"]
    energies = [d for d in record.datasets if d.quantity == "ionization_energy"]
    assert len(energies) == len(rows) == 10
    for dataset in energies:
        (spectrum,) = dataset.abscissae
        row = rows[spectrum]
        assert dataset.values == (row[9],)
        assert dataset.sigmas == (row[11],)
        assert dataset.parameters["final_state"] in {row[7], "the bare nucleus"}
        # A value NIST brackets (semi-empirical) or parenthesises (ab initio) is a
        # check, not a primary datum.
        assert (dataset.use == "check") == (row[8] in {"(", "["})


def test_the_first_ionisation_energy_of_argon():
    (argon,) = [d for d in sources()["nist_asd_2026"].datasets if d.id == "ie_Ar_I"]
    assert argon.scaled_values()[0] == pytest.approx(15.7596119, abs=0)


def test_the_inline_levels_are_nists_strings():
    record = sources()["nist_asd_2026"]
    levels = {
        d.id: d
        for d in record.datasets
        if d.quantity == "level_energy" and d.material.endswith("+")
    }
    assert set(levels) == {f"level_{gas}+_2P12" for gas in IONS}
    for gas in IONS:
        dataset = levels[f"level_{gas}+_2P12"]
        ground, upper, *_ = _levels(gas)
        # The ground term of the ion: np5 2P*, J = 3/2 at zero and J = 1/2 above it.
        assert (ground["Term"], ground["J"]) == ("2P*", "3/2")
        assert float(ground["Level (cm-1)"]) == 0.0
        assert ground["Configuration"] == upper["Configuration"]
        parameters = dataset.parameters
        assert parameters["configuration"] == upper["Configuration"]
        assert (parameters["term"], parameters["j"]) == (upper["Term"], upper["J"])
        assert (parameters["term"], parameters["j"]) == ("2P*", "1/2")
        assert dataset.material == f"{gas}+"
        assert dataset.abscissae == (f"{gas} II",)
        assert dataset.unit == "cm-1"
        assert (upper["Prefix"], upper["Suffix"]) == ("", "")
        assert dataset.values == (upper["Level (cm-1)"],)
        assert dataset.sigmas == (upper["Uncertainty (cm-1)"],)
        # Where NIST gives no uncertainty (Kr II, Xe II) none is recorded.
        stated = upper["Uncertainty (cm-1)"] != ""
        assert dataset.uncertainty.stated == stated
        assert np.isfinite(dataset.standard_uncertainty()[0]) == stated
        if upper["Reference"]:
            assert upper["Reference"] in dataset.notes


ATOMS = ("He", *IONS)
# The neutral levels that bound the off-resonant third-order response (TC-862): per
# atom, the dataset suffix and the level it must be.
NEUTRAL = {
    "He": {
        "1s2p_1P1": ("1s.2p", "1P*", "1"),
        "1s3p_1P1": ("1s.3p", "1P*", "1"),
        "1s2s_1S0": ("1s.2s", "1S", "0"),
    },
}


def _value(row: dict[str, str]) -> float:
    return float(row["Level (cm-1)"])


def _lowest(rows: list[dict[str, str]], *, odd: bool, j: str) -> dict[str, str]:
    return next(
        row
        for row in rows
        if row["J"] == j and ("*" in row["Term"]) == odd and _value(row) > 0
    )


def test_the_record_holds_the_ions_neutrals_and_energies():
    record = sources()["nist_asd_2026"]
    levels = [d for d in record.datasets if d.quantity == "level_energy"]
    # Four ion levels; three for helium and four for each heavier atom.
    assert len(levels) == len(IONS) + 3 + 4 * len(IONS)
    assert len(record.datasets) == 10 + len(levels)


@pytest.mark.parametrize("gas", ATOMS)
def test_the_neutral_levels_are_nists_strings(gas: str):
    record = sources()["nist_asd_2026"]
    datasets = {
        d.id.removeprefix(f"level_{gas}_"): d
        for d in record.datasets
        if d.quantity == "level_energy" and d.material == gas
    }
    rows = [row for row in _levels(gas, "I") if row["J"] and row["Level (cm-1)"]]
    if gas == "He":
        expected = {
            suffix: next(
                row
                for row in rows
                if (row["Configuration"], row["Term"], row["J"]) == level
            )
            for suffix, level in NEUTRAL["He"].items()
        }
    else:
        expected = {
            "odd_J1": _lowest(rows, odd=True, j="1"),
            **{f"even_J{j}": _lowest(rows, odd=False, j=j) for j in ("0", "1", "2")},
        }
    assert set(datasets) == set(expected)
    for suffix, row in expected.items():
        dataset = datasets[suffix]
        parameters = dataset.parameters
        assert (parameters["configuration"], parameters["term"], parameters["j"]) == (
            row["Configuration"],
            row["Term"],
            row["J"],
        )
        assert dataset.abscissae == (f"{gas} I",)
        assert dataset.values == (row["Level (cm-1)"],)
        assert dataset.sigmas == (row["Uncertainty (cm-1)"],)
        # A level NIST brackets (semi-empirical) is a check, as for the energies.
        assert (dataset.use == "check") == (row["Prefix"] == "[")
        stated = row["Uncertainty (cm-1)"] != ""
        assert dataset.uncertainty.stated == stated
        if row["Reference"]:
            assert row["Reference"] in dataset.notes


@pytest.mark.parametrize("gas", ATOMS)
def test_the_one_photon_level_is_the_stored_resonance_line(gas: str):
    # The resonance-line files give the same upper level as the levels file.
    suffix = "1s2p_1P1" if gas == "He" else "odd_J1"
    (level,) = [
        d
        for d in sources()["nist_asd_2026"].datasets
        if d.id == f"level_{gas}_{suffix}"
    ]
    with (DIRECTORY / f"{gas}_I_lines.tsv").open(newline="", encoding="utf-8") as f:
        lines = list(csv.DictReader(f, delimiter="\t", quotechar='"'))
    upper = {
        float(row["Ek(cm-1)"].strip("[]"))
        for row in lines
        if row["J_k"] == "1" and row["term_k"] == level.parameters["term"]
    }
    # The two files print the same NIST level to different digits: 0.01 cm^-1.
    assert min(abs(e - float(level.values[0])) for e in upper) < 0.01


@pytest.mark.parametrize(
    ("gas", "low", "high"),
    [("Ne", 76, 5000), ("Ar", 110, 5000), ("Kr", 125, 5000), ("Xe", 150, 5000)],
)
def test_the_near_infrared_lines_are_kept_as_returned(gas, low, high):
    # PROVENANCE.md records each file's sha256 as returned (TC-983, TC-985); vacuum
    # wavelengths throughout, the window starting where the resonance-line file's ends.
    import hashlib  # noqa: PLC0415

    path = DIRECTORY / f"{gas}_I_lines_nir.tsv"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    provenance = (DIRECTORY / "PROVENANCE.md").read_text(encoding="utf-8")
    assert f"| `{gas}_I_lines_nir.tsv` | {gas} I | {low}-{high} |" in provenance
    assert digest in provenance
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quotechar='"'))
    assert "ritz_wl_vac(nm)" in rows[0]
    # Every printed wavelength is a vacuum wavelength inside the window: show_av=1
    # wrote wavenumbers above 1000 nm as negative numbers (TC-985). A cell may carry a
    # + or * mark.
    cells = [
        r[k].rstrip("+*?") for r in rows for k in ("obs_wl_vac(nm)", "ritz_wl_vac(nm)")
    ]
    wavelengths = [float(c) for c in cells if c]
    assert len(wavelengths) >= len(rows)
    assert low <= min(wavelengths) and max(wavelengths) <= high

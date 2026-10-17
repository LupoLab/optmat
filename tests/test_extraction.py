# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The committed data files are exactly what the extraction scripts make of the library.

The library-marked test needs the curators' paper library (``$PAPER_LIBRARY``) and runs
with ``poe test-library``; the others check the committed tree alone.
"""

import csv
import importlib.util
import io
import itertools
import os
import sys
from pathlib import Path

import pytest

from optmat.datasets import DATA_DIRECTORY, VALUES_DIRECTORY, data_file, sources

TOOLS = Path(__file__).resolve().parents[1] / "tools" / "extract"
SCRIPTS = sorted(p for p in TOOLS.glob("*.py") if not p.name.startswith("_"))


def _load(script: Path):
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    spec = importlib.util.spec_from_file_location(f"extract_{script.stem}", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _render(relative: str, points) -> str:
    library = _load(TOOLS / "_library.py")
    header = library.columns(relative)
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(header)
    for p in points:
        writer.writerow(
            [p.row, p.x, p.value, p.sigma, p.printed]
            + [p.extra[name] for name in header[len(library.COLUMNS) :]]
        )
    return buffer.getvalue()


@pytest.mark.library
@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.stem)
def test_committed_files_match_a_fresh_extraction(script):
    if "PAPER_LIBRARY" not in os.environ:
        pytest.skip("set PAPER_LIBRARY to the paper library to re-run the extraction")
    library = _load(TOOLS / "_library.py")
    for relative, points in _load(script).extract().items():
        committed = library.path(relative).read_text(encoding="utf-8")
        assert committed == _render(relative, points), relative


def test_every_data_file_belongs_to_exactly_one_dataset():
    # A file lives in the directory of its dataset's family, so this also checks
    # that every file is where its quantity puts it.
    files = {str(path) for path in VALUES_DIRECTORY.rglob("*.csv")}
    referenced = [
        str(data_file(quantity, dataset_file))
        for source in sources().values()
        for quantity, dataset_file in _dataset_files(source.id)
    ]
    assert sorted(referenced) == sorted(set(referenced)), "a file used twice"
    assert set(referenced) == files


def test_every_data_file_has_the_general_columns():
    for path in VALUES_DIRECTORY.rglob("*.csv"):
        header = path.read_text(encoding="utf-8").split("\n", 1)[0].split(",")
        assert header[:5] == ["row", "x", "value", "sigma", "printed"], path.name


def _dataset_files(source_id: str) -> list[tuple[str, str]]:
    import tomllib  # noqa: PLC0415 - only this helper reads the raw records

    record = tomllib.loads(
        (DATA_DIRECTORY / "sources" / f"{source_id}.toml").read_text(encoding="utf-8")
    )
    return [
        (entry["quantity"], entry["file"])
        for entry in record.get("datasets", [])
        if "file" in entry
    ]


def test_every_script_names_an_existing_source():
    for script in SCRIPTS:
        assert script.stem in sources(), f"{script.name} has no source record"


def test_the_trecx_two_photon_files_and_uncertainties_match_the_archive():
    # The two-photon tRecX runs are archived in data/external/trecx/two_photon, not the
    # library, so their extraction runs without it; the stated uncertainties are the
    # script's rule (TC-1048).
    library = _load(TOOLS / "_library.py")
    script = _load(TOOLS / "travers_2026_trecx_two_photon.py")
    for relative, points in script.extract().items():
        committed = library.path(relative).read_text(encoding="utf-8")
        assert committed == _render(relative, points), relative
    script.checks()
    by_id = {d.id: d for d in sources()["travers_2026_trecx_two_photon"].datasets}
    for dataset, rows in script.uncertainties().items():
        stated = by_id[dataset].uncertainty
        for row, value in rows.items():
            expected = stated.relative_rows.get(row, stated.relative)
            assert expected == pytest.approx(value, rel=0.01, abs=0.0), (dataset, row)


def test_the_pyroraptor_files_and_uncertainties_match_the_archived_tables():
    # pyroraptor's helium runs are archived in data/external/pyroraptor, not the
    # library; the stated uncertainties are the script's rule, and its checks are what
    # the record says of the runs (TC-1092).
    library = _load(TOOLS / "_library.py")
    script = _load(TOOLS / "travers_2026_pyroraptor.py")
    for relative, points in script.extract().items():
        committed = library.path(relative).read_text(encoding="utf-8")
        assert committed == _render(relative, points), relative
    script.checks()
    by_id = {d.id: d for d in sources()["travers_2026_pyroraptor"].datasets}
    computed = script.uncertainties()
    assert set(computed) == set(by_id)
    for dataset, rows in computed.items():
        stated = by_id[dataset].uncertainty.relative_rows
        assert set(stated) == set(rows), dataset
        for row, value in rows.items():
            # Stated to three figures (two below 1e-3).
            assert stated[row] == pytest.approx(value, rel=0.01, abs=5e-6), (
                dataset,
                row,
            )


def test_the_trecx_files_and_uncertainties_match_the_archived_tables():
    # The tRecX calculations are archived in data/external/trecx, not the library, so
    # their extraction runs without it; the stated uncertainties of the 6 fs ionised
    # fractions and the static rates are what the runs' rule gives (TC-812, TC-813).
    library = _load(TOOLS / "_library.py")
    script = _load(TOOLS / "travers_2026_trecx.py")
    for relative, points in script.extract().items():
        committed = library.path(relative).read_text(encoding="utf-8")
        assert committed == _render(relative, points), relative
    by_id = {d.id: d for d in sources()["travers_2026_trecx"].datasets}
    for family, (_, target) in script.SIX_FS.items():
        stated = by_id[f"ionised_gauss6_{target}"].uncertainty.relative_rows
        computed = script.six_fs_uncertainty(family)
        assert set(stated) == {row for row, u in computed.items() if u is not None}
        for row, value in stated.items():
            assert value == pytest.approx(computed[row], rel=0.01, abs=1e-5), row
    rules = [
        (f"ionised_gauss10_{nm}nm_{target}", script.six_fs_uncertainty(family, table))
        for (nm, table), (family, (_, target)) in itertools.product(
            script.WAVELENGTH_TABLES.items(), script.WAVELENGTH.items()
        )
    ]
    rules += [
        (f"ionised_gauss30_{nm}nm_{target}", script.six_fs_uncertainty(family, table))
        for (nm, table), (family, (_, target)) in itertools.product(
            script.LONG_TABLES.items(), script.WAVELENGTH.items()
        )
    ]
    rules += [
        (f"ionised_wallace_{phase}_{target}", script.wallace_uncertainty(family, phase))
        for phase, (family, (_, target)) in itertools.product(
            script.WALLACE_TABLES, script.WAVELENGTH.items()
        )
    ]
    rules += [
        (f"ionised_gauss15_{nm}nm_H", script.six_fs_uncertainty("h_m0", table))
        for nm, table in script.HYDROGEN_TABLES.items()
    ]
    rules += [
        (
            f"ionised_gauss6_{target}",
            script.six_fs_uncertainty(family, script.KR_XE_TABLE),
        )
        for family, (_, target) in script.KR_XE.items()
    ]
    rules += [
        (
            f"ionised_wallace_{phase}_Kr_m0",
            script.wallace_uncertainty("kr_m0", phase, script.WALLACE_KR_TABLES),
        )
        for phase in script.WALLACE_KR_TABLES
    ]
    rules += [
        (
            f"ionised_wallace_{phase}_Xe_m0",
            script.wallace_uncertainty("xe_m0", phase, script.WALLACE_XE_TABLES),
        )
        for phase in script.WALLACE_XE_TABLES
    ]
    rules += [
        (
            f"{measure}_sin2_{part}_{target}",
            script.six_fs_uncertainty(family, script.SIN2_TABLES[part], column),
        )
        for (part, family), (_, target) in script.SIN2_FAMILIES.items()
        for measure, column in script.SIN2_MEASURES.items()
    ]
    rules += [
        (
            f"{measure}_ipmatched_{energy}_{target}",
            script.ip_matched_uncertainty(energy, family, column),
        )
        for (energy, _table), (measure, column), (family, (_, target)) in (
            itertools.product(
                script.IP_MATCHED_TABLES.items(),
                script.IP_MATCHED_MEASURES.items(),
                script.WAVELENGTH.items(),
            )
        )
    ]
    two = script.he_two_electron_uncertainty()
    rules += [
        (script.he_two_electron_dataset(setup), dict.fromkeys(runs, two))
        for setup, runs in script.HE_TWO_ELECTRON_RUNS.items()
    ]
    rules += [
        (
            "hacc_normloss_gauss6_390nm_He",
            dict.fromkeys(script.HE_HACC_RUNS, script.he_hacc_uncertainty()),
        )
    ]
    rules += [
        (
            f"{measure}_ipmatched_390nm_{energy}_He",
            script.he_one_electron_uncertainty(energy, column),
        )
        for energy in script.HE_ONE_ELECTRON_TABLES
        for measure, column in script.IP_MATCHED_MEASURES.items()
    ]
    rules += [
        (script.hacc_dataset(setup), script.hacc_uncertainty(setup))
        for setup in script.HACC_RUNS
    ]
    rules += [
        (script.xe_rmt_dataset(m, suffix), script.xe_rmt_uncertainty(m, suffix))
        for suffix in script.XE_RMT_SETS
        for m in (0, 1)
    ]
    rules += [
        (f"static_{target}", script.kr_xe_static_uncertainty(material))
        for material, target in script.KR_XE.values()
    ]
    for name, computed in rules:
        dataset = by_id[name]
        stated = dataset.uncertainty.relative_rows
        assert set(stated) == {row for row, u in computed.items() if u is not None}
        # The rows below the floor are the excluded ones (TC-826).
        assert set(dataset.excluded_rows) == set(computed) - set(stated)
        for row, value in stated.items():
            assert value == pytest.approx(computed[row], rel=0.01, abs=1e-5), (
                name,
                row,
            )
    for m in (0, 1):
        stated = by_id[f"static_Ar_m{m}"].uncertainty.relative_rows
        computed = script.static_uncertainty(m)
        assert set(stated) == set(computed)
        for row, value in stated.items():
            assert value == pytest.approx(computed[row], rel=0.01, abs=1e-5), row

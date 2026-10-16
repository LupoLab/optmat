# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Source contracts with synthetic records, independent of the curated library."""

import dataclasses
import math

import numpy as np
import pytest
from scipy.special import ndtri

import optmat.datasets as records
from optmat.conditions import State
from optmat.constants import C
from optmat.datasets import Uncertainty, data_file, load_source, parse_state
from optmat.exceptions import DataError
from optmat.vocabulary import Quantity

# These values are invented to test record semantics, not material measurements.
HEADER = 'citation = "synthetic test"\nkind = "theory"\n'
ENTRY = (
    '[[datasets]]\nid = "sample"\nmaterial = "Ar"\n'
    'quantity = "polarizability"\nunit = "au"\nabscissa_unit = "nm"\n'
    'points = [["a", "500", "10.00", "0.3"], ["b", "inf", "20.0", ""]]\n'
)


def _source(tmp_path, text=HEADER + ENTRY):
    path = tmp_path / "synthetic.toml"
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def sample(tmp_path):
    return load_source(_source(tmp_path)).datasets[0]


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"values": ("1",)}, "columns of unequal length"),
        ({"sigmas": ()}, "columns of unequal length"),
        ({"columns": {"temperature": ("1",)}}, "columns of unequal length"),
        ({"rows": ("a", "a")}, "duplicate row ids"),
        ({"wavelength_medium": "air"}, "unknown wavelength medium"),
        ({"use": "fit"}, "unknown use"),
        ({"use": "reproduced"}, "must name their source"),
        ({"quantity": "imaginary_quantity"}, "unknown quantity"),
        ({"unit": "m"}, "unit 'm'"),
        ({"quantity": "refractivity", "unit": "1"}, "needs the state"),
        ({"columns": {"temperature": ("20", "20")}}, "extra column needs a unit"),
        (
            {"columns": {"x_upper": ("600", "inf")}, "column_units": {"x_upper": "m"}},
            "x_upper is in the abscissa's unit",
        ),
        ({"abscissa": "duration"}, "is tabulated against"),
        ({"abscissa_unit": "furlong"}, "unknown wavelength unit"),
        (
            {"abscissa_unit": "eV", "wavelength_medium": "standard_air"},
            "only wavelengths and wavenumbers",
        ),
        ({"uncertainty": Uncertainty(relative_rows={"absent": 0.1})}, "unknown rows"),
        ({"parameters": {"form": "powder"}}, "form 'powder'"),
        ({"parameters": {"assumed": "yes"}}, "assumed 'yes'"),
        ({"parameters": {"base": "2"}}, "base '2'"),
        ({"parameters": {"site": "interface"}}, "site 'interface'"),
        ({"parameters": {"index_medium": "nitrogen"}}, "index_medium 'nitrogen'"),
    ],
)
def test_inconsistent_dataset_metadata_is_rejected(sample, changes, message):
    with pytest.raises(DataError, match=message):
        dataclasses.replace(sample, **changes)


def test_reproduced_values_retain_their_origin(sample):
    copied = dataclasses.replace(
        sample, use="reproduced", reproduced_from=("original",)
    )
    assert copied.reproduced_from == ("original",)
    assert copied.use == "reproduced"
    assert copied.values == sample.values


@pytest.mark.parametrize("part", ["not_a_part", "zero_point_vibrational"])
def test_third_order_parts_cannot_be_polarizability_corrections(sample, part):
    with pytest.raises(DataError, match="unknown part"):
        dataclasses.replace(
            sample,
            quantity="hyperpolarizability",
            parameters={"process": "kerr", "convention": "T", "part": part},
        )


def test_a_polarizability_correction_uses_its_own_vocabulary(sample):
    correction = dataclasses.replace(
        sample, quantity="polarizability_correction", parameters={"part": "thermal"}
    )
    assert correction.parameters["part"] == "thermal"
    with pytest.raises(DataError, match="unknown part"):
        dataclasses.replace(correction, parameters={"part": "electronic"})


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ("ambient", "unknown state"),
        ({"temperature_k": 300}, "needs a pressure"),
        ({"pressure_pa": 1e5}, "needs a pressure"),
        ([], "cannot read state"),
    ],
)
def test_incomplete_states_are_rejected(spec, message):
    with pytest.raises(DataError, match=message):
        parse_state(spec, "synthetic/sample")


def test_explicit_state_units_and_si_precedence():
    assert parse_state({"pressure_atm": 2.0, "temperature_c": 20.0}, "test") == State(
        202650.0, 293.15
    )
    assert parse_state(
        {
            "pressure_pa": 2e5,
            "pressure_atm": 1.0,
            "temperature_k": 310.0,
            "temperature_c": 20.0,
        },
        "test",
    ) == State(2e5, 310.0)
    assert parse_state(None, "test") is None


@pytest.mark.parametrize("confidence", [90, 95, 99])
def test_confidence_intervals_convert_using_normal_quantiles(confidence):
    spec = Uncertainty(absolute=1.0, level=f"{confidence} % confidence")
    # Independent inverse-normal calculation; the stored decimal factors are rounded.
    expected = 1 / ndtri((1 + confidence / 100) / 2)
    assert spec.to_standard == pytest.approx(expected, rel=2e-15)


def test_uncertainties_combine_in_quadrature_with_scale_and_row_overrides(sample):
    changed = dataclasses.replace(
        sample,
        scale=2.0,
        uncertainty=Uncertainty(
            column=True,
            relative=0.1,
            relative_rows={"b": 0.2},
            absolute=0.4,
            level="2 sigma",
        ),
    )
    # First row: column 0.6, relative 2, absolute 0.8. Second: no column, relative 8.
    expected = np.sqrt([0.6**2 + 2**2 + 0.8**2, 8**2 + 0.8**2]) / 2
    np.testing.assert_allclose(
        changed.standard_uncertainty(), expected, rtol=2e-15, atol=0
    )
    column = dataclasses.replace(sample, uncertainty=Uncertainty(column=True))
    assert column.standard_uncertainty()[0] == 0.3
    assert np.isnan(column.standard_uncertainty()[1])
    assert np.isnan(sample.standard_uncertainty()).all()


def test_unknown_confidence_is_not_silently_treated_as_one_sigma():
    with pytest.raises(DataError, match="unknown uncertainty level"):
        Uncertainty(level="roughly")
    assert Uncertainty(level="not stated").to_standard == 1
    assert Uncertainty(level="bound").to_standard == 1 / math.sqrt(3)


def test_dataset_accessors_keep_printed_cells_and_static_limit(sample):
    assert sample.values == ("10.00", "20.0")
    np.testing.assert_array_equal(sample.abscissa_values(), [500, np.inf])
    # 500 nm corresponds to 2*pi*c/lambda; the static wavelength has zero frequency.
    np.testing.assert_allclose(
        sample.angular_frequency(), [2 * math.pi * C / 500e-9, 0], rtol=2e-15, atol=0
    )
    with pytest.raises(DataError, match="no column"):
        sample.column_values("missing")
    bad = dataclasses.replace(sample, values=("not numeric", "20"))
    with pytest.raises(DataError, match="not a number"):
        bad.scaled_values()


@pytest.mark.parametrize("abscissa", ["label", "none"])
def test_labels_and_standalone_values_are_not_numeric_axes(sample, abscissa):
    labelled = dataclasses.replace(
        sample,
        quantity="ionization_energy",
        unit="eV",
        abscissa=abscissa,
        abscissa_unit="",
        abscissae=("", ""),
        parameters={"kind": "adiabatic"},
    )
    with pytest.raises(DataError, match="has no values"):
        labelled.abscissa_values()


def test_temperatures_convert_to_kelvin_and_cannot_be_not_stated(sample):
    temperature = dataclasses.replace(
        sample,
        quantity="molar_kerr_constant",
        unit="cm3/mol",
        abscissa="temperature",
        abscissa_unit="C",
        abscissae=("0", "20"),
    )
    np.testing.assert_array_equal(temperature.abscissa_values(), [273.15, 293.15])
    with pytest.raises(DataError, match="only a wavelength"):
        dataclasses.replace(temperature, abscissae=("not stated", "20"))


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        ("columns = []", "columns must be a table"),
        ("uncertainty = []", "uncertainty must be a table"),
        ("parameters = []", "parameters must be a table"),
        ('density = "liquid"', "unknown density convention"),
        ("excluded_rows = []", "excluded_rows must map"),
    ],
)
def test_record_tables_are_validated_at_the_loader(tmp_path, extra, message):
    with pytest.raises(DataError, match=message):
        load_source(_source(tmp_path, HEADER + ENTRY + extra))


def test_source_ids_and_dataset_ids_are_unambiguous(tmp_path):
    with pytest.raises(DataError, match="not the file name"):
        load_source(_source(tmp_path, 'id = "different"\n' + HEADER + ENTRY))
    with pytest.raises(DataError, match="duplicate dataset ids"):
        load_source(_source(tmp_path, HEADER + ENTRY + ENTRY))
    with pytest.raises(DataError, match="a point is"):
        load_source(
            _source(
                tmp_path,
                (HEADER + ENTRY).replace('"10.00", "0.3"', '"10.00", "0.3", "extra"'),
            )
        )


def _erratum(**changes):
    fields = {
        "dataset": "sample",
        "row": "a",
        "column": "value",
        "printed": "10.00",
        "corrected": "12.0",
        "reason": "synthetic typo",
        "evidence": "test fixture",
    }
    fields.update(changes)
    return "\n[[errata]]\n" + "\n".join(f'{k} = "{v}"' for k, v in fields.items())


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"dataset": "absent"}, "erratum for unknown dataset"),
        ({"row": "absent"}, "erratum for missing row"),
        ({"column": "absent"}, "erratum for unknown column"),
        ({"column": "row"}, "erratum for unknown column"),
        ({"printed": "10"}, "does not match the extracted"),
        ({"row": "", "column": "scale", "printed": "2"}, "dataset-wide correction"),
        ({"row": "", "column": "value", "printed": "1"}, "dataset-wide correction"),
    ],
)
def test_errata_must_identify_the_exact_printed_target(tmp_path, changes, message):
    with pytest.raises(DataError, match=message):
        load_source(_source(tmp_path, HEADER + ENTRY + _erratum(**changes)))


def test_file_errata_preserve_the_original_and_apply_before_exclusion(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(records, "VALUES_DIRECTORY", tmp_path)
    csv = tmp_path / "refractive" / "synthetic.csv"
    csv.parent.mkdir()
    printed = (
        "row,x,value,sigma,printed\na,500,10.00,0.3,original\nb,inf,20.0,,original\n"
    )
    csv.write_text(printed, encoding="utf-8")
    entry = ENTRY[: ENTRY.index("points =")] + 'file = "synthetic.csv"\n'
    text = (
        HEADER + entry + 'excluded_rows = { b = "synthetic exclusion" }\n' + _erratum()
    )
    path = _source(tmp_path, text)
    source = load_source(path)
    assert source.datasets[0].values == ("12.0",)
    assert source.datasets[0].rows == ("a",)
    assert source.errata[0].printed == "10.00"
    assert csv.read_text(encoding="utf-8") == printed
    assert path.read_text(encoding="utf-8") == text
    # The missing optional sigma column means no stated per-row uncertainty.
    csv.write_text("row,x,value\na,500,10.00\nb,inf,20.0\n", encoding="utf-8")
    assert load_source(path).datasets[0].sigmas == ("",)
    csv.write_text("row,value\na,10.00\n", encoding="utf-8")
    with pytest.raises(DataError, match="lacks columns"):
        load_source(path)
    csv.unlink()
    with pytest.raises(DataError, match=r"data file .* is missing"):
        load_source(path)


def test_data_files_need_a_known_quantity():
    with pytest.raises(DataError, match="unknown quantity"):
        data_file("imaginary", "test.csv")


@pytest.mark.parametrize(
    ("formula", "message"),
    [
        ("coefficients = []", "coefficients must be a table"),
        ('coefficients = { A = "2" }\ncorrections = {}', "corrections a list"),
        (
            'coefficients = { A = "2" }\ncorrections = [{ term = "A" }]',
            "correction lacks",
        ),
        (
            'coefficients = { A = "2" }\n'
            'corrections = [{ term = "B", printed = "2", corrected = "3",'
            ' reason = "test", evidence = "test" }]',
            "missing coefficient",
        ),
        (
            'coefficients = { A = "2" }\n'
            'corrections = [{ term = "A", printed = "2.0", corrected = "3",'
            ' reason = "test", evidence = "test" }]',
            "not the recorded",
        ),
    ],
)
def test_formula_corrections_require_evidence_and_exact_printed_values(
    tmp_path, formula, message
):
    text = HEADER + '\n[[formulas]]\nid = "synthetic"\n' + formula
    with pytest.raises(DataError, match=message):
        load_source(_source(tmp_path, text))


def test_quantity_definitions_reject_unknown_families_and_axes():
    with pytest.raises(ValueError, match="unknown family"):
        Quantity("synthetic", "imaginary")
    with pytest.raises(ValueError, match="unknown abscissae"):
        Quantity("synthetic", "refractive", abscissae=frozenset({"distance"}))

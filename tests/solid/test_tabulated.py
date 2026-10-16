# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The stitched table models (TC-544), tested on curated silica tables."""

import math

import numpy as np
import pytest

from optmat.constants import C
from optmat.datasets import datasets
from optmat.exceptions import DataError, ModelRangeWarning
from optmat.solid.models.tabulated import load_table

# Philipp's n and k of silica glass from Lamy's reflectance (49.6-108.8 nm), the one
# curated pair of n and k tables on shared abscissae; any pair exercises the model.
_N = "philipp_1985/philipp_lamy_n"
_K = "philipp_1985/philipp_lamy_k"


def _config(tmp_path, segments):
    path = tmp_path / "SiO2.toml"
    body = 'material = "SiO2"\n'
    for real, imaginary, low, high in segments:
        body += (
            f'[[segments]]\nreal = "{real}"\nimaginary = "{imaginary}"\n'
            f'range_nm = [{low}, {high}]\nreason = "test"\n'
        )
    path.write_text(body)
    return path


def _dataset(key):
    return next(d for d in datasets() if d.key == key)


def test_the_table_passes_through_its_points(tmp_path):
    model = load_table(_config(tmp_path, [(_N, _K, 40, 120)]))
    n, k = _dataset(_N), _dataset(_K)
    index = model.refractive_index(n.vacuum_wavelength())
    # The interpolant passes through its knots: only rounding in exp(log) remains.
    np.testing.assert_allclose(np.real(index), n.scaled_values(), rtol=1e-12)
    np.testing.assert_allclose(np.imag(index), k.scaled_values(), rtol=1e-12)
    assert model.name == "table/SiO2"
    assert model.kind == "tabulated"


def test_the_derivative_is_that_of_the_interpolant(tmp_path):
    model = load_table(_config(tmp_path, [(_N, _K, 40, 120)]))
    wavelength = np.array([60e-9, 75.3e-9, 101e-9])
    omega = 2 * math.pi * C / wavelength
    h = 1e-7 * omega
    high = model.permittivity(2 * math.pi * C / (omega + h))
    low = model.permittivity(2 * math.pi * C / (omega - h))
    # A central difference of step 1e-7 omega of a C1 piecewise cubic: 1e-6 relative.
    np.testing.assert_allclose(
        model.dpermittivity_domega(wavelength), (high - low) / (2 * h), rtol=1e-6
    )


def test_a_range_selects_points_and_the_ends_are_held(tmp_path):
    model = load_table(_config(tmp_path, [(_N, _K, 60, 90)]))
    assert 60e-9 <= model.wavelength_range[0] < model.wavelength_range[1] <= 90e-9
    with pytest.warns(ModelRangeWarning):
        outside = model.refractive_index(np.array([30e-9, 200e-9]))
    ends = model.refractive_index(np.array(model.wavelength_range))
    np.testing.assert_allclose(outside, ends, rtol=1e-12)


def test_overlapping_ranges_and_bad_pairs_are_refused(tmp_path):
    with pytest.raises(DataError, match="overlap"):
        load_table(_config(tmp_path, [(_N, _K, 40, 80), (_N, _K, 70, 120)]))
    with pytest.raises(DataError, match="cannot pair"):
        load_table(_config(tmp_path, [(_N, _N, 40, 120)]))
    with pytest.raises(DataError, match="no dataset"):
        load_table(_config(tmp_path, [(_N, "philipp_1985/none", 40, 120)]))


def test_a_decrement_pairs_with_k_on_the_wavelengths_they_share(tmp_path):
    # The DESY report prints n to 120 eV and 1 - n above it, and k throughout.
    path = tmp_path / "Al.toml"
    path.write_text(
        'material = "Al"\n[[segments]]\nreal = "hagemann_1974/delta_Al"\n'
        'imaginary = "hagemann_1974/k_Al"\nrange_nm = [1, 11]\nreason = "test"\n'
    )
    model = load_table(path)
    delta, k = _dataset("hagemann_1974/delta_Al"), _dataset("hagemann_1974/k_Al")
    wavelength = delta.vacuum_wavelength()
    inside = (wavelength >= 1e-9) & (wavelength <= 11e-9)
    np.testing.assert_allclose(model.knots, np.sort(wavelength[inside]), rtol=1e-12)
    shared = np.isin(k.vacuum_wavelength(), wavelength[inside])
    assert len(model.knots) == np.count_nonzero(shared)
    index = model.refractive_index(wavelength[inside])
    np.testing.assert_allclose(
        np.real(index), 1 - delta.scaled_values()[inside], rtol=1e-12
    )


def test_henke_segment_follows_the_x_ray_limit(tmp_path):
    # Henke's f1, f2 with Rakic's atom density against Rakic's aluminium at 1 keV,
    # where the independent-atom limit holds: the two agree to a few per cent in k.
    path = tmp_path / "Al.toml"
    path.write_text(
        'material = "Al"\n[[segments]]\nhenke = "al"\n'
        'atom_density = "rakic_1995/atom_density_Al"\nrange_nm = [0.5, 3]\n'
        'reason = "test"\n'
    )
    model = load_table(path)
    assert "Henke" in model.citation
    rakic_k = _dataset("rakic_1995/k_Al")
    at = np.argmin(np.abs(rakic_k.abscissa_values() - 1000.0))
    wavelength = rakic_k.vacuum_wavelength()[at]
    henke = model.refractive_index(np.array([wavelength]))[0]
    assert np.imag(henke) == pytest.approx(rakic_k.scaled_values()[at], rel=0.05)
    with pytest.raises(DataError, match="atom_density"):
        path.write_text(
            'material = "Al"\n[[segments]]\nhenke = "al"\nrange_nm = [0.5, 3]\n'
            'reason = "test"\n'
        )
        load_table(path)

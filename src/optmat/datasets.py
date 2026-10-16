# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Curated datasets: literature values with their provenance and conditions.

A *source* is one paper (or one compilation, theory table or database); its record,
``data/sources/<source>.toml``, holds the citation, where the values are in the paper
library, how they were measured or computed, and one entry per *dataset*: a set of
values of one quantity for one material. A dataset's values are either extracted from
a library table into ``data/values/<family>/<material>/<source>__<dataset>.csv`` (the
family is the quantity's; :mod:`optmat.vocabulary`) by that source's script under
``tools/extract/``, or, when they exist only in the text, listed in the record with
their page.

The files keep the values as printed. Everything that turns them into physics happens
here, once:

- units and printed scale factors are applied;
- wavelengths (and wavenumbers) in standard air become vacuum wavelengths;
- corrections of printed values (errata) are applied, each checked against the value
  it corrects;
- quantities that determine the polarizability (a refractivity at a stated state, a
  polarizability, a molar refractivity) are converted to the per-particle polarizability
  through Lorentz-Lorenz with the number density the source meant.

Quantities that constrain dispersion only indirectly (phase-matching density ratios,
group delays, oscillator strengths, Cauchy moments) are loaded and kept, but have no
polarizability of their own; the fitting layer gives them forward models.
"""

import csv
import functools
import math
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final, Literal

import numpy as np
from jaxtyping import Float

from optmat.conditions import (
    ONE_ATM_0C,
    ONE_ATM_15C,
    ONE_ATM_20C,
    ONE_BAR_0C,
    State,
    vacuum_wavelength_from_standard_air,
)
from optmat.constants import (
    AU_ANGULAR_FREQUENCY,
    AU_FIELD,
    HBAR,
    K_B,
    N_A,
    STD_ATMOSPHERE,
    ZERO_CELSIUS,
    C,
)
from optmat.exceptions import DataError
from optmat.materials import as_material, materials
from optmat.vocabulary import (
    ABSCISSAE,
    CONVENTIONS,
    PARTS,
    POLARIZABILITY_CORRECTIONS,
    PROCESSES,
    QUANTITIES,
)

__all__ = [
    "DATA_DIRECTORY",
    "NOT_STATED",
    "VALUES_DIRECTORY",
    "Dataset",
    "Erratum",
    "Source",
    "Uncertainty",
    "data_file",
    "datasets",
    "load_source",
    "parse_state",
    "sources",
]

DATA_DIRECTORY: Final = Path(__file__).parent / "data"
"""The curated data shipped with the package."""

VALUES_DIRECTORY: Final = DATA_DIRECTORY / "values"
"""The extracted data files, ``<family>/<material>/<source>__<dataset>.csv``."""

_STATES: Final[dict[str, State]] = {
    "1atm_0C": ONE_ATM_0C,
    "1bar_0C": ONE_BAR_0C,
    "1atm_15C": ONE_ATM_15C,
    "1atm_20C": ONE_ATM_20C,
}


# Factors from a stated confidence level to one standard deviation, for a normal law.
_LEVELS: Final[dict[str, float]] = {
    "1 sigma": 1.0,
    "rms": 1.0,
    "standard": 1.0,
    "probable error": 1.482602218505602,  # 0.6745 sigma is the probable error
    "90 % confidence": 1 / 1.6448536269514722,
    "95 % confidence": 1 / 1.959963984540054,
    "99 % confidence": 1 / 2.5758293035489004,
    "2 sigma": 0.5,
    "3 sigma": 1 / 3,
    "bound": 1 / math.sqrt(3),  # a stated bound read as a uniform distribution
}

_ABSCISSA_SCALES: Final[dict[str, float]] = {
    "1": 1.0,
    "Pa": 1.0,
    "bar": 1e5,
    "atm": STD_ATMOSPHERE,
    "torr": STD_ATMOSPHERE / 760,
    "kPa": 1e3,
    "MPa": 1e6,
    "psia": 0.45359237 * 9.80665 / 0.0254**2,  # lbf/in^2, from the exact pound and inch
    "m-3": 1.0,
    "cm-3": 1e6,
    "amagat": STD_ATMOSPHERE / (K_B * ZERO_CELSIUS),
    "mol/m3": N_A,
    "mol/cm3": N_A * 1e6,
    "mol/L": N_A * 1e3,
    "W/m2": 1.0,
    "W/cm2": 1e4,
    "V/m": 1.0,
    "V/cm": 100.0,
    "au": AU_FIELD,
    "K": 1.0,
    "s": 1.0,
    "ns": 1e-9,
    "ps": 1e-12,
    "fs": 1e-15,
}
"""Factors from an abscissa unit to SI (a temperature in C is offset instead)."""

# The closed vocabularies of a solid dataset's parameters (TC-319): its form, whether
# the value is a model input the source assumes rather than a result, the logarithm
# base of an absorption coefficient and the site of laser damage.
_SOLID_PARAMETERS: Final[dict[str, frozenset[object]]] = {
    "form": frozenset({"bulk", "film", "fibre"}),
    "assumed": frozenset({True, False}),
    "base": frozenset({"e", "10"}),
    "site": frozenset({"bulk", "surface"}),
    # What a solid's index is relative to: vacuum (absolute) or the room air of the
    # measurement (TC-398).
    "index_medium": frozenset({"vacuum", "air", "not stated"}),
}

NOT_STATED: Final = "not stated"
"""The abscissa cell of a wavelength the source does not give."""

Density = Literal["real", "ideal", "unknown"]


@dataclass(frozen=True)
class Uncertainty:
    """How the uncertainty of a dataset's values is stated.

    Attributes
    ----------
    relative : float or None
        Per-point relative uncertainty, at ``level``.
    absolute : float or None
        Per-point absolute uncertainty in the printed unit (before ``scale``), at
        ``level``.
    column : bool
        True if the data file carries a per-point ``sigma`` column in the printed unit.
    level : str
        The confidence level as stated ("1 sigma", "probable error", "90 % confidence",
        "95 % confidence", "99 % confidence", "2 sigma", "3 sigma", "rms", "bound"),
        or "not stated".
    systematic_relative : float or None
        A stated relative systematic uncertainty common to every point (a scale error),
        as one standard deviation; the prior of the dataset's scale in a fit.
    relative_rows : dict of str to float
        Per-row relative uncertainties that replace ``relative`` for those rows (a
        source that states a larger error near a resonance, say).
    """

    relative: float | None = None
    absolute: float | None = None
    column: bool = False
    level: str = "not stated"
    systematic_relative: float | None = None
    relative_rows: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Reject an unknown confidence level."""
        if self.level != "not stated" and self.level not in _LEVELS:
            msg = f"unknown uncertainty level {self.level!r}"
            raise DataError(msg)

    @property
    def stated(self) -> bool:
        """Whether any per-point uncertainty is given.

        Returns
        -------
        bool
            True if a column, a relative or an absolute uncertainty is stated.
        """
        return (
            self.column
            or self.relative is not None
            or self.absolute is not None
            or bool(self.relative_rows)
        )

    @property
    def to_standard(self) -> float:
        """The factor from the stated level to one standard deviation.

        Returns
        -------
        float
            1 when the level is one sigma or not stated.
        """
        return _LEVELS.get(self.level, 1.0)


@dataclass(frozen=True)
class Erratum:
    """A correction of one printed value, with its justification.

    Attributes
    ----------
    dataset : str
        The dataset id within the source.
    row : str
        The row id in the data file; "" for a correction of the whole dataset.
    column : str
        "x" (the abscissa), "value", "sigma" or an extra column of the file; or, with
        ``row = ""``, "scale" (a misprinted column factor).
    printed : str
        The value as printed (and as extracted or recorded); checked on loading.
    corrected : str
        The corrected value, in the printed unit.
    reason : str
        Why the printed value is wrong.
    evidence : str
        How the correction was established.
    """

    dataset: str
    row: str
    column: str
    printed: str
    corrected: str
    reason: str
    evidence: str


@dataclass(frozen=True)
class Dataset:
    """Values of one quantity for one material from one source, as printed.

    Attributes
    ----------
    source : str
        The source id.
    id : str
        The dataset id, unique within its source.
    material : str
        The material's id (:mod:`optmat.materials`).
    quantity : str
        What the values are: "refractivity" (``n - 1``), "polarizability",
        "molar_refractivity", "refractivity_per_density" (``(n - 1) / N``), or a
        quantity with no polarizability of its own (e.g. "phase_match_density_ratio",
        "group_delay_dispersion", "oscillator_strength", "cauchy_moment").
    unit : str
        The unit of the values after ``scale``: "1" for dimensionless, a polarizability
        unit ("au", "angstrom3", "cm3", "SI"), "cm3/mol", "cm3" (for a refractivity per
        density) or a quantity-specific unit.
    scale : float
        The printed values times ``scale`` are in ``unit``.
    abscissa : str
        What the ``x`` column is (:data:`~optmat.vocabulary.ABSCISSAE`): "wavelength"
        (a wavelength, wavenumber or photon energy; "inf" is the static limit and "not
        stated" a wavelength the source does not give), "order", "temperature",
        "pressure", "density", "intensity", "duration", "field", "label" or "none".
    abscissa_unit : str
        The unit of the ``x`` column: for a wavelength "m", "nm", "um", "angstrom", a
        wavenumber ("cm-1", "um-1") or a photon energy ("eV", "hartree", "rydberg");
        "" for a label or none.
    wavelength_medium : str
        "vacuum", or "standard_air" (15 degrees C, 101 325 Pa, Edlen 1966) for
        wavelengths and wavenumbers measured in air.
    state : State or None
        The state the values refer to or are reduced to; None for per-particle
        quantities.
    density : str
        How the source related the values to its state: "real" (real-gas density),
        "ideal" (ideal-gas density) or "unknown".
    uncertainty : Uncertainty
        How the uncertainty is stated.
    use : str
        "primary" (enters fits), "reproduced" (a copy of another source's values;
        never fitted beside the original), "check" (compared, not fitted) or "excluded".
    reproduced_from : tuple of str
        Source ids (or citations, for sources outside the library) of reproduced values.
    table : str
        The library table id, or "text" with the page in ``notes``, or, for values
        read from a figure, the library figure id.
    notes : str
        Anything else a reader needs.
    parameters : dict
        Quantity-specific metadata (the reference material of a ratio, the harmonic of a
        phase-matching measurement, the length a group delay refers to, ...); see the
        data-format page of the manual. ``beyond_resonance = true`` marks values on the
        far side of a resonance, where n - 1 may be negative.
    excluded_rows : dict of str to str
        Rows of the file left out of the dataset, each with the reason (the paper's
        own exclusions, points marked extrapolated, ...).
    rows : tuple of str
        Row ids.
    abscissae : tuple of str
        The ``x`` column, as printed (after errata).
    values : tuple of str
        Value column, as printed (after errata).
    sigmas : tuple of str
        Per-point uncertainty column, as printed, or empty strings.
    columns : dict of str to tuple of str
        Extra per-row columns, as printed (after errata): ``x_upper``, the upper edge
        of an interval, in the abscissa's unit; or a per-row condition in the unit of
        ``column_units``.
    column_units : dict of str to str
        The unit of each extra column (``x_upper``'s is the abscissa's).
    digitised : str
        "" for values printed as numbers; "vector" for values read from a figure's
        vector drawing in the PDF (the coordinates of its markers and lines, the axes
        calibrated from the drawing's own ticks); "raster" for values read from a
        figure's image. Figure data are used only where the source prints no number;
        they state their digitising uncertainty, and raster data never enter a fit.
    """

    source: str
    id: str
    material: str
    quantity: str
    unit: str
    scale: float
    abscissa: str
    abscissa_unit: str
    wavelength_medium: str
    state: State | None
    density: Density
    uncertainty: Uncertainty
    use: str
    reproduced_from: tuple[str, ...]
    table: str
    notes: str
    parameters: Mapping[str, object] = field(repr=False)
    excluded_rows: Mapping[str, str] = field(repr=False)
    rows: tuple[str, ...] = field(repr=False)
    abscissae: tuple[str, ...] = field(repr=False)
    values: tuple[str, ...] = field(repr=False)
    sigmas: tuple[str, ...] = field(repr=False)
    columns: Mapping[str, tuple[str, ...]] = field(default_factory=dict, repr=False)
    column_units: Mapping[str, str] = field(default_factory=dict, repr=False)
    digitised: str = ""

    def __post_init__(self) -> None:
        """Check that the columns agree and the metadata are meaningful."""
        lengths = {len(self.rows), len(self.abscissae), len(self.values)}
        lengths |= {len(column) for column in self.columns.values()}
        if len(lengths) != 1 or len(self.sigmas) != len(self.rows):
            msg = f"{self.key}: columns of unequal length"
            raise DataError(msg)
        if len(set(self.rows)) != len(self.rows):
            msg = f"{self.key}: duplicate row ids"
            raise DataError(msg)
        if self.wavelength_medium not in {"vacuum", "standard_air"}:
            msg = f"{self.key}: unknown wavelength medium {self.wavelength_medium!r}"
            raise DataError(msg)
        if self.use not in {"primary", "reproduced", "check", "excluded"}:
            msg = f"{self.key}: unknown use {self.use!r}"
            raise DataError(msg)
        if self.use == "reproduced" and not self.reproduced_from:
            msg = f"{self.key}: reproduced values must name their source"
            raise DataError(msg)
        self._check_digitised()
        self._check_quantity()
        self._check_abscissa()
        self._check_parameters()
        unknown_rows = set(self.uncertainty.relative_rows) - set(self.rows)
        if unknown_rows:
            msg = (
                f"{self.key}: uncertainty given for unknown rows {sorted(unknown_rows)}"
            )
            raise DataError(msg)

    def _check_digitised(self) -> None:
        """Check the rules of values read from a figure (TC-699).

        Raises
        ------
        DataError
            If the kind is unknown, the figure is not named, the digitising uncertainty
            is not stated, or raster data are to be fitted.
        """
        if self.digitised not in {"", "vector", "raster"}:
            msg = f"{self.key}: unknown digitisation {self.digitised!r}"
            raise DataError(msg)
        if not self.digitised:
            return
        if not self.table.startswith("fig_"):
            msg = f"{self.key}: digitised values must name the library figure"
            raise DataError(msg)
        if not self.uncertainty.stated:
            msg = f"{self.key}: digitised values must state their uncertainty"
            raise DataError(msg)
        if self.digitised == "raster" and self.use == "primary":
            msg = f"{self.key}: values read from a figure's image are never fitted"
            raise DataError(msg)

    def _check_quantity(self) -> None:
        """Check the quantity, its unit, its state and its extra columns.

        Raises
        ------
        DataError
            If any of them does not match the vocabulary.
        """
        spec = QUANTITIES.get(self.quantity)
        if spec is None:
            msg = f"{self.key}: unknown quantity {self.quantity!r}"
            raise DataError(msg)
        if spec.units is not None and self.unit not in spec.units:
            msg = (
                f"{self.key}: unit {self.unit!r} is not one of a {self.quantity}'s "
                f"{sorted(spec.units)}"
            )
            raise DataError(msg)
        # A state (pressure and temperature) sets a gas's density; a solid has none, and
        # records its temperature, where stated, as parameters.temperature_k (TC-318).
        solid = materials()[self.material].is_solid
        if solid and (self.state is not None or self.density != "unknown"):
            msg = f"{self.key}: a solid has no gas state or density convention"
            raise DataError(msg)
        if spec.needs_state and self.state is None and not solid:
            msg = f"{self.key}: a {self.quantity} needs the state it refers to"
            raise DataError(msg)
        missing = set(spec.columns) - set(self.columns)
        if missing:
            msg = f"{self.key}: a {self.quantity} needs the columns {sorted(missing)}"
            raise DataError(msg)
        if set(self.column_units) != set(self.columns):
            msg = f"{self.key}: every extra column needs a unit"
            raise DataError(msg)
        if self.column_units.get("x_upper", self.abscissa_unit) != self.abscissa_unit:
            msg = f"{self.key}: x_upper is in the abscissa's unit"
            raise DataError(msg)

    def _check_abscissa(self) -> None:
        """Check the abscissa, its unit, its medium and its cells.

        Raises
        ------
        DataError
            If the abscissa does not suit the quantity or its cells are malformed.
        """
        spec = QUANTITIES[self.quantity]
        if self.abscissa not in spec.abscissae:
            msg = (
                f"{self.key}: a {self.quantity} is tabulated against "
                f"{sorted(spec.abscissae)}, not {self.abscissa!r}"
            )
            raise DataError(msg)
        if self.abscissa_unit not in ABSCISSAE[self.abscissa]:
            msg = (
                f"{self.key}: unknown {self.abscissa} unit {self.abscissa_unit!r}; "
                f"known: {sorted(ABSCISSAE[self.abscissa])}"
            )
            raise DataError(msg)
        if self.wavelength_medium == "standard_air" and self.abscissa_unit not in {
            "m",
            "nm",
            "um",
            "angstrom",
            "cm-1",
            "um-1",
        }:
            msg = f"{self.key}: only wavelengths and wavenumbers can be given in air"
            raise DataError(msg)
        if self.abscissa == "none" and any(self.abscissae):
            msg = f"{self.key}: an abscissa of none has empty cells"
            raise DataError(msg)
        if NOT_STATED in self.abscissae and self.abscissa != "wavelength":
            msg = f"{self.key}: only a wavelength can be {NOT_STATED!r}"
            raise DataError(msg)

    def _check_parameters(self) -> None:
        """Check the parameters the quantity requires and the vocabulary they use.

        Raises
        ------
        DataError
            If a required parameter is missing, or a process, part, correction,
            coefficient or convention is unknown.
        """
        required = QUANTITIES[self.quantity].required
        missing = [name for name in required if name not in self.parameters]
        if missing:
            msg = f"{self.key}: a {self.quantity} must state parameters {missing}"
            raise DataError(msg)
        allowed = {
            "process": PROCESSES,
            "convention": frozenset(CONVENTIONS),
            "coefficient": frozenset({"self", "probe"}),
        }
        for name, values in allowed.items():
            if name in required and self.parameters[name] not in values:
                msg = f"{self.key}: unknown {name} {self.parameters[name]!r}"
                raise DataError(msg)
        if "part" in required:
            parts = (
                POLARIZABILITY_CORRECTIONS
                if self.quantity == "polarizability_correction"
                else PARTS
            )
            if self.parameters["part"] not in parts:
                msg = f"{self.key}: unknown part {self.parameters['part']!r}"
                raise DataError(msg)
        axes = materials()[self.material].axes
        if "axis" in self.parameters and self.parameters["axis"] not in axes:
            msg = f"{self.key}: axis {self.parameters['axis']!r} is not one of "
            msg += f"{self.material}'s {list(axes)}"
            raise DataError(msg)
        for name, values in _SOLID_PARAMETERS.items():
            if name in self.parameters and self.parameters[name] not in values:
                msg = f"{self.key}: {name} {self.parameters[name]!r} is not one of "
                msg += f"{sorted(values, key=str)}"
                raise DataError(msg)

    @property
    def key(self) -> str:
        """The dataset's global id.

        Returns
        -------
        str
            ``source/dataset``.
        """
        return f"{self.source}/{self.id}"

    def __len__(self) -> int:
        """Return the number of values.

        Returns
        -------
        int
            The number of rows.
        """
        return len(self.rows)

    def vacuum_wavelength(self) -> Float[np.ndarray, " n"]:
        """Return the vacuum wavelengths.

        Returns
        -------
        numpy.ndarray
            Vacuum wavelength of each value, m; ``inf`` for a static value, NaN where
            the source does not state it.

        Raises
        ------
        DataError
            If the abscissa is not a wavelength.
        """
        if self.abscissa != "wavelength":
            msg = f"{self.key}: the abscissa is a {self.abscissa}, not a wavelength"
            raise DataError(msg)
        return _vacuum_wavelength(
            np.array(
                [math.nan if x == NOT_STATED else _number(x) for x in self.abscissae]
            ),
            self.abscissa_unit,
            self.wavelength_medium,
        )

    def abscissa_values(self) -> Float[np.ndarray, " n"]:
        """Return the abscissa as numbers in SI (a wavelength as printed).

        Returns
        -------
        numpy.ndarray
            A temperature in K, a pressure in Pa, a density in m^-3, an intensity in
            W/m^2, a field in V/m, an order as printed; a wavelength, wavenumber or
            energy in its printed unit (use :meth:`vacuum_wavelength` for the vacuum
            wavelength).

        Raises
        ------
        DataError
            If the abscissa is a label or none.
        """
        if self.abscissa in {"label", "none"}:
            msg = f"{self.key}: an abscissa of {self.abscissa} has no values"
            raise DataError(msg)
        numbers = np.array(
            [math.nan if x == NOT_STATED else _number(x) for x in self.abscissae]
        )
        if self.abscissa == "wavelength":
            return numbers
        if self.abscissa_unit == "C":
            return numbers + ZERO_CELSIUS
        return numbers * _ABSCISSA_SCALES[self.abscissa_unit]

    def column_values(self, name: str) -> Float[np.ndarray, " n"]:
        """Return an extra column as numbers, as printed.

        Parameters
        ----------
        name : str
            The column (``x_upper``, or a per-row condition).

        Returns
        -------
        numpy.ndarray
            The column, in its unit (``x_upper`` in the abscissa's).

        Raises
        ------
        DataError
            If the dataset has no such column.
        """
        if name not in self.columns:
            msg = f"{self.key}: no column {name!r}"
            raise DataError(msg)
        return np.array([_number(x) for x in self.columns[name]])

    def angular_frequency(self) -> Float[np.ndarray, " n"]:
        """Return the angular frequencies ``2 pi c / lambda_vac``.

        Returns
        -------
        numpy.ndarray
            Angular frequency of each value, rad/s; 0 for a static value.
        """
        return 2 * math.pi * C / self.vacuum_wavelength()

    def scaled_values(self) -> Float[np.ndarray, " n"]:
        """Return the values times ``scale``, in ``unit``.

        Returns
        -------
        numpy.ndarray
            The values.
        """
        return np.array([_number(v) for v in self.values]) * self.scale

    def standard_uncertainty(self) -> Float[np.ndarray, " n"]:
        """Return each value's stated uncertainty as one standard deviation.

        Per-point, relative and absolute statements are combined in quadrature.

        Returns
        -------
        numpy.ndarray
            One-sigma uncertainty in ``unit``; NaN where none is stated.
        """
        values = self.scaled_values()
        spec = self.uncertainty
        if not spec.stated:
            return np.full(len(self), np.nan)
        variance = np.zeros(len(self))
        # A row without any statement that covers it has no uncertainty (NaN), not 0.
        covered = np.zeros(len(self), dtype=bool)
        if spec.column:
            sigmas = np.array(
                [_number(s) if s.strip() else math.nan for s in self.sigmas]
            )
            given = np.isfinite(sigmas)
            variance += np.where(given, sigmas * self.scale, 0.0) ** 2
            covered |= given
        if spec.relative is not None or spec.relative_rows:
            relative = np.array(
                [spec.relative_rows.get(row, spec.relative or 0.0) for row in self.rows]
            )
            variance += (relative * values) ** 2
            covered |= np.array(
                [
                    spec.relative is not None or row in spec.relative_rows
                    for row in self.rows
                ]
            )
        if spec.absolute is not None:
            variance += (spec.absolute * self.scale) ** 2
            covered[:] = True
        return np.where(covered, np.sqrt(variance) * spec.to_standard, math.nan)


@dataclass(frozen=True)
class Source:
    """One source of values: its citation, provenance and datasets.

    Attributes
    ----------
    id : str
        The source id, the record's file stem.
    citation : str
        The citation as it should be printed.
    doi : str
        The DOI, or "" if the work has none.
    library : str
        The paper-library directory, or "" if the work is not in the library.
    identity : str
        The library's bibliographic status of the record ("VALIDATED", ...).
    kind : str
        "measurement", "theory", "compilation", "fit" or "database".
    technique : str
        How the values were obtained.
    notes : str
        Anything else a reader needs.
    datasets : tuple of Dataset
        The datasets, errata applied.
    errata : tuple of Erratum
        The corrections applied.
    formulas : tuple of dict
        The source's published formulas, as recorded (``[[formulas]]``); the ones with
        a ``form`` are built into models by :mod:`optmat.gas.models.formulas`.
    """

    id: str
    citation: str
    doi: str
    library: str
    identity: str
    kind: str
    technique: str
    notes: str
    datasets: tuple[Dataset, ...]
    errata: tuple[Erratum, ...]
    formulas: tuple[Mapping[str, object], ...] = ()


def _number(text: str) -> float:
    """Parse one printed number.

    Parameters
    ----------
    text : str
        The number as printed, or "inf" for the static limit.

    Returns
    -------
    float
        The number.

    Raises
    ------
    DataError
        If the text is not a number.
    """
    text = text.strip()
    if text == "inf":
        return math.inf
    try:
        return float(text)
    except ValueError as exc:
        msg = f"not a number: {text!r}"
        raise DataError(msg) from exc


def _vacuum_wavelength(
    numbers: Float[np.ndarray, " n"], unit: str, medium: str
) -> Float[np.ndarray, " n"]:
    """Convert a wavelength, wavenumber or photon-energy column to vacuum wavelengths.

    Parameters
    ----------
    numbers : numpy.ndarray
        The column's values.
    unit : str
        The column's unit.
    medium : str
        "vacuum" or "standard_air"; only wavelengths and wavenumbers can be given in
        air (a wavenumber in air is the inverse of the air wavelength).

    Returns
    -------
    numpy.ndarray
        Vacuum wavelength in m; ``inf`` for the static limit.

    Notes
    -----
    The caller's Dataset has already validated the unit and wavelength medium.
    """
    lengths = {"m": 1.0, "nm": 1e-9, "um": 1e-6, "angstrom": 1e-10}
    with np.errstate(divide="ignore"):
        if unit in lengths:
            wavelength = numbers * lengths[unit]
        elif unit == "cm-1":
            wavelength = 1e-2 / numbers
        elif unit == "um-1":
            wavelength = 1e-6 / numbers
        elif unit == "eV":
            wavelength = 2 * math.pi * HBAR * C / (numbers * 1.602176634e-19)
        elif unit == "hartree":  # photon energy in atomic units
            wavelength = 2 * math.pi * C / (numbers * AU_ANGULAR_FREQUENCY)
        else:  # rydberg (R_infinity hc, half a hartree)
            wavelength = 2 * math.pi * C / (numbers * AU_ANGULAR_FREQUENCY / 2)
    if medium == "standard_air":
        finite = np.isfinite(wavelength)
        wavelength = wavelength.copy()
        wavelength[finite] = vacuum_wavelength_from_standard_air(wavelength[finite])
    return wavelength


def parse_state(spec: object, where: str) -> State | None:
    """Parse a dataset's state.

    Parameters
    ----------
    spec : object
        A state name, a table with a pressure and a temperature, or None.
    where : str
        The dataset, for error messages.

    Returns
    -------
    State or None
        The state, or None if none is given.

    Raises
    ------
    DataError
        If the state is unknown or incomplete.
    """
    if spec is None:
        return None
    if isinstance(spec, str):
        if spec not in _STATES:
            msg = f"{where}: unknown state {spec!r}; known: {', '.join(_STATES)}"
            raise DataError(msg)
        return _STATES[spec]
    if isinstance(spec, Mapping):
        pressure = spec.get("pressure_pa")
        if pressure is None and "pressure_atm" in spec:
            pressure = spec["pressure_atm"] * STD_ATMOSPHERE
        temperature = spec.get("temperature_k")
        if temperature is None and "temperature_c" in spec:
            temperature = spec["temperature_c"] + ZERO_CELSIUS
        if pressure is None or temperature is None:
            msg = f"{where}: a state needs a pressure and a temperature"
            raise DataError(msg)
        return State(float(pressure), float(temperature))
    msg = f"{where}: cannot read state {spec!r}"
    raise DataError(msg)


_Columns = dict[str, list[str]]
"""A data file's columns as printed: row, x, value, sigma and any extra columns."""


def _read_points(path: Path, extra: Sequence[str], where: str) -> _Columns:
    """Read the columns of an extracted file.

    Parameters
    ----------
    path : pathlib.Path
        The data file.
    extra : sequence of str
        The extra columns the dataset declares.
    where : str
        The dataset, for error messages.

    Returns
    -------
    dict of str to list of str
        The ``row``, ``x``, ``value``, ``sigma`` and extra columns, as printed.

    Raises
    ------
    DataError
        If the file is missing or lacks a column.
    """
    if not path.is_file():
        msg = f"{where}: data file {path} is missing"
        raise DataError(msg)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = {"row", "x", "value", *extra} - set(reader.fieldnames or ())
        if missing:
            msg = f"{where}: {path.name} lacks columns {sorted(missing)}"
            raise DataError(msg)
        records = list(reader)
    columns = {name: [r[name] for r in records] for name in ("row", "x", "value")}
    columns["sigma"] = [r.get("sigma") or "" for r in records]
    for name in extra:
        columns[name] = [r[name] for r in records]
    return columns


def _inline_points(
    points: Sequence[Sequence[object]], extra: Sequence[str], where: str
) -> _Columns:
    """Read the inline points of a source record.

    Parameters
    ----------
    points : sequence of sequence
        ``[row, x, value]``, ``[row, x, value, sigma]``, or, with extra columns,
        ``[row, x, value, sigma, extra...]`` in the declared order (sigma may be "").
    extra : sequence of str
        The extra columns the dataset declares.
    where : str
        The dataset, for error messages.

    Returns
    -------
    dict of str to list of str
        The columns, as text.

    Raises
    ------
    DataError
        If an entry has the wrong length.
    """
    names = ["row", "x", "value", "sigma", *extra]
    columns: _Columns = {name: [] for name in names}
    lengths = {4 + len(extra)} if extra else {3, 4}
    for point in points:
        if len(point) not in lengths:
            msg = f"{where}: a point is [{', '.join(names)}] (sigma optional alone)"
            raise DataError(msg)
        cells = [str(cell) for cell in point] + [""] * (len(names) - len(point))
        for name, cell in zip(names, cells, strict=True):
            columns[name].append(cell)
    return columns


def _apply_errata(columns: _Columns, scale: str, errata: Sequence[Erratum]) -> str:
    """Apply a dataset's errata, checking each printed value.

    Parameters
    ----------
    columns : dict of str to list of str
        The columns, changed in place.
    scale : str
        The dataset's scale as recorded.
    errata : sequence of Erratum
        This dataset's errata.

    Returns
    -------
    str
        The scale, corrected if an erratum corrects it.

    Raises
    ------
    DataError
        If an erratum names a missing row or column, or its printed value is not the
        recorded one.
    """
    rows = columns["row"]
    for erratum in errata:
        if not erratum.row:
            if erratum.column != "scale" or erratum.printed != scale:
                msg = (
                    f"erratum of {erratum.dataset!r}: a dataset-wide correction is of "
                    f"the scale, printed as recorded ({scale!r})"
                )
                raise DataError(msg)
            scale = erratum.corrected
            continue
        if erratum.row not in rows:
            msg = f"erratum for missing row {erratum.row!r} of {erratum.dataset!r}"
            raise DataError(msg)
        if erratum.column not in columns or erratum.column == "row":
            msg = f"erratum for unknown column {erratum.column!r}"
            raise DataError(msg)
        index = rows.index(erratum.row)
        current = columns[erratum.column][index]
        if current != erratum.printed:
            msg = (
                f"erratum for {erratum.dataset!r} row {erratum.row!r}: printed value "
                f"{erratum.printed!r} does not match the extracted {current!r}"
            )
            raise DataError(msg)
        columns[erratum.column][index] = erratum.corrected
    return scale


def data_file(quantity: str, file: str) -> Path:
    """Return where a dataset's file lives.

    Parameters
    ----------
    quantity : str
        The dataset's quantity, whose family is the directory.
    file : str
        The record's ``file``, ``<material>/<source>__<dataset>.csv``.

    Returns
    -------
    pathlib.Path
        ``data/values/<family>/<file>``.

    Raises
    ------
    DataError
        If the quantity is unknown.
    """
    if quantity not in QUANTITIES:
        msg = f"unknown quantity {quantity!r}"
        raise DataError(msg)
    return VALUES_DIRECTORY / QUANTITIES[quantity].family / file


def _material_id(entry: Mapping[str, object], where: str) -> str:
    """Return the dataset's material id, which must be a registered id (not an alias).

    Parameters
    ----------
    entry : mapping
        The dataset's entry in the source record.
    where : str
        The dataset key, for messages.

    Returns
    -------
    str
        The material id.

    Raises
    ------
    DataError
        If the material is missing or is not a registered id.
    """
    name = str(entry.get("material", ""))
    if name not in materials():
        msg = f"{where}: material {name!r} is not a registered id (data/materials.toml)"
        raise DataError(msg)
    return name


def _dataset(
    source_id: str, entry: Mapping[str, object], errata: Sequence[Erratum]
) -> Dataset:
    """Build one dataset from its record entry, reading its file and applying errata.

    Parameters
    ----------
    source_id : str
        The source.
    entry : mapping
        The dataset's entry in the source record.
    errata : sequence of Erratum
        The source's errata.

    Returns
    -------
    Dataset
        The dataset.

    Raises
    ------
    DataError
        If the entry is malformed.
    """
    get = entry.get
    dataset_id = str(entry["id"])
    where = f"{source_id}/{dataset_id}"
    quantity = str(entry["quantity"])
    declared = get("columns", {})
    if not isinstance(declared, Mapping):
        msg = f"{where}: columns must be a table of units"
        raise DataError(msg)
    extra = [str(name) for name in declared]
    points = get("points")
    if points is not None:
        columns = _inline_points(
            points,  # pyright: ignore[reportArgumentType] - TOML array of arrays
            extra,
            where,
        )
    else:
        columns = _read_points(data_file(quantity, str(entry["file"])), extra, where)
    scale = _apply_errata(
        columns,
        str(get("scale", "1")),
        [e for e in errata if e.dataset == dataset_id],
    )
    rows = columns["row"]
    excluded = get("excluded_rows", {})
    if not isinstance(excluded, Mapping) or not set(excluded) <= set(rows):
        msg = f"{where}: excluded_rows must map existing row ids to reasons"
        raise DataError(msg)
    kept = [i for i, row in enumerate(rows) if row not in excluded]
    columns = {name: [column[i] for i in kept] for name, column in columns.items()}
    uncertainty = get("uncertainty", {})
    if not isinstance(uncertainty, Mapping):
        msg = f"{where}: uncertainty must be a table"
        raise DataError(msg)
    density = str(get("density", "unknown"))
    if density not in {"real", "ideal", "unknown"}:
        msg = f"{where}: unknown density convention {density!r}"
        raise DataError(msg)
    reproduced = get("reproduced_from", [])
    parameters = get("parameters", {})
    if not isinstance(parameters, Mapping):
        msg = f"{where}: parameters must be a table"
        raise DataError(msg)
    return Dataset(
        source=source_id,
        id=dataset_id,
        material=_material_id(entry, where),
        quantity=quantity,
        unit=str(get("unit", "1")),
        scale=_number(scale),
        abscissa=str(get("abscissa", "wavelength")),
        abscissa_unit=str(get("abscissa_unit", "")),
        wavelength_medium=str(get("wavelength_medium", "vacuum")),
        state=parse_state(get("state"), where),
        density=density,  # pyright: ignore[reportArgumentType] - checked above
        uncertainty=Uncertainty(**uncertainty),  # pyright: ignore[reportArgumentType]
        use=str(get("use", "primary")),
        reproduced_from=tuple(reproduced),  # pyright: ignore[reportArgumentType]
        table=str(get("table", "text")),
        notes=str(get("notes", "")),
        parameters=dict(parameters),
        excluded_rows={str(k): str(v) for k, v in excluded.items()},
        rows=tuple(columns["row"]),
        abscissae=tuple(columns["x"]),
        values=tuple(columns["value"]),
        sigmas=tuple(columns["sigma"]),
        columns={name: tuple(columns[name]) for name in extra},
        column_units={str(k): str(v) for k, v in declared.items()},
        digitised=str(get("digitised", "")),
    )


def _corrected_formula(entry: Mapping[str, object], where: str) -> Mapping[str, object]:
    """Apply the corrections of a formula's named coefficients.

    A correction with a ``term`` and an ``index`` corrects a coefficient list of a
    formula with a ``form`` and is applied when the model is built
    (:mod:`optmat.gas.models.formulas`); one with a ``term`` alone corrects the named
    entry of ``coefficients`` here.

    Parameters
    ----------
    entry : mapping
        The ``[[formulas]]`` entry.
    where : str
        The formula, for error messages.

    Returns
    -------
    mapping
        The entry, with its named coefficients corrected.

    Raises
    ------
    DataError
        If a correction lacks its justification, names a missing coefficient, or its
        printed value is not the recorded one.
    """
    corrections = entry.get("corrections", [])
    raw = entry.get("coefficients", {})
    if not isinstance(raw, Mapping) or not isinstance(corrections, list):
        msg = f"{where}: coefficients must be a table and corrections a list"
        raise DataError(msg)
    coefficients = {str(name): str(value) for name, value in raw.items()}
    for correction in corrections:
        missing = {"term", "printed", "corrected", "reason", "evidence"}
        missing -= set(correction)
        if missing:
            msg = f"{where}: correction lacks {sorted(missing)}"
            raise DataError(msg)
        if "index" in correction:
            continue
        term = str(correction["term"])
        if term not in coefficients:
            msg = f"{where}: correction of missing coefficient {term!r}"
            raise DataError(msg)
        if coefficients[term] != correction["printed"]:
            msg = (
                f"{where}: correction's printed {term} {correction['printed']!r} is "
                f"not the recorded {coefficients[term]!r}"
            )
            raise DataError(msg)
        coefficients[term] = str(correction["corrected"])
    if not coefficients:
        return entry
    return {**entry, "coefficients": coefficients}


def load_source(path: Path) -> Source:
    """Load one source record and its datasets.

    Parameters
    ----------
    path : pathlib.Path
        The record, ``data/sources/<source>.toml``.

    Returns
    -------
    Source
        The source, with its datasets read and errata applied.

    Raises
    ------
    DataError
        If the record or a data file is malformed or inconsistent.
    """
    with path.open("rb") as handle:
        record = tomllib.load(handle)
    source_id = path.stem
    if record.get("id", source_id) != source_id:
        msg = f"{path.name}: id {record['id']!r} is not the file name"
        raise DataError(msg)
    errata = tuple(Erratum(**entry) for entry in record.get("errata", []))
    dataset_ids = [str(entry["id"]) for entry in record.get("datasets", [])]
    for erratum in errata:
        if erratum.dataset not in dataset_ids:
            msg = f"{source_id}: erratum for unknown dataset {erratum.dataset!r}"
            raise DataError(msg)
    if len(set(dataset_ids)) != len(dataset_ids):
        msg = f"{source_id}: duplicate dataset ids"
        raise DataError(msg)
    return Source(
        id=source_id,
        citation=record["citation"],
        doi=record.get("doi", ""),
        library=record.get("library", ""),
        identity=record.get("identity", ""),
        kind=record["kind"],
        technique=record.get("technique", ""),
        notes=record.get("notes", ""),
        datasets=tuple(
            _dataset(source_id, entry, errata) for entry in record.get("datasets", [])
        ),
        errata=errata,
        formulas=tuple(
            _corrected_formula(entry, f"{source_id}/{entry.get('id', '?')}")
            for entry in record.get("formulas", [])
        ),
    )


@functools.cache
def sources() -> dict[str, Source]:
    """Return every curated source, by id.

    Returns
    -------
    dict of str to Source
        The sources, loaded once per process.

    Examples
    --------
    >>> "huber_1974" in sources()
    True
    """
    return {
        path.stem: load_source(path)
        for path in sorted((DATA_DIRECTORY / "sources").glob("*.toml"))
    }


def datasets(
    material: str | None = None, *, quantity: str | None = None
) -> tuple[Dataset, ...]:
    """Return the curated datasets, optionally for one material or one quantity.

    Parameters
    ----------
    material : str, optional
        Only this material: its id or an alias (a gas species such as ``"He"``, or
        ``"SiO2"``).
    quantity : str, optional
        Only this quantity (e.g. "refractivity").

    Returns
    -------
    tuple of Dataset
        The matching datasets, in source order.

    Examples
    --------
    >>> all(d.material == "He" for d in datasets("He"))
    True
    """
    wanted = None if material is None else as_material(material)
    return tuple(
        dataset
        for source in sources().values()
        for dataset in source.datasets
        if (wanted is None or dataset.material == wanted)
        and (quantity is None or dataset.quantity == quantity)
    )

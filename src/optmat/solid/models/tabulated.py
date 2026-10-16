# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Stitched tables of optical constants as models: a metal's complex index.

A metal has no single model from the ultraviolet to the far infrared: its optical
constants are measured piecewise, each range by a different technique (electron
energy-loss spectra, ellipsometry, reflectance with Kramers-Kronig analysis,
far-infrared reflectance). A table model joins curated datasets, one per wavelength
range, as its configuration ``data/solid/tables/<material>.toml`` chooses them, each
range with its reason (TC-544)::

    material = "Ag"

    [[segments]]
    real = "babar_2015/eps1_Ag"          # refractive_index or permittivity_real
    imaginary = "babar_2015/eps2_Ag"     # extinction_coefficient or ...imaginary
    range_nm = [207.0, 12400.0]
    reason = "..."

The model takes the wavelengths the two datasets of a segment share (an x-ray table may
print n only to some energy and 1 - n above it, as ``refractive_index_decrement``,
while its k runs on). In the extreme
ultraviolet and x-ray range a segment may instead be Henke et al.'s atomic scattering
factors (``henke = "<element>"``, CXRO's file in ``data/external/henke/``) with a
curated ``atom_density`` dataset: n + ik = 1 - r_e lambda^2 N (f1 - i f2) / (2 pi)
(Henke et al. 1993, Eq. 3; their n = 1 - delta - i beta is this package's n + ik with
k = beta), the independent-atom limit, which holds well above the valence excitations
and away from absorption edges (TC-645)::

    [[segments]]
    henke = "ag"
    atom_density = "<source>/<atom_density dataset>"
    range_nm = [1.0, 12.0]
    reason = "..."

The model
interpolates log n and log k against log omega with a monotone piecewise-cubic Hermite
interpolant (PCHIP), which neither overshoots between points nor rings at a join
between sources (a natural cubic spline does both); its first derivative is continuous
and exact for the interpolant, so epsilon and d epsilon / d omega are consistent.
Outside the tabulated range the constants are held at the end points, with a
:class:`~optmat.exceptions.ModelRangeWarning`.
"""

import dataclasses
import functools
import itertools
import math
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np
from jaxtyping import Complex, Float

from optmat import _catalogue
from optmat.constants import CLASSICAL_ELECTRON_RADIUS, E_CHARGE, HBAR, C
from optmat.datasets import DATA_DIRECTORY, Dataset
from optmat.exceptions import DataError, StateError
from optmat.materials import as_material
from optmat.solid import SOLID_DATA_DIRECTORY
from optmat.solid.models._defaults import curated
from optmat.solid.models.base import Model, complex_index

if TYPE_CHECKING:
    from scipy.interpolate import PchipInterpolator

__all__ = [
    "HENKE_DIRECTORY",
    "TABLES_DIRECTORY",
    "Segment",
    "TabulatedModel",
    "load_table",
    "tabulated_models",
]

TABLES_DIRECTORY = SOLID_DATA_DIRECTORY / "tables"
"""The table models' configurations, one per material."""

_PAIRS = {
    "refractive_index": "extinction_coefficient",
    "refractive_index_decrement": "extinction_coefficient",
    "permittivity_real": "permittivity_imaginary",
}
"""The quantity of a segment's real dataset, and the one its imaginary dataset needs."""

HENKE_DIRECTORY = DATA_DIRECTORY / "external" / "henke"
"""CXRO's tables of Henke et al.'s atomic scattering factors, ``<element>.nff``."""

_DENSITY_UNITS = {"m-3": 1.0, "cm-3": 1e6}
"""Factors from an atom density's unit to m^-3."""


@dataclass(frozen=True)
class Segment:
    """One range of a table model and the datasets that supply it.

    Attributes
    ----------
    real : str
        The key of the dataset of n or epsilon_1 (empty for a Henke segment).
    imaginary : str
        The key of the dataset of k or epsilon_2, on the same abscissae (empty for a
        Henke segment).
    range_nm : tuple of float
        The vacuum wavelengths taken from these datasets, nm.
    reason : str
        Why these datasets supply this range.
    henke : str
        The element of a Henke segment (CXRO's file name), or empty.
    atom_density : str
        The key of the ``atom_density`` dataset of a Henke segment, or empty.
    """

    real: str
    imaginary: str
    range_nm: tuple[float, float]
    reason: str
    henke: str = ""
    atom_density: str = ""

    @property
    def label(self) -> str:
        """Return what supplies the segment, for messages.

        Returns
        -------
        str
            ``henke/<element>`` or the key of the real dataset.
        """
        return f"henke/{self.henke}" if self.henke else self.real

    @property
    def cited(self) -> tuple[str, ...]:
        """Return the sources the segment rests on.

        Returns
        -------
        tuple of str
            The source ids.
        """
        if self.henke:
            return ("henke_1993", self.atom_density.split("/")[0])
        return (self.real.split("/")[0],)


@dataclass(frozen=True)
class TabulatedModel(Model):
    """A complex index interpolated from curated tables, joined by range.

    Attributes
    ----------
    segments : tuple of Segment
        The ranges, in increasing wavelength.
    knots : numpy.ndarray
        The tabulated vacuum wavelengths, m, increasing.
    index : numpy.ndarray
        The complex index n + ik at the knots.
    out_of_range : {"hold", "raise"}
        Outside the knots the table holds its end values, with a
        :class:`~optmat.exceptions.ModelRangeWarning` (``"hold"``, the default), or
        refuses with :class:`~optmat.exceptions.StateError` (``"raise"``;
        :meth:`with_out_of_range`).
    """

    segments: tuple[Segment, ...] = ()
    knots: Float[np.ndarray, " n"] = field(
        default_factory=lambda: np.zeros(0), repr=False
    )
    index: Complex[np.ndarray, " n"] = field(
        default_factory=lambda: np.zeros(0, dtype=complex), repr=False
    )
    out_of_range: str = field(default="hold", kw_only=True)

    EVALUATION_VERSION: ClassVar[int] = 1
    INPUTS: ClassVar[tuple[str, ...]] = ("index", "intervals", "knots")
    DERIVED: ClassVar[tuple[str, ...]] = ("segments", "wavelength_range")

    @classmethod
    def from_arrays(
        cls,
        material: str,
        wavelength: Float[np.ndarray, " n"],
        n: Float[np.ndarray, " n"],
        k: Float[np.ndarray, " n"],
        *,
        citation: str = "",
        name: str = "",
    ) -> "TabulatedModel":
        """Return a table model from a caller's arrays.

        Parameters
        ----------
        material : str
            The material, the identity's (a curated id or a caller's name).
        wavelength : numpy.ndarray
            Vacuum wavelengths, m, increasing; the range is their span.
        n : numpy.ndarray
            The real index at them, > 0.
        k : numpy.ndarray
            The extinction coefficient, > 0 (the interpolation is in log k).
        citation : str, optional
            Where the table comes from.
        name : str, optional
            The model's name; ``user/<material>`` by default.

        Returns
        -------
        TabulatedModel
            The model, status ``user``.

        Raises
        ------
        ValueError
            If the arrays differ in length, the wavelengths do not increase, or n or
            k is not positive.

        Examples
        --------
        >>> import numpy as np
        >>> table = TabulatedModel.from_arrays(
        ...     "Au", np.array([500e-9, 600e-9, 700e-9]), np.array([0.9, 0.25, 0.16]),
        ...     np.array([1.9, 3.0, 3.9]))
        >>> table.wavelength_range, table.identity.status
        ((5e-07, 7e-07), 'user')
        """
        knots = np.asarray(wavelength, dtype=float)
        real, imaginary = np.asarray(n, dtype=float), np.asarray(k, dtype=float)
        _check_table(knots, real + 1j * imaginary)
        return cls(
            material=material,
            name=name or f"user/{material}",
            citation=citation,
            wavelength_range=(float(knots[0]), float(knots[-1])),
            kind="user",
            knots=knots,
            index=real + 1j * imaginary,
        )

    def _override(self, changes: Mapping[str, Any]) -> "TabulatedModel":
        """Return the table with some inputs changed; its range follows the knots.

        Parameters
        ----------
        changes : Mapping
            Inputs to new values.

        Returns
        -------
        TabulatedModel
            The rebuilt model, without the curated segments.
        """
        changed: dict[str, Any] = dict(changes)
        knots = np.asarray(changed.get("knots", self.knots), dtype=float)
        index = np.asarray(changed.get("index", self.index), dtype=complex)
        _check_table(knots, index)
        changed |= {"knots": knots, "index": index, "segments": ()}
        if "intervals" not in changed:
            changed["intervals"] = ((float(knots[0]), float(knots[-1])),)
        return super()._override(changed)

    def __post_init__(self) -> None:
        """Check the out-of-range rule.

        Raises
        ------
        DataError
            If it is not ``hold`` or ``raise``.
        """
        super().__post_init__()
        if self.out_of_range not in {"hold", "raise"}:
            msg = f"out_of_range is 'hold' or 'raise', not {self.out_of_range!r}"
            raise DataError(msg)

    def with_out_of_range(self, mode: str) -> "TabulatedModel":
        """Return the table with another rule outside its knots, its status kept.

        Not an override: within the table the numbers are the same.

        Parameters
        ----------
        mode : {"hold", "raise"}
            Hold the end values (with a warning), or refuse.

        Returns
        -------
        TabulatedModel
            The model with that rule.

        Examples
        --------
        >>> from optmat.solid.models.registry import model
        >>> strict = model("table/Ag").with_out_of_range("raise")
        >>> strict.out_of_range, strict.status
        ('raise', 'validated')
        """
        return dataclasses.replace(self, out_of_range=mode)

    @functools.cached_property
    def _interpolants(
        self,
    ) -> "tuple[PchipInterpolator, PchipInterpolator, float, float]":
        """Return the interpolants of log n and log k against log omega.

        Returns
        -------
        tuple
            The interpolants and the lowest and highest log omega.
        """
        # Imported on first use: scipy.interpolate is slow to import.
        from scipy.interpolate import PchipInterpolator  # noqa: PLC0415

        log_omega = np.log(2 * math.pi * C / self.knots)[::-1]
        log_n = np.log(np.real(self.index))[::-1]
        log_k = np.log(np.imag(self.index))[::-1]
        return (
            PchipInterpolator(log_omega, log_n),
            PchipInterpolator(log_omega, log_k),
            float(log_omega[0]),
            float(log_omega[-1]),
        )

    def _index_and_slope(self, omega: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return n + ik and d(n + ik)/d omega, held constant outside the table.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        tuple of numpy.ndarray
            The complex index and its derivative, s/rad.
        """
        log_n, log_k, lowest, highest = self._interpolants
        x = np.log(np.asarray(omega, dtype=float))
        inside = (x >= lowest) & (x <= highest)
        if self.out_of_range == "raise" and np.any(np.isfinite(x) & ~inside):
            low, high = self.knots[0], self.knots[-1]
            msg = (
                f"{self.name} refuses wavelengths outside its table, "
                f"{low * 1e9:.4g}-{high * 1e9:.4g} nm (out_of_range='raise')"
            )
            raise StateError(msg)
        clipped = np.clip(x, lowest, highest)
        n = np.exp(log_n(clipped))
        k = np.exp(log_k(clipped))
        # d/d omega = d/d(log omega) / omega; zero where the end point is held.
        scale = np.where(inside, 1 / np.asarray(omega, dtype=float), 0.0)
        dn = n * log_n(clipped, 1) * scale
        dk = k * log_k(clipped, 1) * scale
        return np.asarray(n + 1j * k), np.asarray(dn + 1j * dk)

    def _permittivity(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return epsilon = (n + ik)^2 from the interpolated index.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The permittivity.
        """
        index, _ = self._index_and_slope(omega)
        return np.asarray(index**2)

    def _dpermittivity_domega(
        self, omega: Float[np.ndarray, "*shape"]
    ) -> Complex[np.ndarray, "*shape"]:
        """Return d epsilon / d omega = 2 (n + ik) d(n + ik)/d omega.

        Parameters
        ----------
        omega : numpy.ndarray
            Angular frequency, rad/s.

        Returns
        -------
        numpy.ndarray
            The derivative, s/rad.
        """
        index, slope = self._index_and_slope(omega)
        return np.asarray(2 * index * slope)


def _dataset(key: str, where: str) -> Dataset:
    """Return a curated dataset by key.

    Parameters
    ----------
    key : str
        ``<source>/<dataset>``.
    where : str
        The configuration, for messages.

    Returns
    -------
    Dataset
        The dataset.

    Raises
    ------
    DataError
        If there is no such dataset.
    """
    found = _catalogue.by_key().get(key)
    if found is None:
        msg = f"{where}: no dataset {key!r}"
        raise DataError(msg)
    return found


def _shared(
    real: Dataset, imaginary: Dataset, where: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the wavelengths two datasets share, and their values there.

    Parameters
    ----------
    real, imaginary : Dataset
        The segment's datasets.
    where : str
        The configuration, for messages.

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelengths, m, and the two datasets' scaled values at them.

    Raises
    ------
    DataError
        If they share no wavelength.
    """
    first_wavelength = real.vacuum_wavelength()
    second_wavelength = imaginary.vacuum_wavelength()
    first, second = real.scaled_values(), imaginary.scaled_values()
    order = np.argsort(second_wavelength)
    position = np.clip(
        np.searchsorted(second_wavelength[order], first_wavelength),
        0,
        len(order) - 1,
    )
    partner = order[position]
    matched = np.isclose(
        second_wavelength[partner], first_wavelength, rtol=1e-12, atol=0.0
    )
    if not np.any(matched):
        msg = f"{where}: {real.key} and {imaginary.key} share no wavelength"
        raise DataError(msg)
    return (
        first_wavelength[matched],
        first[matched],
        second[partner[matched]],
    )


def _segment_points(
    segment: Segment, material: str, where: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return one segment's wavelengths and complex index inside its range.

    Parameters
    ----------
    segment : Segment
        The segment.
    material : str
        The model's material, which both datasets must be of.
    where : str
        The configuration, for messages.

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelengths, m, increasing, and n + ik there.

    Raises
    ------
    DataError
        If the datasets do not pair, or the range holds no usable point.
    """
    real = _dataset(segment.real, where)
    imaginary = _dataset(segment.imaginary, where)
    if _PAIRS.get(real.quantity) != imaginary.quantity:
        msg = f"{where}: {real.key} ({real.quantity}) cannot pair with "
        msg += f"{imaginary.key} ({imaginary.quantity})"
        raise DataError(msg)
    if {real.material, imaginary.material} != {material}:
        msg = f"{where}: {real.key} and {imaginary.key} are not both {material}"
        raise DataError(msg)
    wavelength, first, second = _shared(real, imaginary, where)
    if real.quantity == "permittivity_real":
        index = complex_index(np.asarray(first + 1j * second))
    elif real.quantity == "refractive_index_decrement":
        index = np.asarray(1 - first + 1j * second)
    else:
        index = np.asarray(first + 1j * second)
    low, high = segment.range_nm
    keep = (wavelength >= low * 1e-9) & (wavelength <= high * 1e-9)
    keep &= np.isfinite(first) & np.isfinite(second)
    if not np.all((np.real(index[keep]) > 0) & (np.imag(index[keep]) > 0)):
        msg = f"{where}: {real.key} has n or k <= 0 in {low}-{high} nm"
        raise DataError(msg)
    if not np.any(keep):
        msg = f"{where}: no point of {real.key} lies in {low}-{high} nm"
        raise DataError(msg)
    order = np.argsort(wavelength[keep])
    return wavelength[keep][order], index[keep][order]


def _henke_points(
    segment: Segment, material: str, where: str
) -> tuple[np.ndarray, np.ndarray]:
    """Return a Henke segment's wavelengths and complex index inside its range.

    Parameters
    ----------
    segment : Segment
        The segment, with ``henke`` and ``atom_density``.
    material : str
        The model's material, which the density must be of.
    where : str
        The configuration, for messages.

    Returns
    -------
    tuple of numpy.ndarray
        Vacuum wavelengths, m, increasing, and n + ik there.

    Raises
    ------
    DataError
        If the density is unusable, or the range holds no usable point.
    """
    density = _dataset(segment.atom_density, where)
    if (
        density.quantity != "atom_density"
        or density.material != material
        or density.parameters.get("counts") != "atoms"
        or len(density) != 1
        or density.unit not in _DENSITY_UNITS
    ):
        msg = f"{where}: {density.key} is not one atom density of {material}"
        raise DataError(msg)
    atoms = float(density.scaled_values()[0]) * _DENSITY_UNITS[density.unit]
    table = np.loadtxt(HENKE_DIRECTORY / f"{segment.henke}.nff", skiprows=1)
    energy_ev, f1, f2 = table[:, 0], table[:, 1], table[:, 2]
    given = f1 > -9000  # CXRO writes -9999 where no f1 is given
    wavelength = 2 * math.pi * HBAR * C / (energy_ev[given] * E_CHARGE)
    factor = CLASSICAL_ELECTRON_RADIUS * wavelength**2 * atoms / (2 * math.pi)
    index = 1 - factor * f1[given] + 1j * factor * f2[given]
    low, high = segment.range_nm
    keep = (wavelength >= low * 1e-9) & (wavelength <= high * 1e-9)
    if not np.any(keep):
        msg = f"{where}: no point of henke/{segment.henke} lies in {low}-{high} nm"
        raise DataError(msg)
    if not np.all((np.real(index[keep]) > 0) & (np.imag(index[keep]) > 0)):
        msg = f"{where}: henke/{segment.henke} has n or k <= 0 in {low}-{high} nm"
        raise DataError(msg)
    order = np.argsort(wavelength[keep])
    return wavelength[keep][order], index[keep][order]


def load_table(path: Path) -> TabulatedModel:
    """Build a table model from its configuration.

    Parameters
    ----------
    path : pathlib.Path
        The configuration, ``data/solid/tables/<material>.toml``.

    Returns
    -------
    TabulatedModel
        The model, named ``table/<material>``.

    Raises
    ------
    DataError
        If the configuration is malformed, the ranges overlap or a dataset is unusable.
    """
    where = path.name
    config = tomllib.loads(path.read_text(encoding="utf-8"))
    material = as_material(str(config["material"]))
    raw = config.get("segments")
    if not isinstance(raw, list) or not raw:
        msg = f"{where}: a table needs [[segments]]"
        raise DataError(msg)
    segments = []
    for entry in raw:
        if not isinstance(entry, Mapping) or not str(entry.get("reason", "")).strip():
            msg = f"{where}: every segment needs real, imaginary, range_nm and a reason"
            raise DataError(msg)
        low, high = (float(v) for v in entry["range_nm"])
        if "henke" in entry and "atom_density" not in entry:
            msg = f"{where}: a Henke segment needs an atom_density dataset"
            raise DataError(msg)
        segments.append(
            Segment(
                real=str(entry.get("real", "")),
                imaginary=str(entry.get("imaginary", "")),
                range_nm=(low, high),
                reason=str(entry["reason"]),
                henke=str(entry.get("henke", "")),
                atom_density=str(entry.get("atom_density", "")),
            )
        )
    segments.sort(key=lambda s: s.range_nm[0])
    for before, after in itertools.pairwise(segments):
        if after.range_nm[0] <= before.range_nm[1]:
            msg = f"{where}: the ranges of {before.label} and {after.label} overlap"
            raise DataError(msg)
    parts = [
        _henke_points(s, material, where)
        if s.henke
        else _segment_points(s, material, where)
        for s in segments
    ]
    knots = np.concatenate([p[0] for p in parts])
    index = np.concatenate([p[1] for p in parts])
    if np.any(np.diff(knots) <= 0):
        msg = f"{where}: repeated wavelengths in the joined table"
        raise DataError(msg)
    cited = dict.fromkeys(source for s in segments for source in s.cited)
    citation = "; ".join(_catalogue.source(source).citation for source in cited)
    return TabulatedModel(
        material=material,
        name=f"table/{material}",
        citation=f"optmat table of {material} (TC-544), from: {citation}",
        wavelength_range=(float(knots[0]), float(knots[-1])),
        kind="tabulated",
        segments=tuple(segments),
        knots=knots,
        index=index,
    )


def _check_table(knots: np.ndarray, index: np.ndarray) -> None:
    """Check a table's knots and index.

    Parameters
    ----------
    knots : numpy.ndarray
        Vacuum wavelengths, m.
    index : numpy.ndarray
        The complex index at them.

    Raises
    ------
    ValueError
        If they differ in length, the knots do not increase, or n or k is not
        positive.
    """
    if (
        knots.shape != index.shape
        or knots.size < 2
        or not np.all(np.diff(knots) > 0)
        or not np.all((np.real(index) > 0) & (np.imag(index) > 0))
    ):
        msg = "a table needs increasing wavelengths and positive n and k, alike in size"
        raise ValueError(msg)


@functools.cache
def tabulated_models() -> dict[str, TabulatedModel]:
    """Return every configured table model, named ``table/<material>``.

    Returns
    -------
    dict of str to TabulatedModel
        The models.
    """
    out = {}
    for path in sorted(TABLES_DIRECTORY.glob("*.toml")):
        table = load_table(path)
        out[table.name] = table
    return curated(out)

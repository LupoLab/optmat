# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract van der Hart and Bingham (2005): R-matrix Floquet two- and three-photon He.

The paper's three tables give the positions and widths of He's autoionising states
(no quantity of the vocabulary holds them; the record's notes); its generalised two-
and three-photon ionisation cross sections are given only in figures, all vector
drawings on the page (no image, no form; the PDF's one image is the download banner of
page 1), read from the drawing (data rule 8, TC-699). The PDF's page 1 is that banner:
journal p. 207 is PDF page 2, so the journal page is the PDF page plus 205.

* **Fig. 2** (PDF p. 6, journal p. 211, library figure ``fig_d93cfccc63b5``): the
  generalised two-photon cross section from the total ionisation rate at 1e13 W/cm2
  (solid), 0.453 to 0.840 au, with Saenz and Lambropoulos's perturbative one (dashed).
* **Fig. 9** (p. 14/219, ``fig_f2b3c82b4132``): the same for three photons, 0.302 to
  0.420 au.
* **Fig. 1** (p. 5/210, ``fig_bf4e94aaacb8``): the total ionisation rates (au) at 1, 2,
  5 and 10 x 10^13 W/cm2, 0.455 to 0.899 au. Not a quantity of the vocabulary (a rate
  against the photon energy at one intensity); read here only for the checks: the 1e13
  curve converted to a cross section must reproduce Fig. 2, and the curves must scale
  as I^2 off resonance, as the text says.
* **Fig. 7** (p. 12/217, ``fig_821c1e0be004``): the same in the three-photon regime,
  0.30 to 0.45 au, read only for the check that its 1e13 curve reproduces Fig. 9.

How the drawing is read:

* every path is stroked with one width per figure (0.315 to 0.325 pt), without Bezier
  segments or forms. The frame: the one closed rectangle of five vertices. The ticks:
  the two-vertex subpaths of the undashed paths, perpendicular to an edge with one end
  on it (within :data:`EDGE`), pointing inwards, short (about 1.05 pt) or long (about
  2.1 pt, longer than :data:`LONG`); the top and right edges repeat the bottom and left
  ones (checked, within :data:`MIRROR`). The tick labels are text objects; their values
  are written in :data:`FIGURES` as the text prints them (checked against the page's
  text). x linear: the long ticks are the labels, the short ones divide each interval
  evenly. y logarithmic: the long ticks are the labelled decades, the short ones 2, ...,
  9 times each decade (Figs. 1 and 2); Figs. 7 and 9 have only their labelled ticks,
  every second decade. Every tick lies on the fitted scale within
  :data:`TICK_TOLERANCE`, the drawing's quantum (about 0.105 pt).
* the curves: the paths with more than :data:`CURVE_VERTICES` vertices, one per curve,
  told apart by their dash patterns (:data:`FIGURES`): solid (none), dot-dashed (dash,
  gap, dot, gap), long dashes and short dashes, as the captions name them. Each path is
  a few polyline subpaths, each starting exactly where the last ended (checked), and the
  joined line never steps back in frequency (checked). Every vertex inside the frame
  is taken, less exact repeats; none lies on the frame's top or bottom edge. The lines
  end on the frame's right edge (0.840 au in Fig. 2, 0.420 au in Fig. 9), where the
  frame cuts them: the value there is the drawn line's.
* the values: rows are numbered in drawing order; each ``printed`` cell gives the
  vertex's position in page points. The abscissa is the frequency in atomic units, the
  photon energy in hartree.

**Checks** (:func:`checks`), against every number the paper prints for these curves
(pp. 209-211, 219-220), between Figs. 1 and 2 and between Figs. 7 and 9; the script
fails if one is missed by more than its stated margin; and the dashed lines against
the curves of saenz_1999, which they reproduce (TC-873). The cross sections that Sato
et al. (2011) quote from this calculation are compared and reported
(:func:`sato_2011`).
"""

import csv
import ctypes
import functools
import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf
from _library import DATA, Point, write
from pypdfium2 import raw

DIRECTORY = "Hart_2005_TwoThreePhoton"
SOURCE = "vanderhart_2005"
EDGE = 0.02
"""A tick's end lies on its frame edge within this, pt."""
LONG = 1.6
"""Ticks longer than this are the long, labelled ones, pt (short 1.05, long 2.1)."""
TICK_MAX = 3.0
"""Ticks are shorter than this, pt."""
MIRROR = 0.02
"""A top or right tick repeats a bottom or left one within this, pt."""
TICK_TOLERANCE = 0.11
"""Every tick lies on its fitted scale within this, pt: the drawing's quantum."""
CURVE_VERTICES = 250
"""Paths with more vertices than this are curves (the tick paths have at most 208)."""

SOLID: tuple[float, ...] = ()
DOT_DASHED = "dot-dashed"
LONG_DASHES = "long dashes"
SHORT_DASHES = "short dashes"


@dataclass(frozen=True)
class Figure:
    """One plot: where it is, what its tick labels print, and its curves' dashes."""

    name: str
    page: int
    x_labels: tuple[float, ...]
    """The frequencies the x labels print, au."""
    x_text: str
    """The x labels' text object, as the page's text layer has it."""
    y_decades: tuple[float, ...]
    """The decades the y labels print (10^-k)."""
    curves: dict[str, int]
    """Each curve's name and its dash pattern's number of entries (0: solid)."""


FIGURES = {
    "fig1": Figure(
        "Fig. 1",
        5,
        (0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
        "0.4 0.5 0.6 0.7 0.8 0.9",
        (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3),
        {"1e13": 0, "2e13": 4, "5e13": 2, "1e14": -2},
    ),
    "fig2": Figure(
        "Fig. 2",
        6,
        (0.45, 0.55, 0.65, 0.75),
        "0.45 0.55 0.65 0.75",
        (1e-52, 1e-51, 1e-50, 1e-49, 1e-48, 1e-47),
        {"rmf": 0, "pt": 2},
    ),
    "fig7": Figure(
        "Fig. 7",
        12,
        (0.30, 0.35, 0.40, 0.45),
        "0.3 0.35 0.4 0.45",
        (1e-10, 1e-8, 1e-6, 1e-4, 1e-2),
        {"1e13": 0, "2e13": 4, "5e13": 2, "1e14": -2},
    ),
    "fig9": Figure(
        "Fig. 9",
        14,
        (0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42),
        "0.3 0.32 0.34 0.36 0.38 0.4 0.42",
        (1e-86, 1e-84, 1e-82, 1e-80, 1e-78),
        {"rmf": 0, "pt": 2},
    ),
}
"""The curves' dash patterns: none (solid), four entries (dash, gap, dot, gap: the
dot-dashed line), two of unequal length (dash 2 gap: long dashes, about 4.3 and 2.15
pt) and two equal (short dashes, 2.15 and 2.15 pt), written here as 0, 4, 2 and -2."""

FILES = {
    ("fig2", "rmf"): "He/vanderhart_2005__fig2_sigma2_rmf_He.csv",
    ("fig2", "pt"): "He/vanderhart_2005__fig2_sigma2_saenz_He.csv",
    ("fig9", "rmf"): "He/vanderhart_2005__fig9_sigma3_rmf_He.csv",
    ("fig9", "pt"): "He/vanderhart_2005__fig9_sigma3_saenz_He.csv",
}


@dataclass(frozen=True)
class Stroke:
    """A stroked path object: its subpaths (page points) and its dash pattern."""

    dash: tuple[float, ...]
    subpaths: tuple[np.ndarray, ...]

    @property
    def kind(self) -> int:
        """The dash pattern as :data:`FIGURES` writes it: 0, 4, 2 or -2."""
        if len(self.dash) == 2 and math.isclose(*self.dash, rel_tol=0.01):
            return -2
        return len(self.dash)


@functools.cache
def _strokes(page: int) -> tuple[Stroke, ...]:
    """Return a page's stroked paths, split into subpaths at each move.

    Raises if a path lies in a form or has a Bezier segment (none does; checked).
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        fill, stroked = ctypes.c_int(), ctypes.c_int()
        raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill), ctypes.byref(stroked))
        if not stroked.value:
            continue
        if obj.container is not None:
            msg = f"p. {page}: a stroked path inside a form, which is not mapped here"
            raise ValueError(msg)
        count = raw.FPDFPageObj_GetDashCount(obj.raw)
        dashes = (ctypes.c_float * max(count, 1))()
        if count > 0:
            raw.FPDFPageObj_GetDashArray(obj.raw, dashes, count)
        matrix = obj.get_matrix()
        subpaths: list[list[tuple[float, float]]] = []
        for index in range(raw.FPDFPath_CountSegments(obj.raw)):
            segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
            kind = raw.FPDFPathSegment_GetType(segment)
            if kind == raw.FPDF_SEGMENT_BEZIERTO:
                msg = f"p. {page}: a Bezier segment, which this script does not read"
                raise ValueError(msg)
            x, y = ctypes.c_float(), ctypes.c_float()
            raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
            point = matrix.on_point(x.value, y.value)
            if kind == raw.FPDF_SEGMENT_MOVETO or not subpaths:
                subpaths.append([])
            subpaths[-1].append((float(point[0]), float(point[1])))
        found.append(
            Stroke(
                tuple(round(float(dashes[k]), 3) for k in range(count)),
                tuple(np.array(s, dtype=float).reshape(-1, 2) for s in subpaths),
            )
        )
    return tuple(found)


@functools.cache
def _frame(key: str) -> tuple[float, float, float, float]:
    """Return a plot's frame (l, b, r, t, pt): its one stroked closed rectangle."""
    found = {
        (
            round(float(sub[:, 0].min()), 3),
            round(float(sub[:, 1].min()), 3),
            round(float(sub[:, 0].max()), 3),
            round(float(sub[:, 1].max()), 3),
        )
        for stroke in _strokes(FIGURES[key].page)
        for sub in stroke.subpaths
        if len(sub) == 5 and np.allclose(sub[0], sub[-1])
    }
    if len(found) != 1:
        msg = f"{FIGURES[key].name}: {len(found)} frames"
        raise ValueError(msg)
    return found.pop()


@dataclass(frozen=True)
class Ticks:
    """A plot's tick positions along one edge, pt: all of them, and the long ones."""

    every: tuple[float, ...]
    long: tuple[float, ...]


def _edge_ticks(key: str) -> dict[str, list[tuple[float, float]]]:
    """Return the ticks on each edge of a plot's frame: (position, length), pt."""
    left, bottom, right, top = _frame(key)
    edges: dict[str, list[tuple[float, float]]] = {
        "bottom": [],
        "top": [],
        "left": [],
        "right": [],
    }
    segments = [
        sub
        for stroke in _strokes(FIGURES[key].page)
        if not stroke.dash
        for sub in stroke.subpaths
        if len(sub) == 2 and math.dist(sub[0], sub[1]) <= TICK_MAX
    ]
    for (x0, y0), (x1, y1) in segments:
        length = math.hypot(x1 - x0, y1 - y0)
        if abs(x0 - x1) < EDGE and left - EDGE <= x0 <= right + EDGE:
            for name, edge in (("bottom", bottom), ("top", top)):
                if abs(y0 - edge) < EDGE:
                    edges[name].append((x0, length))
        if abs(y0 - y1) < EDGE and bottom - EDGE <= y0 <= top + EDGE:
            for name, edge in (("left", left), ("right", right)):
                if abs(x0 - edge) < EDGE:
                    edges[name].append((y0, length))
    return edges


def _ticks(key: str) -> tuple[Ticks, Ticks]:
    """Return the x ticks (bottom edge) and y ticks (left edge) of a plot.

    Checks that the top and right edges repeat them within :data:`MIRROR`.
    """
    edges = _edge_ticks(key)
    out = []
    for near, far in (("bottom", "top"), ("left", "right")):
        every = sorted({round(p, 3) for p, _ in edges[near]})
        long = sorted({round(p, 3) for p, length in edges[near] if length > LONG})
        mirror = sorted({round(p, 3) for p, _ in edges[far]})
        if len(mirror) != len(every) or any(
            abs(m - p) > MIRROR for m, p in zip(mirror, every, strict=True)
        ):
            msg = f"{FIGURES[key].name}: the {far} ticks do not repeat the {near} ones"
            raise ValueError(msg)
        out.append(Ticks(tuple(every), tuple(long)))
    return out[0], out[1]


def _linear_values(ticks: Ticks, labels: Sequence[float], name: str) -> list[float]:
    """Return the values of every tick of a linear axis, from its labelled long ones."""
    if len(ticks.long) != len(labels):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(labels)} labels"
        raise ValueError(msg)
    first = ticks.every.index(ticks.long[0])
    per = (ticks.every.index(ticks.long[-1]) - first) / (len(labels) - 1)
    if per != int(per):
        msg = f"{name}: the short ticks do not divide the labels' intervals evenly"
        raise ValueError(msg)
    step = (labels[1] - labels[0]) / per
    return [labels[0] + (k - first) * step for k in range(len(ticks.every))]


def _log_values(ticks: Ticks, decades: Sequence[float], name: str) -> list[float]:
    """Return the values of every tick of a log axis.

    The long ticks are the decades; the short ones, where drawn, 2, ..., 9 times each.
    """
    if len(ticks.long) != len(decades):
        msg = f"{name}: {len(ticks.long)} long ticks for {len(decades)} decades"
        raise ValueError(msg)
    if ticks.every == ticks.long:
        return list(decades)
    values = [d * k for d in decades[:-1] for k in range(1, 10)] + [decades[-1]]
    if len(ticks.every) != len(values):
        msg = f"{name}: {len(ticks.every)} ticks for {len(values)} values"
        raise ValueError(msg)
    return values


@functools.cache
def _labels_on_page(key: str) -> tuple[str, ...]:
    """Return the text objects of a figure's page (to check the x labels against)."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    page = document[FIGURES[key].page - 1]
    text_page = page.get_textpage()
    return tuple(
        obj.extract().strip()
        for obj in page.get_objects(max_depth=15, textpage=text_page)
        if isinstance(obj, pdfium.PdfTextObj)
    )


@functools.cache
def axes(key: str) -> tuple[Axis, Axis]:
    """Return a plot's frequency (au) and value axes, checked on every tick."""
    figure = FIGURES[key]
    if figure.x_text not in _labels_on_page(key):
        msg = f"{figure.name}: the x labels {figure.x_text!r} are not on the page"
        raise ValueError(msg)
    x_ticks, y_ticks = _ticks(key)
    x_axis = Axis(
        x_ticks.every,
        tuple(_linear_values(x_ticks, figure.x_labels, figure.name)),
        tolerance=TICK_TOLERANCE,
    )
    y_axis = Axis(
        y_ticks.every,
        tuple(_log_values(y_ticks, figure.y_decades, figure.name)),
        log=True,
        tolerance=TICK_TOLERANCE,
    )
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def residuals(key: str) -> tuple[float, float]:
    """Return the worst tick residual of a plot's x and y axes, pt."""
    out = []
    for axis in axes(key):
        y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
        slope, intercept = np.polyfit(axis.positions, y, 1)
        out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
    return out[0], out[1]


@functools.cache
def curve(key: str, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a curve: its vertices (pt), frequencies (au) and values, drawing order."""
    figure = FIGURES[key]
    found = [
        stroke
        for stroke in _strokes(figure.page)
        if stroke.kind == figure.curves[name]
        and sum(len(s) for s in stroke.subpaths) > CURVE_VERTICES
    ]
    if len(found) != 1:
        msg = f"{figure.name}: {len(found)} paths for the curve {name!r}"
        raise ValueError(msg)
    subpaths = found[0].subpaths
    for before, after in itertools.pairwise(subpaths):
        if not np.array_equal(before[-1], after[0]):
            msg = f"{figure.name} {name}: a subpath does not start where the last ended"
            raise ValueError(msg)
    line = np.vstack(subpaths)
    keep = np.ones(len(line), dtype=bool)
    keep[1:] = np.any(np.diff(line, axis=0) != 0, axis=1)
    line = line[keep]
    if np.any(np.diff(line[:, 0]) < 0):
        msg = f"{figure.name} {name}: the line steps back in frequency"
        raise ValueError(msg)
    left, bottom, right, top = _frame(key)
    inside = (
        (line[:, 0] >= left - EDGE)
        & (line[:, 0] <= right + EDGE)
        & (line[:, 1] > bottom + EDGE)
        & (line[:, 1] < top - EDGE)
    )
    if not inside.all():
        msg = f"{figure.name} {name}: {int((~inside).sum())} vertices off the frame"
        raise ValueError(msg)
    x_axis, y_axis = axes(key)
    return (
        line,
        np.asarray(x_axis(line[:, 0]), dtype=float),
        np.asarray(y_axis(line[:, 1]), dtype=float),
    )


def _points(key: str, name: str) -> list[Point]:
    """Return a curve as data points, one per vertex."""
    line, frequency, value = curve(key, name)
    return [
        Point(
            f"v{k:04d}",
            f"{w:.6f}",
            f"{v:.4e}",
            "",
            f"vertex {k} at ({x:.3f}, {y:.3f}) pt",
        )
        for k, ((x, y), w, v) in enumerate(zip(line, frequency, value, strict=True), 1)
    ]


def at(key: str, name: str, frequency: float) -> float:
    """Return a curve's value at a frequency (au), log-interpolated between vertices."""
    _, frequencies, values = curve(key, name)
    if not frequencies[0] <= frequency <= frequencies[-1]:
        msg = f"{FIGURES[key].name} {name}: {frequency} au is outside the curve"
        raise ValueError(msg)
    above = max(int(np.argmax(frequencies >= frequency)), 1)
    w0, w1 = frequencies[above - 1], frequencies[above]
    v0, v1 = values[above - 1], values[above]
    t = 0.0 if w1 == w0 else (frequency - w0) / (w1 - w0)
    return float(10 ** ((1 - t) * math.log10(v0) + t * math.log10(v1)))


def _integral(key: str, name: str, low: float, high: float) -> float:
    """Return a curve's integral over frequency (au) between two frequencies.

    The trapezoidal rule over the drawn vertices, with the ends interpolated.
    """
    _, frequencies, values = curve(key, name)
    inside = (frequencies > low) & (frequencies < high)
    w = np.concatenate([[low], frequencies[inside], [high]])
    v = np.concatenate([[at(key, name, low)], values[inside], [at(key, name, high)]])
    return float(np.sum(np.diff(w) * (v[1:] + v[:-1]) / 2))


# CODATA 2018: the hartree in eV, the atomic unit of time, the electron volt in J.
HARTREE_EV = 27.211386245988
TIME_AU = 2.4188843265857e-17
ELECTRONVOLT = 1.602176634e-19


def rate_to_cross_section(
    rate_au: float, frequency: float, intensity: float, photons: int = 2
) -> float:
    """Return sigma^(N) = Gamma / Phi^N, cm^2N s^(N-1), from an N-photon rate in au.

    Phi = I / (hbar omega) is the photon flux in cm^-2 s^-1, I in W/cm2 and omega the
    frequency in au (the definition the record's notes establish for Figs. 2 and 9).
    """
    flux = intensity / (frequency * HARTREE_EV * ELECTRONVOLT)
    return rate_au / TIME_AU / flux**photons


def _saenz(name: str) -> tuple[np.ndarray, np.ndarray]:
    """Return a committed saenz_1999 curve: photon energies (eV) and values (cm2N s)."""
    path = DATA / "values" / "ionisation" / "He" / f"saenz_1999__{name}_He.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(r["x"]) for r in rows]),
        np.array([float(r["value"]) for r in rows]),
    )


def saenz_total(names: Sequence[str], energy: float) -> float:
    """Return the sum of saenz_1999's partial cross sections at a photon energy (eV).

    Each partial curve is linearly interpolated between its drawn vertices.
    """
    total = 0.0
    for name in names:
        energies, values = _saenz(name)
        total += float(np.interp(energy, energies, values))
    return total


# How far the drawing may stray from the numbers the text prints, and the two figures
# that draw one calculation from each other: the drawing's quantum, 0.105 pt, is 0.7 %
# of a value on Fig. 2's axis (36.3 pt per decade) and 1.1 % on Fig. 9's (22.1 pt per
# decade); between Figs. 1 and 2 the two readings add, and Fig. 1's rate is drawn
# with 36.1 pt per decade.
SAME = 0.03
"""Fig. 1's 1e13 rate as a cross section against Fig. 2, off resonance."""
SCALING = 0.05
"""Fig. 1's curves scaled by I^2 against each other, off resonance ('very close
agreement'): the 1e14 W/cm2 curve falls 2 to 4 % below I^2 times the 1e13 one from
0.47 to 0.70 au; a curve taken for the wrong intensity would be off by 4 to 25 times."""
TEXT = 0.02
"""Beyond the text's 'generally less than 10%' and '15%': two readings, 0.7 to 1.1 %
each."""
STRENGTH = 0.10
"""The resonance strengths' ratio against the text's factor (30 and 5.6), which has
two figures: the integral of a drawn line through its peak depends on the calculation's
energy grid, which the paper does not give (read: 29.6 and 5.75)."""
DASHED = 0.02
"""The dashed lines against saenz_1999's partial curves summed (each reading 0.2 to
0.7 %, and saenz_1999's interpolated between its vertices)."""


@dataclass
class _Checks:
    """The checks' results, and the failures among them."""

    out: dict[str, object]
    failures: list[str]

    def ratio(self, label: str, read: float, against: float, margin: float) -> None:
        """Record read / against - 1, failing beyond a relative margin."""
        self.out[label] = read / against - 1
        if abs(read / against - 1) > margin:
            self.failures.append(f"{label}: {read:.4g} against {against:.4g}")


def _conversions(result: _Checks) -> None:
    """Figs. 1 and 7's 1e13 rates as cross sections against Figs. 2 and 9."""
    # Figs. 2 and 9 are 'obtained through transformation of the total ionization rates
    # obtained at an intensity of 1e13 W/cm2' (captions, pp. 211, 219): Fig. 1's and
    # Fig. 7's 1e13 curves, as Gamma / Phi^N, are their solid lines, off resonance.
    for rates, sigma, photons, frequencies in (
        ("fig1", "fig2", 2, (0.47, 0.50, 0.55, 0.60, 0.65, 0.70)),
        ("fig7", "fig9", 3, (0.31, 0.32, 0.34, 0.36, 0.40, 0.41)),
    ):
        for w in frequencies:
            converted = rate_to_cross_section(at(rates, "1e13", w), w, 1e13, photons)
            result.ratio(
                f"{FIGURES[rates].name} 1e13 as sigma over {FIGURES[sigma].name} at "
                f"{w} au",
                converted,
                at(sigma, "rmf", w),
                SAME,
            )
    # p. 210: 'an I^2-scaling of the present results shows very close agreement
    # between calculations at different intensities, except near resonances'.
    for w in (0.47, 0.50, 0.55, 0.60, 0.65, 0.70):
        base = at("fig1", "1e13", w)
        for name, factor in (("2e13", 2.0), ("5e13", 5.0), ("1e14", 10.0)):
            result.ratio(
                f"Fig. 1 {name} / I^2 over 1e13 at {w} au",
                at("fig1", name, w) / factor**2,
                base,
                SCALING,
            )


def _worst(key: str, low: float, high: float) -> float:
    """Return the largest |RMF / PT - 1| of a figure between two frequencies (au)."""
    return max(
        abs(at(key, "rmf", w) / at(key, "pt", w) - 1)
        for w in np.linspace(low, high, 48)
    )


def _peaks(key: str, window: float) -> tuple[float, float, float, float, float]:
    """Return a figure's perturbative peak and the RMF line's highest point near it.

    (PT peak, au; its value; the RMF line's highest point within a window, au, of the
    PT peak, au; its value; PT over RMF.)
    """
    _, w_pt, v_pt = curve(key, "pt")
    _, w_rmf, v_rmf = curve(key, "rmf")
    peak = int(np.argmax(v_pt))
    near = np.abs(w_rmf - w_pt[peak]) <= window
    top = int(np.argmax(np.where(near, v_rmf, 0.0)))
    return (
        float(w_pt[peak]),
        float(v_pt[peak]),
        float(w_rmf[top]),
        float(v_rmf[top]),
        float(v_pt[peak] / v_rmf[top]),
    )


def _text(result: _Checks) -> None:
    """Check the readings against the text's numbers for Figs. 2 and 9."""
    # p. 211: 'very good agreement ... when resonances do not influence
    # photoionization. Differences are generally less than 10%.' Off resonance, from
    # where both lines are drawn (0.465 au) to the rise of the 1s2p wing (0.70 au);
    # p. 219: 'Differences are generally less than 15%' (both drawn from 0.325 au to
    # the wing of the 1s2s two-photon resonance, 0.38 au).
    for key, low, high, printed in (
        ("fig2", 0.465, 0.70, 0.10),
        ("fig9", 0.325, 0.36, 0.15),
    ):
        worst = _worst(key, low, high)
        result.out[f"{FIGURES[key].name} worst |RMF / PT - 1|, {low}-{high} au"] = worst
        if worst > printed + TEXT:
            result.failures.append(f"{FIGURES[key].name}'s lines differ by {worst:.3f}")
    # p. 211: integrating 'from a photon energy of 0.71 to 0.82 au, we find that the
    # resonance strength in the perturbative regime is a factor 30 larger'; p. 220:
    # 'from a frequency of 0.355 to 0.405 au, we find a difference in resonance
    # strength ... of a factor of 5.6'.
    for key, low, high, printed in (
        ("fig2", 0.71, 0.82, 30.0),
        ("fig9", 0.355, 0.405, 5.6),
    ):
        result.ratio(
            f"{FIGURES[key].name} resonance strength PT / RMF, {low}-{high} au "
            f"(text: {printed})",
            _integral(key, "pt", low, high) / _integral(key, "rmf", low, high),
            printed,
            STRENGTH,
        )
    # p. 211: near the 1s2p resonance 'the difference amounts to three orders of
    # magnitude'; p. 219: 'Near the 1s2s resonance, there is a factor of 10
    # difference': the drawn peaks (reported). At 1s2p the RMF line at the PT peak's
    # frequency (within the drawing's quantum, 0.0002 au; nearby are the RMF line's
    # three-photon spikes); at 1s2s its own peak, 0.0003 au higher (window 0.001 au).
    result.out["Fig. 2 PT peak, RMF at it (au, value each), PT/RMF"] = _peaks(
        "fig2", 0.0002
    )
    result.out["Fig. 9 PT peak, RMF peak (au, value each), PT/RMF"] = _peaks(
        "fig9", 0.001
    )


def checks() -> dict[str, object]:
    """Return the checks of the reading against the text and between figures.

    Raises if one fails its margin (:data:`SAME`, :data:`SCALING`,
    :data:`STRENGTH`, :data:`DASHED`, and the text's own 10 % and 15 % with
    :data:`TEXT`).
    """
    result = _Checks({}, [])
    _conversions(result)
    _text(result)
    # The dashed lines are Saenz and Lambropoulos's perturbative cross sections ('Dr A
    # Saenz ... providing his results in numerical form', p. 220): the sum of their
    # partial curves (saenz_1999), on their axes (photon energy, eV, over the hartree).
    for label, (dashed, summed) in saenz_1999().items():
        result.ratio(f"dashed over saenz_1999 summed, {label}", dashed, summed, DASHED)
    if result.failures:
        raise ValueError("; ".join(result.failures))
    return result.out


SATO = {61.4: 8.5e-52, 58.4: 63e-52, 56.0: 3.5e-52, 53.4: 6.9e-52}
"""Sato et al. (2011), Table 2, 'Theory (R-matrix Floquet method)' at 1e13 W/cm2,
cm4 s, at the wavelengths (nm) of their measurement: this calculation, supplied by
van der Hart (sato_2011, TC-876)."""
NM_EV = 1239.84198
"""hc in eV nm (CODATA 2018)."""


def sato_2011() -> dict[float, tuple[float, float, float]]:
    """Return, per wavelength, (frequency au, this paper's value, Sato's quoted one).

    From Fig. 2 where it is drawn (to 0.840 au); at 53.4 nm (0.853 au) from Fig. 1's
    1e13 rate converted as Fig. 2 is.
    """
    out = {}
    for nm, quoted in SATO.items():
        w = NM_EV / nm / HARTREE_EV
        if w <= curve("fig2", "rmf")[1][-1]:
            value = at("fig2", "rmf", w)
        else:
            value = rate_to_cross_section(at("fig1", "1e13", w), w, 1e13)
        out[nm] = (w, value, quoted)
    return out


def saenz_1999() -> dict[str, tuple[float, float]]:
    """Return the dashed lines and saenz_1999's partial curves summed, per frequency.

    At frequencies off the resonances, on saenz_1999's axes as drawn (photon energy,
    eV, the hartree times the frequency): (dashed line, sum of the partials).
    """
    out = {}
    for w in (0.47, 0.50, 0.55, 0.60, 0.65, 0.70, 0.74, 0.82):
        energy = w * HARTREE_EV
        out[f"two-photon at {w} au"] = (
            at("fig2", "pt", w),
            saenz_total(("fig1a_sigma2_1S", "fig1b_sigma2_1D"), energy),
        )
    for w in (0.33, 0.34, 0.35, 0.36, 0.37, 0.39, 0.40, 0.41):
        energy = w * HARTREE_EV
        out[f"three-photon at {w} au"] = (
            at("fig9", "pt", w),
            saenz_total(("fig2a_sigma3_1P", "fig2b_sigma3_1F"), energy),
        )
    return out


def extract() -> dict[str, list[Point]]:
    """Return every taken curve, after the checks."""
    checks()
    return {relative: _points(*key) for key, relative in FILES.items()}


def main() -> None:
    """Write the data files, and print the checks and the comparisons."""
    write(extract())
    for key, name in [*FILES, ("fig1", "1e13"), ("fig1", "1e14")]:
        line, frequencies, values = curve(key, name)
        x_res, y_res = residuals(key)
        print(
            f"{FIGURES[key].name} {name}: {len(line)} vertices, {frequencies[0]:.4f} "
            f"to {frequencies[-1]:.4f} au, {values.min():.3e} to {values.max():.3e}; "
            f"ticks within {x_res:.3f} (x) and {y_res:.3f} (y) pt"
        )
    for label, result in checks().items():
        print(f"check {label}: {result}")
    for nm, (w, value, quoted) in sato_2011().items():
        print(
            f"sato_2011 {nm} nm ({w:.5f} au): read {value:.3e}, quoted {quoted:.3e}, "
            f"ratio {quoted / value:.3f}"
        )


if __name__ == "__main__":
    main()

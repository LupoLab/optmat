# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Read values from a paper's figures (data rule 8, TC-699).

A figure drawn as vectors in the PDF is read from the drawing: :func:`paths` returns
the path objects of a region of a page with their vertices in page points, and
:class:`Axis` maps a coordinate to a value, calibrated from the positions of the
axis's own tick marks and checked against every one of them. :func:`painted` does the
same as :func:`paths` inside form XObjects too, with each path's colours, and
:func:`texts` returns the text objects of a region. A figure embedded as an
image is read from its pixels at their native resolution (:func:`image`); the same
:class:`Axis` maps a pixel to a value, calibrated from pixel positions written in the
extraction script.

The PDF is the library's original (``<directory>/original/source_01/*.pdf``).
"""

import ctypes
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pypdfium2 as pdfium
from _library import LIBRARY
from pypdfium2 import raw


def pdf(directory: str) -> Path:
    """Return the one PDF of a library entry's original source."""
    found = sorted((LIBRARY / directory / "original" / "source_01").glob("*.pdf"))
    if len(found) != 1:
        msg = f"{directory}: expected one original PDF, found {len(found)}"
        raise ValueError(msg)
    return found[0]


@dataclass(frozen=True)
class Drawn:
    """A vector path: its bounding box and its vertices, in page points."""

    left: float
    bottom: float
    right: float
    top: float
    vertices: np.ndarray

    @property
    def centre(self) -> tuple[float, float]:
        """The centre of the bounding box (a marker's position)."""
        return (self.left + self.right) / 2, (self.bottom + self.top) / 2

    @property
    def size(self) -> float:
        """The larger side of the bounding box."""
        return max(self.right - self.left, self.top - self.bottom)


def _vertices(obj: pdfium.PdfObject) -> np.ndarray:
    """Return a path object's vertices in page coordinates."""
    matrix = obj.get_matrix()
    points = []
    for index in range(raw.FPDFPath_CountSegments(obj.raw)):
        segment = raw.FPDFPath_GetPathSegment(obj.raw, index)
        x, y = ctypes.c_float(), ctypes.c_float()
        raw.FPDFPathSegment_GetPoint(segment, ctypes.byref(x), ctypes.byref(y))
        points.append(matrix.on_point(x.value, y.value))
    return np.array(points, dtype=float).reshape(-1, 2)


def paths(
    directory: str, page: int, region: tuple[float, float, float, float]
) -> list[Drawn]:
    """Return the vector paths of a page (1-based) that touch a region (l, b, r, t, pt).

    A path may extend beyond the region: a curve clipped by the plot's frame keeps its
    vertices outside it, which the caller drops.
    """
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(directory)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        west, south, east, north = obj.get_bounds()
        if east >= left and west <= right and north >= bottom and south <= top:
            found.append(Drawn(west, south, east, north, _vertices(obj)))
    return found


def image(directory: str, page: int, index: int = 0) -> np.ndarray:
    """Return a page's (1-based) embedded image as grey levels 0-255, native size."""
    document = pdfium.PdfDocument(str(pdf(directory)))
    images = [
        obj
        for obj in document[page - 1].get_objects(max_depth=15)
        if isinstance(obj, pdfium.PdfImage)
    ]
    bitmap = images[index].get_bitmap(render=False).to_numpy()
    if bitmap.ndim == 3:  # BGR(A): the mean of the colour channels
        bitmap = bitmap[..., :3].mean(axis=2)
    return np.asarray(bitmap, dtype=float)


@dataclass(frozen=True)
class Axis:
    """A linear or logarithmic axis, calibrated from its tick marks.

    Attributes
    ----------
    positions : tuple of float
        The ticks' coordinates (page points or pixels).
    values : tuple of float
        The values their labels print.
    log : bool
        True for a logarithmic axis.
    tolerance : float
        The largest residual allowed of any tick, in the coordinate's unit.
    """

    positions: tuple[float, ...]
    values: tuple[float, ...]
    log: bool = False
    tolerance: float = 0.05

    def _fit(self) -> tuple[float, float]:
        """Return the slope and intercept of the value (or its log10) at a position."""
        y = np.log10(self.values) if self.log else np.asarray(self.values)
        slope, intercept = np.polyfit(self.positions, y, 1)
        return float(slope), float(intercept)

    def check(self) -> None:
        """Raise unless every tick lies on the calibration within the tolerance."""
        slope, intercept = self._fit()
        y = np.log10(self.values) if self.log else np.asarray(self.values)
        predicted = (y - intercept) / slope
        worst = float(np.max(np.abs(predicted - np.asarray(self.positions))))
        if worst > self.tolerance:
            msg = f"the ticks do not lie on one scale (worst residual {worst:.3g})"
            raise ValueError(msg)

    def __call__(self, position: float | np.ndarray) -> float | np.ndarray:
        """Return the value at a coordinate."""
        slope, intercept = self._fit()
        y = slope * np.asarray(position, dtype=float) + intercept
        return 10**y if self.log else y

    def resolution(self) -> float:
        """Return the value change per unit coordinate: relative for a log axis."""
        slope, _ = self._fit()
        return abs(slope) * math.log(10) if self.log else abs(slope)


def unique(
    points: Sequence[tuple[float, float]], digits: int = 3
) -> list[tuple[float, float]]:
    """Return points without repeats (a marker drawn as fill and outline), sorted."""
    return sorted({(round(x, digits), round(y, digits)) for x, y in points})


@dataclass(frozen=True)
class Painted(Drawn):
    """A vector path with its paint, in page points.

    A figure placed as a form XObject draws its paths in the form's coordinates;
    :func:`painted` maps them to the page through every enclosing form's matrix.
    """

    stroke: tuple[int, int, int] | None = None
    """The stroke colour (RGB, 0-255), or None if the path is not stroked."""
    fill: tuple[int, int, int] | None = None
    """The fill colour (RGB, 0-255), or None if the path is not filled."""
    width: float = 0.0
    """The stroke width, in the path's own coordinates."""


def _to_page(obj: pdfium.PdfObject, points: np.ndarray) -> np.ndarray:
    """Map points from an object's container (a form, or the page) to the page."""
    container = obj.container
    out = np.array(points, dtype=float).reshape(-1, 2)
    while container is not None:
        matrix = container.get_matrix()
        out = np.array([matrix.on_point(x, y) for x, y in out], dtype=float)
        out = out.reshape(-1, 2)
        container = container.container
    return out


def _colour(obj: pdfium.PdfObject, stroke: bool) -> tuple[int, int, int] | None:
    """Return a path's stroke or fill colour, or None if it is not so painted."""
    fill_mode, stroked = ctypes.c_int(), ctypes.c_int()
    raw.FPDFPath_GetDrawMode(obj.raw, ctypes.byref(fill_mode), ctypes.byref(stroked))
    if not (stroked.value if stroke else fill_mode.value):
        return None
    channels = [ctypes.c_uint() for _ in range(4)]
    get = raw.FPDFPageObj_GetStrokeColor if stroke else raw.FPDFPageObj_GetFillColor
    if not get(obj.raw, *(ctypes.byref(c) for c in channels)):
        return None
    return channels[0].value, channels[1].value, channels[2].value


def painted(
    directory: str, page: int, region: tuple[float, float, float, float]
) -> list[Painted]:
    """Return the paths of a page (1-based) touching a region (l, b, r, t, pt), painted.

    As :func:`paths`, but also for paths inside form XObjects (mapped to the page
    through the forms' matrices), and with each path's stroke and fill colours and
    stroke width.
    """
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(directory)))
    found = []
    for obj in document[page - 1].get_objects(max_depth=15):
        if raw.FPDFPageObj_GetType(obj.raw) != raw.FPDF_PAGEOBJ_PATH:
            continue
        west, south, east, north = obj.get_bounds()
        corners = _to_page(obj, np.array([[west, south], [east, north]]))
        west, east = sorted(corners[:, 0])
        south, north = sorted(corners[:, 1])
        if not (east >= left and west <= right and north >= bottom and south <= top):
            continue
        width = ctypes.c_float()
        raw.FPDFPageObj_GetStrokeWidth(obj.raw, ctypes.byref(width))
        found.append(
            Painted(
                float(west),
                float(south),
                float(east),
                float(north),
                _to_page(obj, _vertices(obj)),
                _colour(obj, stroke=True),
                _colour(obj, stroke=False),
                float(width.value),
            )
        )
    return found


def texts(
    directory: str, page: int, region: tuple[float, float, float, float]
) -> list[tuple[str, float, float]]:
    """Return the text objects of a page (1-based) inside a region (l, b, r, t, pt).

    Each as (text, x, y): its text and the centre of its bounding box in page points,
    also inside form XObjects. A figure's glyphs drawn as paths are not text.
    """
    left, bottom, right, top = region
    document = pdfium.PdfDocument(str(pdf(directory)))
    pdf_page = document[page - 1]
    text_page = pdf_page.get_textpage()
    found = []
    for obj in pdf_page.get_objects(max_depth=15, textpage=text_page):
        if not isinstance(obj, pdfium.PdfTextObj):
            continue
        west, south, east, north = obj.get_bounds()
        x, y = _to_page(obj, np.array([[(west + east) / 2, (south + north) / 2]]))[0]
        if left <= x <= right and bottom <= y <= top:
            found.append((obj.extract().strip(), float(x), float(y)))
    return found

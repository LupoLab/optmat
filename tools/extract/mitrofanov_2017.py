# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Check Mitrofanov et al. (2017): the statements the record quotes, and Figure 1d.

The paper prints no table. Its one measured quantity, the differential group delay
Delta tau_g of 60 m of room air at 3.70-4.15 um, is shown only in Figure 1d, an
embedded raster image (page 3, 1650 x 1478 pixels), with a single number in the text:
'Delta tau_g ~ 1.8 ps ... within the wavelength range from 3.70 to 4.15 um' (p. 6).
The record keeps no dataset: optmat's vocabulary has no quantity for a group delay with
an unstated zero, and air is not a registered material (see the record's notes).

This script writes no data file. ``extract()`` asserts that every statement the record
quotes is in the PDF's text layer on the page it cites, and that Figure 1d's axes are
where :func:`calibration` says (the tick marks found in the image, at its native
resolution). :func:`trace` reads the transform-limited probe's curve (navy) by the rule
written there; ``main()`` prints it. The trace is the comparison quoted in the record's
notes, not data: it enters no dataset.
"""

import re
from functools import cache

import numpy as np
import pypdfium2 as pdfium
from _figure import Axis, pdf

DIRECTORY = "Mitrofanov_2017_MappingAnomalousDispersion"

# (pdf page, phrase): each must be in that page's text, whitespace collapsed and the
# minus sign (U+2212) written as a hyphen.
STATEMENTS = (
    (2, "B1=0.05792105 µm-2 , B2=0.00167917µm-2 , C1=238.0185µm-2 , and C2=57.362"),
    (2, "NHO2 ≈5.4∙1016 cm-3"),
    (3, "n0=1.000270232 is the refractive index of air at λ≈3.9 µm"),
    (3, "(d) The differential group delay Δτg induced by 60m of atmospheric air"),
    (4, "thus propagating a path of L≈60 m in atmospheric air"),
    (4, "spectrum stretching from 3.6 to 4.2μm"),
    (4, "the differential group delay Δτg=τg-τ0"),
    (6, "using a standard 7798 CO2 datalogger"),
    (6, "≈1.32 ∙ 1016 cm-3 taken from the readings of the CO2 monitor"),
    (6, "Δτg≈1.8ps are achieved within the wavelength range from 3.70 to 4.15μm"),
)

# Figure 1d's calibration, in pixels of the embedded image (row 0 at the top). The
# tick marks are the black runs (all channels below 70) just outside the axes: x ticks
# on row 1340 below the x axis (rows 1334-1338), y ticks on column 158 left of the
# y axis (columns 160-163). The x axis line ends at column 780, 4.150 um.
X_TICKS = ((161.5, 3.7), (299.0, 3.8), (436.5, 3.9), (574.0, 4.0), (711.5, 4.1))
Y_TICKS = ((1084.5, 1.5), (1152.5, 1.0), (1220.5, 0.5), (1288.5, 0.0))
PANEL_ROWS = (1040, 1333)  # inside panel d, above the x axis
PANEL_COLUMNS = (165, 785)  # inside panel d, short of panel e


def _text(page: int) -> str:
    """Return a page's (1-based) text layer, normalised as STATEMENTS are written.

    Whitespace is collapsed, the line-break hyphens (U+FFFE) are joined and the minus
    sign (U+2212) becomes a hyphen.
    """
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    raw = document[page - 1].get_textpage().get_text_range()
    raw = raw.replace(chr(0xFFFE), "").replace(chr(0x2212), "-")
    return re.sub(r"\s+", " ", raw)


@cache
def _image() -> np.ndarray:
    """Return Figure 1, page 3's one embedded image, as RGB integers, native size."""
    document = pdfium.PdfDocument(str(pdf(DIRECTORY)))
    images = [
        obj
        for obj in document[2].get_objects(max_depth=15)
        if isinstance(obj, pdfium.PdfImage)
    ]
    assert len(images) == 1, len(images)
    bitmap = images[0].get_bitmap(render=False)
    assert (bitmap.width, bitmap.height, bitmap.mode) == (1650, 1478, "BGR")
    return np.asarray(bitmap.to_numpy()[..., ::-1], dtype=int)


def _runs(mask: np.ndarray) -> list[float]:
    """Return the centres of the runs of True in a 1-D mask."""
    index = np.flatnonzero(mask)
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(g[0] + g[-1]) / 2 for g in groups if g.size]


def calibration() -> tuple[Axis, Axis]:
    """Return panel d's wavelength (um) and delay (ps) axes, checked on the image."""
    black = _image().max(axis=2) < 70
    x_found = [150 + c for c in _runs(black[1340, 150:830])]
    y_found = [1040 + r for r in _runs(black[1040:1340, 158])]
    assert x_found == [p for p, _ in X_TICKS], x_found
    assert y_found == [p for p, _ in Y_TICKS], y_found
    assert int(np.flatnonzero(black[1336, :900]).max()) == 780
    x_axis = Axis(tuple(p for p, _ in X_TICKS), tuple(v for _, v in X_TICKS))
    y_axis = Axis(tuple(p for p, _ in Y_TICKS), tuple(v for _, v in Y_TICKS))
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def _colour(name: str) -> np.ndarray:
    """Return the mask of a curve's pixels.

    navy, the transform-limited probe: blue 91-199, red and green below 80 (the
    bright-blue annotation, blue above 200, fails it); green, the chirped probe: green
    above 90, red and blue below 90; magenta, the authors' calculation (the caption's
    'red line'): red and blue above 180, green below 120.
    """
    rgb = _image()
    red, green, blue = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    masks = {
        "navy": (blue > 90) & (blue < 200) & (red < 80) & (green < 80),
        "green": (green > 90) & (red < 90) & (blue < 90),
        "magenta": (red > 180) & (blue > 180) & (green < 120),
    }
    return masks[name]


def trace(name: str = "navy") -> list[tuple[float, float]]:
    """Return one curve of panel d, um and ps, one point per column.

    The rule: per column of the panel the curve is the mean row of its pixels of the
    colour (:func:`_colour`), where there are at least two. Where another curve is
    drawn over it the column has fewer and is skipped.
    """
    x_axis, y_axis = calibration()
    mask = _colour(name)
    top, bottom = PANEL_ROWS
    out = []
    for column in range(*PANEL_COLUMNS):
        rows = np.flatnonzero(mask[top:bottom, column]) + top
        if rows.size >= 2:
            out.append((float(x_axis(column)), float(y_axis(rows.mean()))))
    return out


def at(points: list[tuple[float, float]], wavelength_um: float) -> float:
    """Return the trace's mean within 0.004 um (about 5 columns) of a wavelength.

    NaN where the curve is hidden there.
    """
    values = [v for x, v in points if abs(x - wavelength_um) < 0.004]
    return float(np.mean(values)) if values else float("nan")


def extract() -> dict:
    """Assert the quoted statements and the figure's calibration; return no files."""
    for page, phrase in STATEMENTS:
        assert phrase in _text(page), (page, phrase)
    points = trace()
    # The record quotes the trace relative to 3.70 um; its end, at the frame's corner
    # (4.150 um, 1.5 ps), gives the text's 1.8 ps.
    assert abs((at(points, 4.15) - at(points, 3.70)) - 1.78) < 0.005
    return {}


def main() -> None:
    """Run the checks and print the trace relative to 3.70 um."""
    extract()
    curves = {name: trace(name) for name in ("navy", "green", "magenta")}
    print("um     navy (TL)  green (chirped)  magenta (calc.), ps as read")
    for wavelength in (
        3.70,
        3.75,
        3.80,
        3.85,
        3.90,
        3.95,
        4.00,
        4.05,
        4.08,
        4.10,
        4.12,
        4.13,
        4.14,
        4.15,
    ):
        cells = "  ".join(f"{at(curves[n], wavelength):+8.3f}" for n in curves)
        print(f"{wavelength:.2f}  {cells}")


if __name__ == "__main__":
    main()

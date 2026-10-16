# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Starace and Jiang (1987): transition-matrix two-photon cross sections of Ar.

The paper prints its generalised two-photon ionisation cross sections of the Ar 3p
subshell below the one-photon threshold only in figures (length form, lowest-order
perturbation theory); no table, no supplement. The PDF is a scan: every page is one
embedded 1-bit image, 5169 x 7129 pixels (about 600 dpi). The figures are read from
those native pixels (data rule 8, TC-699). Journal page = PDF page + 1704.

* **Fig. 2** (p. 9/1713, ``fig_3bc749c424e4``): the 3p5 epsilon-f 1D channel, linear
  polarisation, y in units of 10^-51 cm^4 s, 0 to 16: four approximations, HF
  (dash-dot), HF + ground-state correlations GSC (dash-double-dot), + intermediate-state
  interchannel coupling IIC (solid), + electron-scattering interactions ESI (dashed).
* **Fig. 3** (p. 10/1714, ``fig_87f36d953f72``): 3p5 epsilon-p 1D, 0 to 5e-51; HF, IIC
  and ESI (GSC 'not shown separately as they are nearly the same as the IIC results').
* **Fig. 4** (p. 10/1714, ``fig_4f3d33d66e2a``): 3p5 epsilon-p 1S, 0 to 16e-51; all four
  approximations.
* **Fig. 5** (p. 11/1715, ``fig_ca863326b494``): epsilon-f 1D, IIC (solid) and ESI
  (dashed) shifted to the experimental resonance energies; Pindzola and Kelly's 'LGI'
  (dot-dash) is not read.
* **Fig. 6** (p. 12/1716, ``fig_ea1ce7340e25``): the total for circular polarisation,
  0 to 3.5e-50; IIC (solid) and ESI (dashed). McGuire's (dotted), Pindzola and Kelly's
  (dot-dash) and Moccia, Rahman and Rizzo's (dash-double-dot) curves are not read.
* **Fig. 7** (p. 13/1717, ``fig_167e2294aa4c``): the same for linear polarisation.

How the images are read:

* the axes: the tick marks, found ``DEPTH`` pixels inside each frame line and checked
  against the positions written in :data:`FIGURES` (each within ``TICK_FOUND``); each
  side's ticks lie on one linear scale within ``TICK_TOLERANCE`` pixels
  (:meth:`_figure.Axis.check`). Ticks a curve runs through are left out of the lists.
  The scan is rotated by about 0.2 degrees and slightly stretched: the photon energy is
  interpolated between the bottom and top ticks' scales by the row, the value between
  the left and right ticks' scales by the column (:class:`Calibration`).
* the resonance lines: the vertical lines A, B, C (experimental) and a, b, c (HF) of
  Table III, found ``LINE_DEPTH`` pixels above the bottom frame line; their energies on
  the calibration are checked against Table III (:func:`checks`), an independent check
  of the x axis.
* the curves: each stretch of a curve (:data:`SEGMENTS`) is traced column by column from
  a seed pixel on it (:func:`trace`), and its line style is checked from the trace
  (:func:`style`): a solid line has no gap; a dashed one has gaps with nothing in them;
  a dash-dot one one dot per gap, a dash-double-dot one two. This is what tells the
  authors' curves (and which approximation) from the others' and from the resonance
  lines.
* the values: every 0.1 eV within a stretch, the centre of the curve's run at that
  photon energy's column where the run is clear: as thin as one line at its slope, no
  other dark pixel within ``CLEAR`` pixels (about a line width) of it, and so for the
  line's runs ``WINDOW`` columns either side, away from the tick marks (:func:`_read`);
  across a gap of the dashes or a resonance line, the linear interpolation between the
  nearest clear columns. Where another curve comes closer, or the line is steeper than
  ``MAX_SLOPE`` rows per column, the sample is left out: :data:`LEFT_OUT` lists them,
  and the script fails if they change. The ``printed`` cell of each value gives its
  pixel (column, row).

One pixel is 0.0057 to 0.0050 eV in x; in y 0.0020e-50 cm^4 s (Figs. 6, 7), 0.0085e-51
(Fig. 2), 0.0029e-51 (Fig. 3), 0.0093e-51 (Fig. 4), 0.010e-51 (Fig. 5). The lines are
4 to 8 pixels wide.

**Checks** (:func:`checks`, run by :func:`main`; it fails on a miss): the drawn
resonance lines against Table III; Table II's f_q (the circular total is 3/2 the linear
1D channels); Eq. 24's unit against CODATA; Fig. 6 against 3/2 times Figs. 2 + 3, Fig. 7
against Figs. 2 + 3 + 4 (the text's statements of how they are made, pp. 1716-1717).
Reported: the ratios between approximations the text describes, Fig. 7 against McKenna
and van der Hart 2004 (TC-856), Fig. 5 (shifted) against Fig. 2.
"""

import itertools
import json
import math
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import LIBRARY, Point, table, write

from optmat.datasets import datasets

DIRECTORY = "Starace_1987_TransitionMatrixTheory"
DARK = 128


@dataclass(frozen=True)
class Figure:
    """One plot on a scanned page: its frame and tick marks in image pixels.

    The frame lines' rows and columns are those of their centres near the frame's
    middle. Ticks are (pixel, printed value): the bottom and top ones (columns, eV) as
    found ``DEPTH`` pixels inside the bottom and top frame lines, the left and right
    ones (rows, the y axis's printed number) ``DEPTH`` pixels inside the left and right
    ones; the resonance lines (label, column) ``LINE_DEPTH`` pixels above the bottom.
    """

    page: int
    table: str
    unit: float
    top: int
    bottom: int
    left: int
    right: int
    x_bottom: tuple[tuple[float, float], ...]
    x_top: tuple[tuple[float, float], ...]
    y_left: tuple[tuple[float, float], ...]
    y_right: tuple[tuple[float, float], ...]
    lines: tuple[tuple[str, float], ...]


DEPTH = 20
LINE_DEPTH = 60
TICK_TOLERANCE = 3.5
TICK_FOUND = 1.0

FIGURES = {
    # Fig. 2, epsilon-f 1D, 10^-51 cm^4 s: the 4 tick on the left is under the curves.
    "fig2": Figure(
        9,
        "fig_3bc749c424e4",
        1e-51,
        3201,
        5087,
        3191,
        4514,
        (
            (3369.5, 9),
            (3548.5, 10),
            (3723.5, 11),
            (3898.0, 12),
            (4073.0, 13),
            (4252.0, 14),
            (4429.0, 15),
        ),
        (
            (3362.5, 9),
            (3542.0, 10),
            (3716.0, 11),
            (3890.5, 12),
            (4064.5, 13),
            (4243.5, 14),
            (4419.0, 15),
        ),
        (
            (4969.0, 1),
            (4848.0, 2),
            (4731.0, 3),
            (4495.5, 5),
            (4378.5, 6),
            (4260.0, 7),
            (4144.0, 8),
            (4026.0, 9),
            (3908.5, 10),
            (3790.5, 11),
            (3672.5, 12),
            (3557.0, 13),
            (3438.0, 14),
            (3318.5, 15),
        ),
        (
            (4969.0, 1),
            (4849.0, 2),
            (4732.5, 3),
            (4615.0, 4),
            (4496.0, 5),
            (4378.5, 6),
            (4260.0, 7),
            (4142.5, 8),
            (4024.5, 9),
            (3906.0, 10),
            (3788.5, 11),
            (3669.5, 12),
            (3554.5, 13),
            (3435.0, 14),
            (3314.5, 15),
        ),
        (
            ("A", 3844.5),
            ("a", 3936.5),
            ("B", 4270.5),
            ("C", 4285.0),
            ("b", 4332.5),
            ("c", 4354.5),
        ),
    ),
    # Fig. 3, epsilon-p 1D, 10^-51 cm^4 s: the bottom 14 eV tick and the left 2.0 one
    # are under curves.
    "fig3": Figure(
        10,
        "fig_87f36d953f72",
        1e-51,
        902,
        2636,
        963,
        2358,
        (
            (1148.0, 9),
            (1336.0, 10),
            (1519.5, 11),
            (1709.0, 12),
            (1896.0, 13),
            (2271.0, 15),
        ),
        (
            (1142.0, 9),
            (1330.0, 10),
            (1513.5, 11),
            (1704.0, 12),
            (1890.0, 13),
            (2076.0, 14),
            (2265.0, 15),
        ),
        (
            (2463.0, 0.5),
            (2290.5, 1.0),
            (2116.0, 1.5),
            (1769.0, 2.5),
            (1597.5, 3.0),
            (1422.0, 3.5),
            (1249.0, 4.0),
            (1076.0, 4.5),
        ),
        (
            (2460.5, 0.5),
            (2289.0, 1.0),
            (2114.5, 1.5),
            (1939.5, 2.0),
            (1764.0, 2.5),
            (1593.0, 3.0),
            (1418.0, 3.5),
            (1245.5, 4.0),
            (1071.5, 4.5),
        ),
        (
            ("A", 1644.0),
            ("a", 1760.0),
            ("B", 2104.5),
            ("C", 2115.5),
            ("b", 2172.0),
            ("c", 2193.5),
        ),
    ),
    # Fig. 4, epsilon-p 1S, 10^-51 cm^4 s: the left 6 and 7 ticks are under curves.
    "fig4": Figure(
        10,
        "fig_4f3d33d66e2a",
        1e-51,
        4155,
        5877,
        1007,
        2359,
        (
            (1182.5, 9),
            (1362.5, 10),
            (1547.0, 11),
            (1727.0, 12),
            (1911.0, 13),
            (2093.0, 14),
            (2273.5, 15),
        ),
        (
            (1173.0, 9),
            (1353.0, 10),
            (1537.5, 11),
            (1717.0, 12),
            (1902.0, 13),
            (2083.5, 14),
            (2263.5, 15),
        ),
        (
            (5770.5, 1),
            (5664.0, 2),
            (5555.5, 3),
            (5452.5, 4),
            (5343.5, 5),
            (5020.0, 8),
            (4913.5, 9),
            (4806.0, 10),
            (4697.5, 11),
            (4590.5, 12),
            (4482.0, 13),
            (4376.0, 14),
            (4266.0, 15),
        ),
        (
            (5766.5, 1),
            (5660.5, 2),
            (5551.5, 3),
            (5448.0, 4),
            (5339.0, 5),
            (5231.5, 6),
            (5123.5, 7),
            (5015.0, 8),
            (4907.0, 9),
            (4800.0, 10),
            (4691.0, 11),
            (4584.5, 12),
            (4475.0, 13),
            (4370.0, 14),
            (4259.5, 15),
        ),
        (
            ("A", 1665.0),
            ("a", 1776.0),
            ("B", 2104.5),
            ("C", 2119.5),
            ("b", 2175.0),
            ("c", 2191.5),
        ),
    ),
    # Fig. 5, epsilon-f 1D shifted, 10^-51 cm^4 s: no 15 eV tick (the right frame line
    # is about 15.0 eV); the left 4 tick is under the curves. Only A, B and C are drawn;
    # A is the left of two solid vertical lines; the other, 0.04 eV above and not
    # labelled, is probably the IIC profile's steep part.
    "fig5": Figure(
        11,
        "fig_ca863326b494",
        1e-51,
        815,
        2411,
        3151,
        4521,
        (
            (3342.0, 9),
            (3542.0, 10),
            (3734.0, 11),
            (3931.5, 12),
            (4128.5, 13),
            (4322.0, 14),
        ),
        (
            (3341.0, 9),
            (3541.0, 10),
            (3734.5, 11),
            (3931.5, 12),
            (4126.5, 13),
            (4321.0, 14),
        ),
        (
            (2311.5, 1),
            (2211.0, 2),
            (2111.0, 3),
            (1913.0, 5),
            (1813.0, 6),
            (1713.0, 7),
            (1614.0, 8),
            (1514.0, 9),
            (1416.0, 10),
            (1316.0, 11),
            (1217.0, 12),
            (1119.0, 13),
            (1020.0, 14),
            (923.5, 15),
        ),
        (
            (2312.5, 1),
            (2211.0, 2),
            (2113.0, 3),
            (2010.5, 4),
            (1912.0, 5),
            (1812.0, 6),
            (1710.5, 7),
            (1613.0, 8),
            (1512.5, 9),
            (1412.5, 10),
            (1314.0, 11),
            (1214.0, 12),
            (1115.0, 13),
            (1015.5, 14),
            (919.0, 15),
        ),
        (("A", 3861.0), ("B", 4343.5), ("C", 4355.0)),
    ),
    # Fig. 6, circular total, 10^-50 cm^4 s.
    "fig6": Figure(
        12,
        "fig_ea1ce7340e25",
        1e-50,
        897,
        2678,
        3252,
        4542,
        (
            (3428.0, 9),
            (3597.0, 10),
            (3771.0, 11),
            (3942.5, 12),
            (4112.0, 13),
            (4285.0, 14),
            (4459.5, 15),
        ),
        (
            (3423.5, 9),
            (3592.0, 10),
            (3765.5, 11),
            (3936.0, 12),
            (4107.0, 13),
            (4278.5, 14),
            (4452.5, 15),
        ),
        (
            (2424.0, 0.5),
            (2170.0, 1.0),
            (1914.5, 1.5),
            (1660.0, 2.0),
            (1406.5, 2.5),
            (1154.0, 3.0),
        ),
        (
            (2425.0, 0.5),
            (2171.5, 1.0),
            (1915.0, 1.5),
            (1659.5, 2.0),
            (1404.5, 2.5),
            (1152.0, 3.0),
        ),
        (
            ("A", 3872.0),
            ("a", 3979.0),
            ("B", 4300.0),
            ("C", 4314.5),
            ("b", 4367.0),
            ("c", 4384.5),
        ),
    ),
    # Fig. 7, linear total, 10^-50 cm^4 s: the top 14 eV tick is under a curve.
    "fig7": Figure(
        13,
        "fig_167e2294aa4c",
        1e-50,
        887,
        2691,
        1055,
        2364,
        (
            (1234.5, 9),
            (1410.0, 10),
            (1583.5, 11),
            (1757.5, 12),
            (1934.0, 13),
            (2108.5, 14),
            (2284.5, 15),
        ),
        (
            (1224.5, 9),
            (1398.0, 10),
            (1572.5, 11),
            (1745.0, 12),
            (1921.5, 13),
            (2274.0, 15),
        ),
        (
            (2433.0, 0.5),
            (2175.5, 1.0),
            (1918.0, 1.5),
            (1660.5, 2.0),
            (1405.5, 2.5),
            (1147.0, 3.0),
        ),
        (
            (2430.5, 0.5),
            (2173.5, 1.0),
            (1915.0, 1.5),
            (1656.0, 2.0),
            (1399.0, 2.5),
            (1141.0, 3.0),
        ),
        (
            ("A", 1693.5),
            ("a", 1801.5),
            ("B", 2124.5),
            ("C", 2142.5),
            ("b", 2185.5),
            ("c", 2199.5),
        ),
    ),
}


@dataclass(frozen=True)
class Segment:
    """A stretch of one curve: its file, line style, a seed pixel on it, its span (eV).

    The style is checked over ``checked`` (eV) if given, a stretch with dots and dashes
    enough to tell it, or clear of a curve crossing; else over the span.
    """

    file: str
    figure: str
    style: str
    seed: tuple[int, int]
    span: tuple[float, float]
    checked: tuple[float, float] | None = None


def _file(name: str) -> str:
    """Return a dataset's data file."""
    return f"Ar/starace_1987__{name}_Ar.csv"


# Each seed is a pixel on the curve, identified on the image by the curves' order at
# that column and checked by its line style (:func:`style`); each span is where the
# curve is followed and sampled (the steep rises into the resonances and the first
# tenths of an eV above 8 eV, where the curves meet each other or the ticks, are out).
SEGMENTS = (
    Segment(
        _file("fig2_hf_f1D_linear"),
        "fig2",
        "dash-dot",
        (3299, 4411),
        (8.5, 14.0),
        (8.5, 12.9),
    ),
    Segment(
        _file("fig2_gsc_f1D_linear"),
        "fig2",
        "dash-double-dot",
        (3299, 4485),
        (8.5, 13.0),
    ),
    Segment(_file("fig2_iic_f1D_linear"), "fig2", "solid", (3299, 4464), (8.5, 11.6)),
    Segment(_file("fig2_esi_f1D_linear"), "fig2", "dashed", (3299, 4518), (8.5, 11.6)),
    Segment(_file("fig2_esi_f1D_linear"), "fig2", "dashed", (4220, 4678), (13.1, 14.0)),
    Segment(_file("fig3_hf_p1D_linear"), "fig3", "dash-dot", (1240, 2433), (8.1, 10.4)),
    Segment(
        _file("fig3_hf_p1D_linear"), "fig3", "dash-dot", (2030, 2594), (13.0, 14.0)
    ),
    Segment(_file("fig3_iic_p1D_linear"), "fig3", "solid", (1240, 2485), (8.1, 10.6)),
    Segment(_file("fig3_esi_p1D_linear"), "fig3", "dashed", (1240, 2597), (8.1, 10.0)),
    Segment(_file("fig3_esi_p1D_linear"), "fig3", "dashed", (2030, 2508), (13.0, 14.0)),
    Segment(_file("fig4_hf_p1S_linear"), "fig4", "dash-dot", (1310, 5161), (8.2, 10.8)),
    Segment(
        _file("fig4_gsc_p1S_linear"),
        "fig4",
        "dash-double-dot",
        (1310, 5219),
        (8.2, 10.8),
    ),
    Segment(_file("fig4_iic_p1S_linear"), "fig4", "solid", (1310, 5259), (8.2, 10.9)),
    Segment(_file("fig4_esi_p1S_linear"), "fig4", "dashed", (1350, 5487), (8.1, 11.2)),
    Segment(_file("fig4_esi_p1S_linear"), "fig4", "dashed", (2060, 5793), (13.5, 14.0)),
    Segment(
        _file("fig4_hf_p1S_linear"), "fig4", "dash-dot", (2060, 5618), (13.5, 14.05)
    ),
    Segment(
        _file("fig4_iic_p1S_linear"),
        "fig4",
        "solid",
        (2060, 5685),
        (13.6, 14.05),
        (13.0, 14.05),
    ),
    Segment(
        _file("fig5_iic_shifted_f1D_linear"), "fig5", "solid", (3340, 1792), (8.6, 10.3)
    ),
    Segment(
        _file("fig5_iic_shifted_f1D_linear"),
        "fig5",
        "solid",
        (4010, 1635),
        (11.9, 13.5),
    ),
    Segment(
        _file("fig5_esi_shifted_f1D_linear"),
        "fig5",
        "dashed",
        (3380, 1867),
        (8.4, 11.4),
    ),
    Segment(
        _file("fig5_esi_shifted_f1D_linear"),
        "fig5",
        "dashed",
        (4040, 1879),
        (11.9, 13.4),
    ),
    Segment(_file("fig6_iic_circular"), "fig6", "solid", (3300, 2232), (8.1, 11.1)),
    Segment(_file("fig6_iic_circular"), "fig6", "solid", (4200, 2162), (13.05, 14.05)),
    Segment(_file("fig6_esi_circular"), "fig6", "dashed", (3300, 2286), (8.1, 11.2)),
    Segment(_file("fig6_esi_circular"), "fig6", "dashed", (4210, 2346), (13.2, 14.0)),
    Segment(_file("fig7_iic_linear"), "fig7", "solid", (1145, 2030), (8.2, 10.8)),
    Segment(_file("fig7_iic_linear"), "fig7", "solid", (2020, 2256), (13.0, 14.05)),
    Segment(_file("fig7_esi_linear"), "fig7", "dashed", (1145, 2134), (8.1, 9.5)),
    Segment(_file("fig7_esi_linear"), "fig7", "dashed", (2050, 2427), (13.4, 14.0)),
)

# Samples (eV) that cannot be read cleanly (another curve within CLEAR pixels, or the
# line too steep), by file.
LEFT_OUT: dict[str, tuple[float, ...]] = {
    _file("fig2_hf_f1D_linear"): (13.0, 13.1, 13.2),
    _file("fig2_iic_f1D_linear"): (8.5, 11.6),
    _file("fig2_esi_f1D_linear"): (11.6, 13.1, 13.2, 13.3, 13.4, 13.5),
    _file("fig3_hf_p1D_linear"): (8.1, 13.0, 14.0),
    _file("fig3_iic_p1D_linear"): (8.1, 8.2, 8.3, 10.4, 10.5, 10.6),
    _file("fig3_esi_p1D_linear"): (
        8.1,
        8.2,
        8.3,
        9.8,
        9.9,
        10.0,
        13.0,
        13.1,
        13.2,
        13.3,
        13.4,
        13.9,
        14.0,
    ),
    _file("fig4_hf_p1S_linear"): (10.8,),
    _file("fig4_gsc_p1S_linear"): (8.2, 8.3, 10.8),
    _file("fig4_iic_p1S_linear"): (8.2, 8.3, 10.9),
    _file("fig4_esi_p1S_linear"): (8.1, 8.2, 8.3, 11.2, 13.5),
    _file("fig5_iic_shifted_f1D_linear"): (8.6, 8.7, 10.0, 10.1, 10.2, 10.3, 13.5),
    _file("fig6_iic_circular"): (
        8.7,
        8.8,
        8.9,
        9.0,
        9.6,
        9.7,
        9.8,
        9.9,
        10.0,
        10.1,
        10.2,
        11.1,
        13.0,
        13.1,
        14.0,
    ),
    _file("fig6_esi_circular"): (11.2, 13.2),
    _file("fig7_iic_linear"): (8.2, 8.3, 13.0, 14.0),
    _file("fig7_esi_linear"): (9.5, 13.4),
}

# The tracing (:func:`trace`), pixels: the centres a prediction is fitted through, the
# fewest for a fit, how far a run may lie from the prediction, the widest run that
# moves the trace.
FIT = 30
MIN_FIT = 6
STEP = 6
MAX_WIDTH = 12
# The line style (:func:`style`), columns: the longest dot, the longest dash of a
# dashed line, the shortest solid stretch, how near a resonance line a column says
# nothing of the style.
DOT = 10
DASH = 80
MIN_SOLID = 160
NEAR_LINE = 8
# The reading (:func:`_read`), pixels: the clearance from any other dark pixel, a line's
# width and slack, the columns either side that must be clear too, the widest gap
# interpolated across, the most crowded columns in it, the steepest slope (rows per
# column), how near the frame a tick can be taken for a curve; the samples' step, eV.
CLEAR = 6
LINE = 7
SLACK = 3
WINDOW = 10
GAP = 40
CROWDED = 8
MAX_SLOPE = 3.0
TICK_ZONE = 70
GRID = 0.1

Run = tuple[int, int]


def _spans(mask: np.ndarray) -> list[Run]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


class Calibration:
    """Map a pixel (column, row) to (photon energy, value) through a figure's ticks.

    The photon energy is interpolated between the bottom and top ticks' scales by the
    row, the value between the left and right ticks' scales by the column: the scan's
    small rotation and stretch.
    """

    def __init__(self, figure: Figure) -> None:
        self.figure = figure
        axes = {
            name: Axis(
                tuple(p for p, _ in ticks),
                tuple(float(v) for _, v in ticks),
                tolerance=TICK_TOLERANCE,
            )
            for name, ticks in (
                ("bottom", figure.x_bottom),
                ("top", figure.x_top),
                ("left", figure.y_left),
                ("right", figure.y_right),
            )
        }
        for axis in axes.values():
            axis.check()
        # The scales are linear: value = slope * pixel + intercept.
        self.scales = {
            name: (float(axis(1.0)) - float(axis(0.0)), float(axis(0.0)))
            for name, axis in axes.items()
        }
        self.rows = (figure.top + DEPTH, figure.bottom - DEPTH)
        self.columns = (figure.left + DEPTH, figure.right - DEPTH)

    def _at(self, side: str, pixel: float) -> float:
        """Return one side's scale at a pixel."""
        slope, intercept = self.scales[side]
        return slope * pixel + intercept

    def energy(self, column: float, row: float) -> float:
        """Return the photon energy (eV) at a pixel."""
        top, bottom = self.rows
        t = (row - top) / (bottom - top)
        upper = self._at("top", column)
        return upper + t * (self._at("bottom", column) - upper)

    def value(self, column: float, row: float) -> float:
        """Return the y value, in the axis's printed unit, at a pixel."""
        left, right = self.columns
        s = (column - left) / (right - left)
        west = self._at("left", row)
        return west + s * (self._at("right", row) - west)

    def column(self, energy: float, row: float) -> float:
        """Return the (fractional) column of a photon energy at a row."""
        a = self.energy(0.0, row)
        return (energy - a) / (self.energy(1.0, row) - a)

    def row(self, value: float, column: float) -> float:
        """Return the (fractional) row of a value at a column."""
        a = self.value(column, 0.0)
        return (value - a) / (self.value(column, 1.0) - a)


def _check_ticks(dark: np.ndarray, figure: Figure) -> None:
    """Raise unless every listed tick and resonance line, and the frame, are found."""
    rows = slice(figure.top + 10, figure.bottom - 10)
    columns = slice(figure.left + 10, figure.right - 10)
    lines = tuple((p, 0.0) for _, p in figure.lines)
    probes = (
        ("bottom", dark[figure.bottom - DEPTH, columns], figure.left, figure.x_bottom),
        ("top", dark[figure.top + DEPTH, columns], figure.left, figure.x_top),
        ("left", dark[rows, figure.left + DEPTH], figure.top, figure.y_left),
        ("right", dark[rows, figure.right - DEPTH], figure.top, figure.y_right),
        ("line", dark[figure.bottom - LINE_DEPTH, columns], figure.left, lines),
    )
    for name, pixels, offset, ticks in probes:
        found = [offset + 10 + (a + b) / 2 for a, b in _spans(pixels)]
        for position, _ in ticks:
            if not any(abs(f - position) <= TICK_FOUND for f in found):
                msg = f"p. {figure.page}: no {name} tick at {position}"
                raise ValueError(msg)
    for row in (figure.top, figure.bottom):
        band = dark[row - 8 : row + 9, figure.left + 20 : figure.right - 20]
        if band.any(axis=0).mean() < 0.99:
            msg = f"p. {figure.page}: no frame line at row {row}"
            raise ValueError(msg)


def trace(
    dark: np.ndarray, seed: tuple[int, int], columns: tuple[int, int]
) -> dict[int, Run | None]:
    """Return the run a curve passes through at each column, traced from a seed.

    From the seed forward, then back (starting from the first ``FIT`` centres traced
    forward), column by column: the position predicted by a parabola through the last
    ``FIT`` traced centres (a straight line through fewer); the run whose nearest pixel
    is within ``STEP`` of it is the curve's (the nearest, if several), and one no wider
    than ``MAX_WIDTH`` moves the trace. Columns with no such run (a gap of a dashed
    line) are None.
    """
    out: dict[int, Run | None] = {}
    first, last = columns
    start = [(float(seed[0]), float(seed[1]))]
    for direction in (1, -1):
        history = list(start)
        column = seed[0]
        while first <= column <= last:
            recent = history[-FIT:]
            if len(recent) >= MIN_FIT:
                xs = np.array([x for x, _ in recent])
                ys = np.array([y for _, y in recent])
                degree = 2 if len(recent) == FIT else 1
                predicted = float(np.polyval(np.polyfit(xs, ys, degree), column))
            else:
                predicted = recent[-1][1]
            low = int(predicted) - STEP - 3 * MAX_WIDTH
            pixels = dark[low : int(predicted) + STEP + 3 * MAX_WIDTH + 1, column]
            near = [
                (abs((a + b) / 2 - predicted), (a, b))
                for a, b in ((low + a, low + b) for a, b in _spans(pixels))
                if abs(min(max(predicted, a), b) - predicted) <= STEP
            ]
            out[column] = min(near)[1] if near else None
            run = out[column]
            if run is not None and run[1] - run[0] + 1 <= MAX_WIDTH:
                history.append((float(column), _centre(run)))
            column += direction
        # The way back starts from the line's slope at the seed, as traced forward.
        start = history[:FIT][::-1]
    return out


def _pieces(found: list[bool]) -> list[tuple[bool, int]]:
    """Return the runs of a sequence of flags, as (flag, length)."""
    pieces: list[tuple[bool, int]] = []
    for flag in found:
        if pieces and pieces[-1][0] == flag:
            pieces[-1] = (flag, pieces[-1][1] + 1)
        else:
            pieces.append((flag, 1))
    return pieces


def style(found: list[bool]) -> str:
    """Return a traced line's style from which of its columns have a run.

    Before its first and after its last run the line is not drawn. ``solid``: no gap,
    over at least ``MIN_SOLID`` columns (longer than any dash). Otherwise, with two
    dashes (runs longer than ``DOT``) or more, by the number of dots (shorter runs)
    between consecutive dashes, which must be the same between every two: 0 ``dashed``
    (with no dash longer than ``DASH``), 1 ``dash-dot``, 2 ``dash-double-dot``. Anything
    else is ``unknown``.
    """
    pieces = _pieces(found)
    while pieces and not pieces[0][0]:
        pieces.pop(0)
    while pieces and not pieces[-1][0]:
        pieces.pop()
    if all(flag for flag, _ in pieces):
        return "solid" if sum(n for _, n in pieces) >= MIN_SOLID else "unknown"
    dashes = [i for i, (flag, n) in enumerate(pieces) if flag and n > DOT]
    counts = {
        sum(1 for flag, n in pieces[a + 1 : b] if flag and n <= DOT)
        for a, b in itertools.pairwise(dashes)
    }
    names = {0: "dashed", 1: "dash-dot", 2: "dash-double-dot"}
    if len(counts) != 1 or not counts <= set(names):
        return "unknown"
    (count,) = counts
    if count == 0 and max(n for flag, n in pieces if flag) > DASH:
        return "unknown"
    return names[count]


def _clear(dark: np.ndarray, column: int, run: Run, slope: float) -> bool:
    """Return True if a run is thin and no other dark pixel is within ``CLEAR`` of it.

    Thin: no wider than a line ``LINE`` pixels wide crosses the column at the curve's
    slope (rows per column), plus ``SLACK``.
    """
    a, b = run
    if b - a + 1 > LINE * np.hypot(1.0, slope) + SLACK:
        return False
    above = dark[a - CLEAR : a, column]
    below = dark[b + 1 : b + 1 + CLEAR, column]
    return not (above.any() or below.any())


def _centre(run: Run) -> float:
    """Return a run's centre."""
    return (run[0] + run[1]) / 2


@dataclass(frozen=True)
class Reading:
    """A value read from a figure: the sample, the pixel, the energy and the value."""

    wanted: float
    column: int
    row: float
    energy: float
    value: float
    printed: str


class _Line:
    """A traced line's usable runs, and the columns where nothing crowds it."""

    def __init__(
        self, dark: np.ndarray, figure: Figure, line: dict[int, Run | None]
    ) -> None:
        ticks = [r for r, _ in figure.y_left + figure.y_right]

        def by_tick(c: int, row: float) -> bool:
            edge = c - figure.left < TICK_ZONE or figure.right - c < TICK_ZONE
            return edge and any(abs(row - r) <= MAX_WIDTH for r in ticks)

        def slope(c: int) -> float:
            before, after = line.get(c - 3), line.get(c + 3)
            if before is None or after is None:
                return 0.0
            return (_centre(after) - _centre(before)) / 6

        clear = {
            c: _clear(dark, c, run, slope(c)) and not by_tick(c, _centre(run))
            for c, run in line.items()
            if run is not None
        }
        # Where a resonance line crosses, the curve's run is not clear, but no other
        # curve is near.
        lines = [p for _, p in figure.lines]
        self.open = {
            c: run is None or clear[c] or any(abs(c - p) <= NEAR_LINE for p in lines)
            for c, run in line.items()
        }
        self.usable = {
            c: _centre(run)
            for c, run in line.items()
            if run is not None
            and clear[c]
            and all(self.open.get(d, True) for d in range(c - WINDOW, c + WINDOW + 1))
        }


def _read(calibration: Calibration, line: _Line, wanted: float) -> Reading | None:
    """Return the curve's value at a photon energy, or None if it cannot be read.

    A run is usable if it, and the line's every run within ``WINDOW`` columns of it, is
    clear (:func:`_clear`; so a gap in a neighbouring dashed line does not make a
    crowded stretch look clear) and not where a tick mark could be taken for it (within
    ``TICK_ZONE`` columns of the left or right frame line and ``MAX_WIDTH`` rows of a
    tick). The value is the usable run's centre at the sample's column, or the linear
    interpolation between the nearest usable columns (at most ``GAP`` apart, with at
    most ``CROWDED`` crowded columns between); none where the line is steeper than
    ``MAX_SLOPE``.
    """
    usable = line.usable
    if not usable:
        return None
    nearest = min(usable, key=lambda c: abs(calibration.energy(c, usable[c]) - wanted))
    column = round(calibration.column(wanted, usable[nearest]))
    window = [c for c in range(column - 8, column + 9) if c in usable]
    if len(window) >= 4:
        slope = abs(float(np.polyfit(window, [usable[c] for c in window], 1)[0]))
        if slope > MAX_SLOPE:
            return None
    if column in usable:
        row = usable[column]
        printed = f"pixel ({column}, {row:.1f})"
    else:
        left = next((c for c in range(column - 1, column - GAP, -1) if c in usable), 0)
        right = next((c for c in range(column + 1, column + GAP) if c in usable), 0)
        if not (left and right) or right - left > GAP:
            return None
        if sum(1 for c in range(left + 1, right) if not line.open[c]) > CROWDED:
            return None
        weight = (column - left) / (right - left)
        row = (1 - weight) * usable[left] + weight * usable[right]
        printed = (
            f"between pixels ({left}, {usable[left]:.1f}) and "
            f"({right}, {usable[right]:.1f})"
        )
    return Reading(
        wanted,
        column,
        row,
        calibration.energy(column, row),
        calibration.value(column, row),
        printed,
    )


def _grid(span: tuple[float, float]) -> list[float]:
    """Return the samples, every ``GRID`` eV within a span."""
    first, last = (round(x / GRID) for x in span)
    return [round(k * GRID, 2) for k in range(first, last + 1)]


def _follow(
    dark: np.ndarray, calibration: Calibration, segment: Segment
) -> dict[int, Run | None]:
    """Return a segment's traced line; raise if its style is not the one written."""
    figure = calibration.figure
    low, high = segment.checked or segment.span
    row = segment.seed[1]
    first = int(calibration.column(min(low, segment.span[0]), row))
    last = int(np.ceil(calibration.column(max(high, segment.span[1]), row)))
    line = trace(dark, segment.seed, (first, last))
    # The resonance lines cross the curves: near one, a column keeps the previous
    # column's state.
    found: list[bool] = []
    for c in range(
        int(calibration.column(low, row)),
        int(np.ceil(calibration.column(high, row))) + 1,
    ):
        near_line = any(abs(c - p) <= NEAR_LINE for _, p in figure.lines)
        found.append(found[-1] if near_line and found else line[c] is not None)
    drawn = style(found)
    if drawn != segment.style:
        msg = f"{segment.file} {segment.span}: a {drawn} line, not {segment.style}"
        raise ValueError(msg)
    return line


def read_all() -> dict[str, list[Reading]]:
    """Return every segment's readings, by data file, in photon-energy order.

    Raises if a segment's line style is not the one written, or the samples left out
    are not exactly those of :data:`LEFT_OUT`.
    """
    pages: dict[int, np.ndarray] = {}
    calibrations: dict[str, Calibration] = {}
    out: dict[str, list[Reading]] = {}
    missing: dict[str, list[float]] = {}
    for segment in SEGMENTS:
        figure = FIGURES[segment.figure]
        if figure.page not in pages:
            pages[figure.page] = image(DIRECTORY, figure.page) < DARK
        dark = pages[figure.page]
        if segment.figure not in calibrations:
            _check_ticks(dark, figure)
            calibrations[segment.figure] = Calibration(figure)
        calibration = calibrations[segment.figure]
        line = _Line(dark, figure, _follow(dark, calibration, segment))
        for wanted in _grid(segment.span):
            reading = _read(calibration, line, wanted)
            if reading is None:
                missing.setdefault(segment.file, []).append(wanted)
            else:
                out.setdefault(segment.file, []).append(reading)
    for name in set(missing) | set(LEFT_OUT):
        if tuple(missing.get(name, ())) != LEFT_OUT.get(name, ()):
            msg = f"{name}: left out {missing.get(name)}, not {LEFT_OUT.get(name)}"
            raise ValueError(msg)
    return {name: sorted(rows, key=lambda r: r.wanted) for name, rows in out.items()}


# Table III (p. 11/1715): the resonance energies the lines A, B, C (experiment) and
# a, b, c (HF) are drawn at; the drawn lines lie within LINE_MARGIN (eV) of them.
TABLE_III = "tbl_3de3df541ce4"
LINE_MARGIN = 0.08
# Table II (p. 9/1713): f_q of Eq. 23.
TABLE_II = "tbl_5e868ff91477"
# Eq. 24 (p. 9/1713): 'L^4 T (a.u.) = 1.8967 x 10^-50 cm^4 sec'; a0 (cm) and the
# atomic unit of time (s), CODATA 2018.
EQ_24 = 1.8967e-50
BOHR_CM = 5.29177210903e-9
TIME_AU_S = 2.4188843265857e-17
# How far the totals may stray from the sums of the channels read from Figs. 2-4.
SUMS = 0.03
# The comparison with McKenna and van der Hart 2004 below this photon energy (eV), off
# the rise into the 4s resonance.
MCKENNA_BELOW = 10.4
# The text's 'region between threshold and the ... 4s(1P) intermediate-state
# resonance' (p. 10/1714), for the reported ratios between approximations.
BELOW_4S = 11.5
CHANNELS = ("fig2_{}_f1D_linear", "fig3_{}_p1D_linear", "fig4_{}_p1S_linear")


def _resonances() -> dict[str, float]:
    """Return Table III's energies (eV) by line label, from the table's cells."""
    path = LIBRARY / DIRECTORY / "tables" / f"{TABLE_III}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    grid = {(c["row"], c["column"]): c["text"] for c in cells}
    labels = {"4s": ("A", "a"), "5s": ("B", "b"), "3d": ("C", "c")}
    out: dict[str, float] = {}
    for row in (2, 3, 4):
        key = next(k for k in labels if k in grid[(row, 0)])
        experiment, hartree_fock = labels[key]
        out[experiment] = float(grid[(row, 1)])
        out[hartree_fock] = float(grid[(row, 2)])
    return out


def _fraction(cell: str) -> float:
    r"""Return a Table II cell's number: ``\frac{a}{b}`` or a plain number."""
    text = cell.strip()
    if text.startswith("\\frac{"):
        numerator, denominator = (
            text.removeprefix("\\frac{").removesuffix("}").split("}{")
        )
        return float(numerator) / float(denominator)
    return float(text)


def _values(readings: dict[str, list[Reading]]) -> dict[str, dict[float, float]]:
    """Return each dataset's values, cm^4 s, by sample, keyed by its name."""
    units = {segment.file: FIGURES[segment.figure].unit for segment in SEGMENTS}
    return {
        name.removeprefix("Ar/starace_1987__").removesuffix("_Ar.csv"): {
            r.wanted: r.value * units[name] for r in rows
        }
        for name, rows in readings.items()
    }


def _ratios(
    top: dict[float, float], bottom: dict[float, float], below: float = math.inf
) -> tuple[float, float, float, float]:
    """Return the first and last shared samples below an energy, and the ratios."""
    shared = sorted(e for e in top if e in bottom and e < below)
    ratios = [top[e] / bottom[e] for e in shared]
    return shared[0], shared[-1], min(ratios), max(ratios)


def _sum(parts: list[dict[float, float]], factor: float = 1.0) -> dict[float, float]:
    """Return the sum of datasets times a factor where all are read."""
    shared = set(parts[0]).intersection(*parts[1:])
    return {e: factor * sum(p[e] for p in parts) for e in shared}


def checks(readings: dict[str, list[Reading]]) -> dict[str, object]:
    """Return the checks of the reading; raise if one fails its margin.

    The resonance lines against Table III (:data:`LINE_MARGIN`); Table II's 3/2 between
    circular and linear 1D channels; Eq. 24 against CODATA; Fig. 6 against 3/2 times
    Figs. 2 + 3 and Fig. 7 against Figs. 2 + 3 + 4 where all are read (:data:`SUMS`).
    Reported only: the ratios between the approximations the text describes (p.
    10/1714), Fig. 7 against McKenna and van der Hart 2004, and the shifted Fig. 5
    against Fig. 2.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    table_iii = _resonances()
    for name, figure in FIGURES.items():
        calibration = Calibration(figure)
        for label, column in figure.lines:
            drawn = calibration.energy(column, figure.bottom - LINE_DEPTH)
            out[f"{name} line {label}"] = (round(drawn, 3), table_iii[label])
            if abs(drawn - table_iii[label]) > LINE_MARGIN:
                failures.append(f"{name}: line {label} at {drawn:.3f} eV")
    f_q = {
        (q.strip(), term.strip()): _fraction(value)
        for q, term, value in table(DIRECTORY, TABLE_II)[1:]
    }
    ratio = f_q[("\u00b11", "^{1}D")] / f_q[("0", "^{1}D")]
    out["Table II circular/linear 1D"] = ratio
    if abs(ratio - 1.5) > 1e-12 or f_q[("\u00b11", "^{1}S")] != 0:
        failures.append("Table II: circular is not 3/2 linear for 1D")
    unit = BOHR_CM**4 * TIME_AU_S
    out["Eq. 24 against CODATA"] = (EQ_24, unit)
    if abs(EQ_24 / unit - 1) > 1e-4:
        failures.append(f"Eq. 24: {EQ_24} against {unit:.5e}")
    values = _values(readings)
    for method in ("iic", "esi"):
        f1d, p1d, p1s = (values[c.format(method)] for c in CHANNELS)
        for name, parts, factor in (
            (f"fig6_{method}_circular", [f1d, p1d], 1.5),
            (f"fig7_{method}_linear", [f1d, p1d, p1s], 1.0),
        ):
            first, last, low, high = _ratios(values[name], _sum(parts, factor))
            out[f"{name} / the channels' sum"] = (first, last, low, high)
            if max(abs(low - 1), abs(high - 1)) > SUMS:
                failures.append(f"{name}: {low:.3f} to {high:.3f} of the channels")
    for channel in CHANNELS:
        hf, iic, esi = (values[channel.format(m)] for m in ("hf", "iic", "esi"))
        out[channel.format("iic/hf")] = _ratios(iic, hf, BELOW_4S)
        out[channel.format("esi/iic")] = _ratios(esi, iic, BELOW_4S)
    hf = _sum([values[c.format("hf")] for c in CHANNELS])
    for method in ("iic", "esi"):
        total = _sum([values[c.format(method)] for c in CHANNELS])
        out[f"linear sum {method}/hf"] = _ratios(total, hf, BELOW_4S)
    # McKenna and van der Hart 2004 (R-matrix Floquet, linear; TC-856).
    floquet = next(d for d in datasets() if d.key == "mckenna_2004/fig5_sigma2_Ar")
    energy, sigma = floquet.abscissa_values(), floquet.scaled_values()
    for method in ("iic", "esi"):
        total = values[f"fig7_{method}_linear"]
        mckenna = {e: float(np.interp(e, energy, sigma)) for e in total}
        out[f"fig7 {method} / mckenna_2004"] = _ratios(total, mckenna, MCKENNA_BELOW)
    out["fig5 shifted / fig2 IIC"] = _ratios(
        values["fig5_iic_shifted_f1D_linear"], values["fig2_iic_f1D_linear"], BELOW_4S
    )
    if failures:
        raise ValueError("; ".join(failures))
    return out


def _points(readings: dict[str, list[Reading]]) -> dict[str, list[Point]]:
    """Return the data files' points: cm^4 s against photon energy (eV)."""
    units = {segment.file: FIGURES[segment.figure].unit for segment in SEGMENTS}
    return {
        name: [
            Point(
                f"{r.wanted:.2f}",
                f"{r.energy:.3f}",
                f"{r.value * units[name]:.3e}",
                "",
                r.printed,
            )
            for r in rows
        ]
        for name, rows in readings.items()
    }


def extract() -> dict[str, list[Point]]:
    """Return the cross sections read from Figs. 2-7, cm^4 s against photon energy."""
    return _points(read_all())


def main() -> None:
    """Write the data files and print the checks."""
    readings = read_all()
    write(_points(readings))
    for label, result in checks(readings).items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()

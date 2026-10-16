# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pan, Gao and Starace (1990): two-photon ionisation cross sections of Ar.

The paper prints its generalised two-photon ionisation cross sections of the Ar 3p
subshell below the one-photon threshold only in figures (variationally stable
perturbation theory with electron correlations, length and velocity forms); no table,
no supplement. The PDF is a scan: every page is one embedded image, 5060 x 6898
pixels (about 595 dpi). The figures are read from those native pixels (data rule 8,
TC-699). Journal page = PDF page + 6270.

* **Fig. 1(b)** (p. 9/6279, ``fig_41ffe3dfe10f``): the 3p5 epsilon-p 1S channel, linear
  polarisation, y in units of 10^-51 cm^4 s, 0 to 10: length (solid) and velocity
  (dashed) with all the processes of Table I; length (dash-dot) and velocity
  (dash-double-dot) without the intermediate-state shake-up, process (4).
* **Fig. 2** (p. 10/6280, ``fig_08fff07d737e``): 3p5 epsilon-p 1D, linear, 0 to 10e-51;
  length with all processes (solid). The HF length curve (dashed) is not read.
* **Fig. 3** (p. 10/6280, ``fig_7ef7a8747556``): 3p5 epsilon-f 1D, linear, 0 to 10e-51;
  length with all processes (solid). The HF (dash-dot) and processes (1)-(3) (dashed)
  curves are not read.
* **Fig. 4** (p. 10/6280, ``fig_6af979c5d055``): the linear total, 0 to 2e-50; (a) HF
  length (solid) and velocity (dashed); (b) without the shake-up, length (solid) and
  velocity (dashed); (c) all processes, length (solid) and velocity (dashed). Panel
  (a)'s processes (1)-(3) curves (dash-dot, dash-double-dot) are not read.
* **Fig. 5** (p. 11/6281, ``fig_cc071c4e5baf``): the circular total, 0 to 3.5e-50; HF
  length (dash-dot) and velocity (dash-double-dot), all processes length (solid) and
  velocity (dashed).
* **Fig. 6(a)** (p. 11/6281, ``fig_4360febfcc00``): the linear total, length, all
  processes (solid), drawn again beside other authors' curves: read for the checks
  only (not a data file). The other authors' curves are not read.

How the images are read (as starace_1987, the same authors' group and style, TC-878):

* the axes: the tick marks, found ``DEPTH`` (or ``DEEP``) pixels inside each frame line
  and checked against the positions written in :data:`FIGURES` (each within
  ``TICK_FOUND``); the frame lines at the labelled corners (0 at the bottom, the top
  label where it stands at the corner) are measured at the probe columns and added to
  the y ticks. Each side's ticks lie on one linear scale within ``TICK_TOLERANCE``
  pixels (:meth:`_figure.Axis.check`). Ticks a curve runs through are left out of the
  lists. The scan is rotated by about 0.2 degrees: the photon energy is interpolated
  between the bottom and top ticks' scales by the row, the value between the left and
  right ticks' scales by the column (:class:`Calibration`).
* the curves: each stretch of a curve (:data:`SEGMENTS`) is traced column by column from
  a seed pixel on it (:func:`trace`), and its line style is checked from the trace
  (:func:`style`): a solid line has no gap; a dashed one has gaps with nothing in them;
  a dash-dot one one dot per gap, a dash-double-dot one two. This is what tells the
  curves apart.
* the values: every 0.1 eV within a stretch, the centre of the curve's run at that
  photon energy's column where the run is clear: as thin as one line at its slope, no
  other dark pixel within ``CLEAR`` pixels (about a line width) of it, and so for the
  line's runs ``WINDOW`` columns either side, away from the tick marks (:func:`_read`);
  across a gap of the dashes, the linear interpolation between the nearest clear
  columns. Where another curve comes closer, or the line is steeper than ``MAX_SLOPE``
  rows per column, the sample is left out: :data:`LEFT_OUT` lists them. Where two
  curves of a figure are read within ``LINE`` rows of each other, the reading cannot
  tell whose line it is, and the sample is left out of both (:data:`TOGETHER`). The
  script fails if either list changes. The ``printed`` cell of each value gives its
  pixel (column, row).

One pixel is 0.0046 eV in x; in y 0.0097e-51 cm^4 s (Figs. 1-3), 0.0019e-50 (Fig. 4),
0.0033e-50 (Fig. 5). The lines are 5 to 8 pixels wide.

**Checks** (:func:`checks`, run by :func:`main` and :func:`extract`; it fails on a
miss): Fig. 4(c) against Figs. 1(b) + 2 + 3 (the total is the sum of the channels, Eq.
49); Fig. 5 against 3/2 times Figs. 2 + 3 (Table II's f_q: the circular total has no 1S
part and 3/2 the linear 1D channels); Fig. 6(a) against Fig. 4(c) (the same curve drawn
twice); Table II's f_q; Eq. 54's unit against CODATA. Reported: the pixel sizes; the
shake-up's
reduction of the p 1S channel and of the total ('by approximately 0.1 x 10^-50', p.
12/6282) from Figs. 1(b) and 4(b, c); length over velocity; Fig. 4(c) against McKenna
and van der Hart 2004 (TC-856), Petrov et al. 2016 (TC-877) and Starace and Jiang 1987
(TC-878) (:func:`comparisons`).
"""

import itertools
import math
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, table, write

from optmat.datasets import DATA_DIRECTORY, Dataset, load_source

DIRECTORY = "Pan_1990_TwoPhotonIonization"
SOURCE = "pan_1990b"
DARK = 128


@dataclass(frozen=True)
class Figure:
    """One plot on a scanned page: its frame and tick marks in image pixels.

    The frame lines' rows and columns are those of their centres near the frame's
    middle. Ticks are (pixel, printed value): the bottom and top ones (columns, eV) as
    found ``DEPTH`` or ``DEEP`` pixels inside the bottom and top frame lines, the left
    and right ones (rows, the y axis's printed number) inside the left and right ones.
    ``corners`` are the y values of the bottom and top frame lines where a label stands
    at the corner (None where it does not).
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
    corners: tuple[float | None, float | None]


DEPTH = 10
DEEP = 20
TICK_TOLERANCE = 4.5
TICK_FOUND = 1.5
# The frame line's centre at a probe column is searched within this many rows of the
# frame's row near the middle.
FRAME_SEARCH = 12

FIGURES = {
    # Fig. 1(b), epsilon-p 1S, 10^-51 cm^4 s: the bottom 12 and 13 eV ticks, the top 12
    # eV one, the left 6 and the right 7-9 are under curves.
    "fig1b": Figure(
        9,
        "fig_41ffe3dfe10f",
        1e-51,
        4223,
        5256,
        2956,
        4458,
        ((3171.5, 9), (3384.5, 10), (3601.0, 11), (4245.5, 14)),
        ((3172.0, 9), (3384.5, 10), (3602.0, 11), (4030.5, 13), (4244.0, 14)),
        (
            (4321.0, 9),
            (4425.5, 8),
            (4530.5, 7),
            (4739.0, 5),
            (4843.5, 4),
            (4945.0, 3),
            (5048.5, 2),
            (5154.0, 1),
        ),
        (
            (4638.0, 6),
            (4743.0, 5),
            (4845.5, 4),
            (4947.0, 3),
            (5051.0, 2),
            (5157.0, 1),
        ),
        (0.0, 10.0),
    ),
    # Fig. 2, epsilon-p 1D, 10^-51 cm^4 s: the bottom 14 eV tick and the right 3-5 are
    # under curves.
    "fig2": Figure(
        10,
        "fig_08fff07d737e",
        1e-51,
        678,
        1706,
        786,
        2296,
        ((1007.0, 9), (1222.5, 10), (1440.0, 11), (1653.5, 12), (1867.5, 13)),
        (
            (1003.0, 9),
            (1218.0, 10),
            (1434.5, 11),
            (1649.5, 12),
            (1863.0, 13),
            (2079.5, 14),
        ),
        (
            (785.0, 9),
            (889.0, 8),
            (990.5, 7),
            (1094.5, 6),
            (1193.5, 5),
            (1298.5, 4),
            (1404.0, 3),
            (1505.5, 2),
            (1609.5, 1),
        ),
        (
            (778.5, 9),
            (884.0, 8),
            (985.0, 7),
            (1089.5, 6),
            (1501.5, 2),
            (1607.5, 1),
        ),
        (0.0, 10.0),
    ),
    # Fig. 3, epsilon-f 1D, 10^-51 cm^4 s: the top 12 and 13 eV ticks, the left 4 and
    # the right 4 and 5 are under curves.
    "fig3": Figure(
        10,
        "fig_7ef7a8747556",
        1e-51,
        4546,
        5566,
        756,
        2267,
        (
            (973.0, 9),
            (1190.0, 10),
            (1404.5, 11),
            (1620.0, 12),
            (1836.0, 13),
            (2050.5, 14),
        ),
        ((969.5, 9), (1185.5, 10), (1398.5, 11), (2042.0, 14)),
        (
            (4648.5, 9),
            (4752.5, 8),
            (4854.0, 7),
            (4957.5, 6),
            (5060.0, 5),
            (5264.5, 3),
            (5364.5, 2),
            (5470.0, 1),
        ),
        (
            (4638.0, 9),
            (4744.0, 8),
            (4845.5, 7),
            (4951.0, 6),
            (5261.0, 3),
            (5358.5, 2),
            (5463.0, 1),
        ),
        (0.0, 10.0),
    ),
    # Fig. 4(a), the linear total, 10^-50 cm^4 s: the right 1.0 tick is under a curve.
    "fig4a": Figure(
        10,
        "fig_6af979c5d055",
        1e-50,
        1254,
        2294,
        2904,
        4434,
        (
            (3123.5, 9),
            (3341.5, 10),
            (3560.5, 11),
            (3780.5, 12),
            (4000.0, 13),
            (4215.5, 14),
        ),
        (
            (3119.5, 9),
            (3336.5, 10),
            (3555.5, 11),
            (3775.0, 12),
            (3995.0, 13),
            (4209.5, 14),
        ),
        ((1517.0, 1.5), (1775.0, 1.0), (2039.0, 0.5)),
        ((1507.5, 1.5), (2030.0, 0.5)),
        (0.0, 2.0),
    ),
    # Fig. 4(b).
    "fig4b": Figure(
        10,
        "fig_6af979c5d055",
        1e-50,
        2725,
        3763,
        2928,
        4449,
        (
            (3141.5, 9),
            (3363.5, 10),
            (3582.0, 11),
            (3801.0, 12),
            (4016.0, 13),
            (4235.5, 14),
        ),
        (
            (3140.0, 9),
            (3361.5, 10),
            (3579.5, 11),
            (3797.5, 12),
            (4014.0, 13),
            (4231.5, 14),
        ),
        ((2983.5, 1.5), (3241.5, 1.0), (3507.5, 0.5)),
        ((2983.0, 1.5), (3240.5, 1.0), (3505.5, 0.5)),
        (0.0, 2.0),
    ),
    # Fig. 4(c).
    "fig4c": Figure(
        10,
        "fig_6af979c5d055",
        1e-50,
        4086,
        5130,
        2938,
        4470,
        (
            (3162.0, 9),
            (3381.0, 10),
            (3598.0, 11),
            (3819.5, 12),
            (4038.0, 13),
            (4255.0, 14),
        ),
        (
            (3157.5, 9),
            (3377.0, 10),
            (3592.0, 11),
            (3813.0, 12),
            (4032.0, 13),
            (4248.0, 14),
        ),
        ((4349.5, 1.5), (4609.0, 1.0), (4874.0, 0.5)),
        ((4340.5, 1.5), (4600.0, 1.0), (4868.0, 0.5)),
        (0.0, 2.0),
    ),
    # Fig. 5, the circular total, 10^-50 cm^4 s: the right 1.5 tick is under a curve;
    # the top frame line (3.5) carries no label.
    "fig5": Figure(
        11,
        "fig_cc071c4e5baf",
        1e-50,
        666,
        1724,
        731,
        2260,
        (
            (947.5, 9),
            (1168.0, 10),
            (1390.5, 11),
            (1610.0, 12),
            (1827.5, 13),
            (2047.0, 14),
        ),
        (
            (948.0, 9),
            (1168.0, 10),
            (1390.0, 11),
            (1608.0, 12),
            (1825.5, 13),
            (2044.0, 14),
        ),
        (
            (814.0, 3.0),
            (966.0, 2.5),
            (1117.0, 2.0),
            (1269.0, 1.5),
            (1424.0, 1.0),
            (1575.5, 0.5),
        ),
        ((811.5, 3.0), (965.0, 2.5), (1118.0, 2.0), (1421.0, 1.0), (1572.5, 0.5)),
        (0.0, None),
    ),
    # Fig. 6(a), the linear total beside other calculations, 10^-50 cm^4 s.
    "fig6a": Figure(
        11,
        "fig_4360febfcc00",
        1e-50,
        2810,
        3862,
        2980,
        4518,
        (
            (3202.5, 9),
            (3423.0, 10),
            (3645.5, 11),
            (3864.0, 12),
            (4086.0, 13),
            (4303.5, 14),
        ),
        (
            (3205.5, 9),
            (3425.5, 10),
            (3646.0, 11),
            (3864.5, 12),
            (4086.0, 13),
            (4304.5, 14),
        ),
        ((3073.0, 1.5), (3335.0, 1.0), (3600.0, 0.5)),
        ((3075.5, 1.5), (3334.5, 1.0), (3601.5, 0.5)),
        (0.0, 2.0),
    ),
}


@dataclass(frozen=True)
class Segment:
    """A stretch of one curve: its data file, line style, a seed pixel, its span (eV).

    The style is checked over ``checked`` (eV) if given, a stretch with dots and dashes
    enough to tell it, or clear of a curve crossing; else over the span; ``intervals``
    is the fewest intervals between dashes that must tell it (:func:`style`). A file
    name starting with ``check:`` is read for :func:`checks` only and not written.
    """

    file: str
    figure: str
    style: str
    seed: tuple[int, int]
    span: tuple[float, float]
    checked: tuple[float, float] | None = None
    intervals: int = 2


def _file(name: str) -> str:
    """Return a dataset's data file."""
    return f"Ar/{SOURCE}__{name}_Ar.csv"


SEGMENTS: tuple[Segment, ...] = (
    Segment(
        _file("fig1b_all_L_p1S"),
        "fig1b",
        "solid",
        (3100, 4869),
        (8.1, 10.9),
        (8.1, 11.1),
    ),
    Segment(_file("fig1b_all_V_p1S"), "fig1b", "dashed", (3100, 4921), (8.1, 11.0)),
    Segment(
        _file("fig1b_noshakeup_L_p1S"), "fig1b", "dash-dot", (3100, 4724), (8.1, 11.6)
    ),
    Segment(
        _file("fig1b_noshakeup_V_p1S"),
        "fig1b",
        "dash-double-dot",
        (3100, 5222),
        (8.1, 9.3),
    ),
    Segment(_file("fig2_all_L_p1D"), "fig2", "solid", (900, 1621), (8.1, 9.7)),
    Segment(_file("fig3_all_L_f1D"), "fig3", "solid", (972, 5004), (8.1, 11.8)),
    Segment(_file("fig3_all_L_f1D"), "fig3", "solid", (1746, 5117), (12.3, 14.0)),
    Segment(_file("fig4a_hf_L"), "fig4a", "solid", (3011, 1559), (8.1, 11.5)),
    Segment(_file("fig4a_hf_L"), "fig4a", "solid", (3996, 1769), (12.4, 14.0)),
    Segment(_file("fig4a_hf_V"), "fig4a", "dashed", (3013, 1938), (8.1, 11.6)),
    Segment(_file("fig4a_hf_V"), "fig4a", "dashed", (3997, 2016), (12.4, 14.0)),
    Segment(_file("fig4b_noshakeup_L"), "fig4b", "solid", (3034, 3191), (8.1, 11.5)),
    Segment(_file("fig4b_noshakeup_L"), "fig4b", "solid", (4061, 3543), (12.9, 14.0)),
    Segment(_file("fig4b_noshakeup_V"), "fig4b", "dashed", (3034, 3453), (8.1, 11.5)),
    Segment(_file("fig4b_noshakeup_V"), "fig4b", "dashed", (4061, 3478), (12.9, 14.0)),
    Segment(_file("fig4c_all_L"), "fig4c", "solid", (3051, 4627), (8.1, 9.3)),
    Segment(_file("fig4c_all_L"), "fig4c", "solid", (3488, 4737), (9.9, 11.4)),
    Segment(_file("fig4c_all_L"), "fig4c", "solid", (4035, 4895), (12.7, 14.0)),
    Segment(_file("fig4c_all_V"), "fig4c", "dashed", (3000, 4675), (8.1, 9.3)),
    Segment(_file("fig4c_all_V"), "fig4c", "dashed", (3488, 4715), (9.9, 11.1)),
    Segment(_file("fig4c_all_V"), "fig4c", "dashed", (4035, 4841), (12.7, 14.0)),
    Segment(_file("fig5_hf_L_circular"), "fig5", "dash-dot", (817, 1385), (8.1, 11.5)),
    Segment(
        _file("fig5_hf_L_circular"), "fig5", "dash-dot", (1915, 1337), (12.4, 14.0)
    ),
    Segment(
        _file("fig5_hf_V_circular"),
        "fig5",
        "dash-double-dot",
        (1169, 1468),
        (9.6, 11.0),
        (9.6, 10.4),
    ),
    Segment(
        _file("fig5_hf_V_circular"),
        "fig5",
        "dash-double-dot",
        (1959, 1488),
        (13.5, 14.0),
        (13.5, 13.9),
        1,
    ),
    Segment(_file("fig5_all_L_circular"), "fig5", "solid", (817, 1467), (8.1, 10.9)),
    Segment(_file("fig5_all_L_circular"), "fig5", "solid", (1915, 1581), (12.4, 14.0)),
    Segment(_file("fig5_all_V_circular"), "fig5", "dashed", (1300, 1437), (9.2, 11.3)),
    Segment(
        _file("fig5_all_V_circular"),
        "fig5",
        "dashed",
        (1959, 1508),
        (13.5, 14.0),
        (13.0, 14.0),
    ),
    Segment("check:fig6a_all_L", "fig6a", "solid", (3314, 3404), (8.2, 11.4)),
    Segment("check:fig6a_all_L", "fig6a", "solid", (4085, 3631), (12.7, 13.8)),
)

# Samples (eV) that cannot be read cleanly (another curve within CLEAR pixels, or the
# line too steep), by file.
LEFT_OUT: dict[str, tuple[float, ...]] = {
    _file("fig1b_all_V_p1S"): (8.1,),
    _file("fig1b_noshakeup_L_p1S"): (8.1, 8.2, 10.5, 10.6, 11.4, 11.5, 11.6),
    _file("fig2_all_L_p1D"): (8.2, 8.3, 9.7),
    _file("fig3_all_L_f1D"): (8.1, 11.8, 12.3, 14.0),
    _file("fig4a_hf_L"): (11.3, 11.4, 11.5, 12.4, 13.9, 14.0),
    _file("fig4a_hf_V"): (14.0,),
    _file("fig4b_noshakeup_L"): (11.4, 11.5, 14.0),
    _file("fig4b_noshakeup_V"): (13.9, 14.0),
    _file("fig4c_all_L"): (8.1, 8.2, 8.3, 9.9, 10.0, 11.2, 11.3, 11.4, 12.7, 14.0),
    _file("fig4c_all_V"): (9.9, 10.0, 13.9, 14.0),
    _file("fig5_hf_L_circular"): (8.1, 8.2, 12.4, 12.5, 12.6, 12.7, 12.8, 14.0),
    _file("fig5_hf_V_circular"): (10.1, 10.2, 10.3, 10.4, 10.9, 11.0, 13.5),
    _file("fig5_all_L_circular"): (
        8.6,
        8.7,
        8.8,
        8.9,
        9.3,
        9.4,
        9.5,
        9.6,
        10.1,
        10.2,
        10.3,
        10.4,
        10.9,
        12.4,
        12.5,
        12.6,
        12.7,
        14.0,
    ),
    _file("fig5_all_V_circular"): (9.3, 9.4, 9.5, 9.6, 11.3),
    "check:fig6a_all_L": (8.2, 8.3, 12.7),
}
# Samples read on two curves of one figure within a line width of each other
# (:func:`_together`): left out of both.
TOGETHER: dict[str, tuple[float, ...]] = {
    _file("fig5_all_L_circular"): (9.2, 10.5, 10.6, 10.7, 10.8),
    _file("fig5_all_V_circular"): (9.2,),
    _file("fig5_hf_V_circular"): (10.5, 10.6, 10.7, 10.8),
}

# The tracing (:func:`trace`), pixels: the centres a prediction is fitted through, the
# fewest for a fit, how far a run may lie from the prediction, the widest run that
# moves the trace.
FIT = 30
MIN_FIT = 6
STEP = 6
MAX_WIDTH = 12
# The line style (:func:`style`), columns: the longest dot, the longest dash of a
# dashed line, the shortest solid stretch, the most columns a solid line may miss,
# the shortest gap that hides a dot, the fewest intervals between dashes that tell a
# style and the share of them that must agree.
DOT = 10
DASH = 80
MIN_SOLID = 160
SPECK = 3
LONG_GAP = 14
MIN_INTERVALS = 2
MAJORITY = 0.7
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


def _centre(run: Run) -> float:
    """Return a run's centre."""
    return (run[0] + run[1]) / 2


def _frame_row(dark: np.ndarray, column: int, row: int) -> float:
    """Return the centre row of a frame line crossing a column near a row."""
    low = row - FRAME_SEARCH
    found = [
        (abs(_centre(run) + low - row), _centre(run) + low)
        for run in _spans(dark[low : row + FRAME_SEARCH + 1, column])
        if run[1] - run[0] + 1 <= 2 * FRAME_SEARCH
    ]
    if not found:
        msg = f"no frame line near row {row} at column {column}"
        raise ValueError(msg)
    return min(found)[1]


class Calibration:
    """Map a pixel (column, row) to (photon energy, value) through a figure's ticks.

    The photon energy is interpolated between the bottom and top ticks' scales by the
    row, the value between the left and right ticks' scales by the column: the scan's
    small rotation and stretch. The frame lines at labelled corners are y ticks.
    """

    def __init__(self, figure: Figure, dark: np.ndarray) -> None:
        self.figure = figure
        y_sides = {"left": list(figure.y_left), "right": list(figure.y_right)}
        for side, column in (
            ("left", figure.left + DEPTH),
            ("right", figure.right - DEPTH),
        ):
            for row, value in zip(
                (figure.bottom, figure.top), figure.corners, strict=True
            ):
                if value is not None:
                    y_sides[side].append((_frame_row(dark, column, row), value))
        self.ticks = {
            "bottom": figure.x_bottom,
            "top": figure.x_top,
            "left": tuple(sorted(y_sides["left"])),
            "right": tuple(sorted(y_sides["right"])),
        }
        axes = {
            name: Axis(
                tuple(p for p, _ in ticks),
                tuple(float(v) for _, v in ticks),
                tolerance=TICK_TOLERANCE,
            )
            for name, ticks in self.ticks.items()
        }
        for axis in axes.values():
            axis.check()
        self.axes = axes
        # The scales are linear: value = slope * pixel + intercept.
        self.scales = {
            name: (float(axis(1.0)) - float(axis(0.0)), float(axis(0.0)))
            for name, axis in axes.items()
        }
        self.rows = (figure.top + DEPTH, figure.bottom - DEPTH)
        self.columns = (figure.left + DEPTH, figure.right - DEPTH)

    def residual(self) -> float:
        """Return the worst tick residual of the four sides, pixels."""
        worst = 0.0
        for name, (slope, intercept) in self.scales.items():
            for position, value in self.ticks[name]:
                worst = max(worst, abs((value - intercept) / slope - position))
        return worst

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
    """Raise unless every listed tick, and the frame, are found in the image."""
    rows = slice(figure.top + 5, figure.bottom - 5)
    columns = slice(figure.left + 5, figure.right - 5)
    probes = (
        (
            "bottom",
            lambda d: dark[figure.bottom - d, columns],
            figure.left,
            figure.x_bottom,
        ),
        ("top", lambda d: dark[figure.top + d, columns], figure.left, figure.x_top),
        ("left", lambda d: dark[rows, figure.left + d], figure.top, figure.y_left),
        ("right", lambda d: dark[rows, figure.right - d], figure.top, figure.y_right),
    )
    for name, probe, offset, ticks in probes:
        found = [
            offset + 5 + (a + b) / 2
            for depth in (DEPTH, DEEP)
            for a, b in _spans(probe(depth))
        ]
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


def style(found: list[bool], intervals: int = MIN_INTERVALS) -> str:
    """Return a traced line's style from which of its columns have a run.

    Before its first and after its last run the line is not drawn. ``solid``: gaps
    of at most ``SPECK`` columns in all (the scan drops a pixel column here and
    there), over at least ``MIN_SOLID`` columns (longer than any dash). Otherwise no
    run may be longer than ``DASH`` (a longer one is another, solid line), and between
    consecutive dashes (runs longer than ``DOT``) the dots (shorter runs) are counted
    where every gap is shorter than ``LONG_GAP`` (a longer gap is a dot the trace
    missed, off the line by a pixel or two): the count of at least ``MAJORITY`` of
    these intervals, and of ``intervals`` (``MIN_INTERVALS``) or more, names the
    style, 0 ``dashed``, 1 ``dash-dot``, 2 ``dash-double-dot``. Anything else is
    ``unknown``.
    """
    pieces = _pieces(found)
    while pieces and not pieces[0][0]:
        pieces.pop(0)
    while pieces and not pieces[-1][0]:
        pieces.pop()
    gaps = sum(n for flag, n in pieces if not flag)
    if gaps <= SPECK:
        return "solid" if sum(n for _, n in pieces) >= MIN_SOLID else "unknown"
    if max(n for flag, n in pieces if flag) > DASH:
        return "unknown"
    dashes = [i for i, (flag, n) in enumerate(pieces) if flag and n > DOT]
    counts = [
        sum(1 for flag, n in pieces[a + 1 : b] if flag)
        for a, b in itertools.pairwise(dashes)
        if all(flag or n < LONG_GAP for flag, n in pieces[a + 1 : b])
    ]
    if len(counts) < intervals:
        return "unknown"
    count = max(set(counts), key=counts.count)
    names = {0: "dashed", 1: "dash-dot", 2: "dash-double-dot"}
    if counts.count(count) < MAJORITY * len(counts) or count not in names:
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
        self.open = {c: run is None or clear[c] for c, run in line.items()}
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
    dark: np.ndarray, calibration: Calibration, segment: Segment, strict: bool = True
) -> dict[int, Run | None]:
    """Return a segment's traced line; raise if its style is not the one written."""
    low, high = segment.checked or segment.span
    row = segment.seed[1]
    first = int(calibration.column(min(low, segment.span[0]), row))
    last = int(np.ceil(calibration.column(max(high, segment.span[1]), row)))
    line = trace(dark, segment.seed, (first, last))
    found = [
        line[c] is not None
        for c in range(
            int(calibration.column(low, row)),
            int(np.ceil(calibration.column(high, row))) + 1,
        )
    ]
    drawn = style(found, segment.intervals)
    if drawn != segment.style:
        msg = f"{segment.file} {segment.span}: a {drawn} line, not {segment.style}"
        if strict:
            raise ValueError(msg)
        print(msg, _pieces(found)[:40])
    return line


def calibrations() -> dict[str, Calibration]:
    """Return every figure's calibration, after checking its ticks in the image."""
    pages: dict[int, np.ndarray] = {}
    out = {}
    for name, figure in FIGURES.items():
        if figure.page not in pages:
            pages[figure.page] = image(DIRECTORY, figure.page) < DARK
        _check_ticks(pages[figure.page], figure)
        out[name] = Calibration(figure, pages[figure.page])
    return out


def _together(out: dict[str, list[Reading]]) -> dict[str, list[float]]:
    """Return the samples where two curves of one figure are read within ``LINE`` rows.

    There the two lines run together and the reading cannot say whose it is.
    """
    figure_of = {segment.file: segment.figure for segment in SEGMENTS}
    found: dict[str, set[float]] = {}
    for a, b in itertools.combinations(sorted(out), 2):
        if figure_of[a] != figure_of[b]:
            continue
        rows_b = {r.wanted: r.row for r in out[b]}
        for r in out[a]:
            if r.wanted in rows_b and abs(r.row - rows_b[r.wanted]) < LINE:
                found.setdefault(a, set()).add(r.wanted)
                found.setdefault(b, set()).add(r.wanted)
    return {name: sorted(samples) for name, samples in found.items()}


def read_all(strict: bool = True) -> dict[str, list[Reading]]:
    """Return every segment's readings, by data file, in photon-energy order.

    Raises if a segment's line style is not the one written, or (if ``strict``) the
    samples left out are not exactly those of :data:`LEFT_OUT` and :data:`TOGETHER`.
    """
    pages: dict[int, np.ndarray] = {}
    scales = calibrations()
    out: dict[str, list[Reading]] = {}
    missing: dict[str, list[float]] = {}
    for segment in SEGMENTS:
        figure = FIGURES[segment.figure]
        if figure.page not in pages:
            pages[figure.page] = image(DIRECTORY, figure.page) < DARK
        dark = pages[figure.page]
        calibration = scales[segment.figure]
        line = _Line(dark, figure, _follow(dark, calibration, segment, strict))
        for wanted in _grid(segment.span):
            reading = _read(calibration, line, wanted)
            if reading is None:
                missing.setdefault(segment.file, []).append(wanted)
            else:
                out.setdefault(segment.file, []).append(reading)
    together = _together(out)
    for name, samples in together.items():
        out[name] = [r for r in out[name] if r.wanted not in samples]
    for name in set(missing) | set(LEFT_OUT) | set(together) | set(TOGETHER):
        if strict and tuple(missing.get(name, ())) != LEFT_OUT.get(name, ()):
            msg = f"{name}: left out {missing.get(name)}, not {LEFT_OUT.get(name)}"
            raise ValueError(msg)
        if strict and tuple(together.get(name, ())) != TOGETHER.get(name, ()):
            msg = f"{name}: together {together.get(name)}, not {TOGETHER.get(name)}"
            raise ValueError(msg)
    if not strict:
        print("LEFT_OUT", {k: tuple(v) for k, v in missing.items()})
        print("TOGETHER", {k: tuple(v) for k, v in together.items()})
    return {name: sorted(rows, key=lambda r: r.wanted) for name, rows in out.items()}


# Table II (p. 9/6279): f_q of Eq. 53.
TABLE_II = "tbl_524fc7b9d556"
# Eq. 54 (p. 9/6279): 'L^4 T (a.u.) = 1.8967 x 10^-50 cm^4 sec'; a0 (cm) and the
# atomic unit of time (s), CODATA 2018.
EQ_54 = 1.8967e-50
BOHR_CM = 5.29177210903e-9
TIME_AU_S = 2.4188843265857e-17
# How far the totals may stray from the sums of the channels, and Fig. 6(a) from
# Fig. 4(c).
SUMS = 0.04
SAME = 0.02
# The comparisons: the photon energies (eV) asked for, each taken at the nearest sample
# read, and the region below the 4s resonance's rise.
COMPARE_AT = (8.0, 9.0, 10.0, 11.0)
BELOW_4S = 11.15
OTHERS = {
    "mckenna_2004": "mckenna_2004/fig5_sigma2_Ar",
    "petrov_2016 corr+CP L": "petrov_2016/fig10_corrcp_L",
    "petrov_2016 corr+CP V": "petrov_2016/fig10_corrcp_V",
    "petrov_2016 corr L": "petrov_2016/fig10_corr_L",
    "starace_1987 ESI": "starace_1987/fig7_esi_linear_Ar",
}


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
    """Return each dataset's values, cm^4 s, by sample, keyed by its short name."""
    units = {segment.file: FIGURES[segment.figure].unit for segment in SEGMENTS}
    return {
        name.removeprefix("check:")
        .removeprefix(f"Ar/{SOURCE}__")
        .removesuffix("_Ar.csv"): {r.wanted: r.value * units[name] for r in rows}
        for name, rows in readings.items()
    }


def _ratios(
    top: dict[float, float], bottom: dict[float, float], below: float = math.inf
) -> tuple[float, float, float, float, int]:
    """Return the first and last shared samples below an energy, ratios and count."""
    shared = sorted(e for e in top if e in bottom and e < below)
    ratios = [top[e] / bottom[e] for e in shared]
    return shared[0], shared[-1], min(ratios), max(ratios), len(shared)


def _sum(parts: list[dict[float, float]], factor: float = 1.0) -> dict[float, float]:
    """Return the sum of datasets times a factor where all are read."""
    shared = set(parts[0]).intersection(*parts[1:])
    return {e: factor * sum(p[e] for p in parts) for e in shared}


def _resolution() -> dict[str, tuple[float, float]]:
    """Return each figure's pixel in x (eV) and in y (cm^4 s)."""
    out = {}
    for name, calibration in calibrations().items():
        x = abs(calibration.scales["bottom"][0])
        y = abs(calibration.scales["left"][0]) * FIGURES[name].unit
        out[name] = (round(x, 4), float(f"{y:.2e}"))
    return out


def _dataset(key: str) -> Dataset:
    """Return another record's dataset by its key, ``<source>/<id>``."""
    source, _ = key.split("/")
    record = load_source(DATA_DIRECTORY / "sources" / f"{source}.toml")
    return next(d for d in record.datasets if d.key == key)


def comparisons(values: dict[str, dict[float, float]]) -> dict[str, object]:
    """Return Fig. 4(c)'s totals against the other calculations of the library.

    At the samples nearest :data:`COMPARE_AT` (this paper over the other), and the range
    of the ratio below the 4s resonance's rise (:data:`BELOW_4S`).
    """
    others = {label: _dataset(key) for label, key in OTHERS.items()}
    out: dict[str, object] = {}
    for gauge in ("L", "V"):
        ours = values[f"fig4c_all_{gauge}"]
        for label, other in others.items():
            energy, sigma = other.abscissa_values(), other.scaled_values()
            theirs = {
                e: float(np.interp(e, energy, sigma))
                for e in ours
                if energy.min() <= e <= energy.max()
            }
            at = []
            for wanted in COMPARE_AT:
                e = min(theirs, key=lambda x: abs(x - wanted))
                at.append((e, round(ours[e] / theirs[e], 3)))
            out[f"fig4c {gauge} / {label}"] = (
                at,
                tuple(round(x, 3) for x in _ratios(ours, theirs, BELOW_4S)),
            )
    return out


def checks(readings: dict[str, list[Reading]]) -> dict[str, object]:
    """Return the checks of the reading; raise if one fails its margin.

    Table II's f_q (the circular 1D channels 3/2 the linear, no circular 1S); Eq. 54
    against CODATA; Fig. 4(c) against Figs. 1(b) + 2 + 3 (Eq. 49) and Fig. 5's
    all-process length curve against 3/2 times Figs. 2 + 3 where all are read
    (:data:`SUMS`); Fig. 6(a) against Fig. 4(c) (:data:`SAME`). Reported only: the
    pixel sizes; the shake-up's lowering of the p 1S channel and of the total (the
    text's 'approximately 0.1 x 10^-50', p. 12/6282); length over velocity; and
    :func:`comparisons`.
    """
    out: dict[str, object] = {}
    failures: list[str] = []
    f_q = {
        (q.strip(), term.strip()): _fraction(value)
        for q, term, value in table(DIRECTORY, TABLE_II)[1:]
    }
    ratio = f_q[("±1", "^{1}D")] / f_q[("0", "^{1}D")]
    out["Table II circular/linear 1D"] = ratio
    if abs(ratio - 1.5) > 1e-12 or f_q[("±1", "^{1}S")] != 0:
        failures.append("Table II: circular is not 3/2 linear for 1D")
    unit = BOHR_CM**4 * TIME_AU_S
    out["Eq. 54 against CODATA"] = (EQ_54, float(f"{unit:.5e}"))
    if abs(EQ_54 / unit - 1) > 1e-4:
        failures.append(f"Eq. 54: {EQ_54} against {unit:.5e}")
    out["pixel (eV, cm4 s)"] = _resolution()
    values = _values(readings)
    p1s, p1d, f1d = (
        values[name] for name in ("fig1b_all_L_p1S", "fig2_all_L_p1D", "fig3_all_L_f1D")
    )
    for name, parts, factor in (
        ("fig4c_all_L", [p1s, p1d, f1d], 1.0),
        ("fig5_all_L_circular", [p1d, f1d], 1.5),
    ):
        result = _ratios(values[name], _sum(parts, factor))
        out[f"{name} / the channels' sum"] = tuple(round(x, 3) for x in result)
        if max(abs(result[2] - 1), abs(result[3] - 1)) > SUMS:
            failures.append(
                f"{name}: {result[2]:.3f} to {result[3]:.3f} of the channels"
            )
    result = _ratios(values["fig6a_all_L"], values["fig4c_all_L"])
    out["fig6a / fig4c (L)"] = tuple(round(x, 3) for x in result)
    if max(abs(result[2] - 1), abs(result[3] - 1)) > SAME:
        failures.append(f"fig6a: {result[2]:.3f} to {result[3]:.3f} of fig4c")
    for gauge in ("L", "V"):
        no, all_ = (
            values[f"fig1b_noshakeup_{gauge}_p1S"],
            values[f"fig1b_all_{gauge}_p1S"],
        )
        shared = sorted(set(no) & set(all_))
        drop = [no[e] - all_[e] for e in shared]
        out[f"p 1S {gauge}: without - with shake-up, cm4 s"] = (
            shared[0],
            shared[-1],
            float(f"{min(drop):.2e}"),
            float(f"{max(drop):.2e}"),
        )
        no, all_ = values[f"fig4b_noshakeup_{gauge}"], values[f"fig4c_all_{gauge}"]
        shared = sorted(e for e in set(no) & set(all_) if e < BELOW_4S)
        drop = [no[e] - all_[e] for e in shared]
        out[f"total {gauge}: Fig. 4(b) - Fig. 4(c) below 4s, cm4 s"] = (
            shared[0],
            shared[-1],
            float(f"{min(drop):.2e}"),
            float(f"{max(drop):.2e}"),
        )
    for name in ("fig4a_hf", "fig4b_noshakeup", "fig4c_all"):
        result = _ratios(values[f"{name}_L"], values[f"{name}_V"], BELOW_4S)
        out[f"{name} L/V below 4s"] = tuple(round(x, 3) for x in result)
    out.update(comparisons(values))
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
        if not name.startswith("check:")
    }


def extract() -> dict[str, list[Point]]:
    """Return the cross sections read from Figs. 1(b)-5, cm^4 s against photon eV."""
    readings = read_all()
    checks(readings)
    return _points(readings)


def main() -> None:
    """Write the data files and print the checks."""
    readings = read_all()
    results = checks(readings)
    write(_points(readings))
    for label, result in results.items():
        print(f"{label}: {result}")


if __name__ == "__main__":
    main()

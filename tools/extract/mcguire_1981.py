# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract McGuire (1981): two-photon ionisation cross sections of Ne, Ar, Kr and Xe.

The paper prints no cross section as a number: they are in figures only (data rule 8,
TC-699). **Figs. 1-4** (pp. 842-843, library figures ``fig_ce2be8a0b31b``,
``fig_0ad46010a4ff``, ``fig_79c7e7519da3``, ``fig_b8ffcec227ce``) draw, for Ne, Ar, Kr
and Xe, the two-photon cross section for circularly polarised light sigma_2^C (solid
line, log left axis, cm^4 s) and the ratio R = sigma_2^C / sigma_2^L (crosses joined by
a line, log right axis) against the wavelength in angstrom. The PDF is a scan: each page
is one 1-bit image (4816 x 6367 pixels, about 600 dpi), read here at its native
resolution:

* **the frame**: each side's line found every 20 pixels near the position written in
  ``Plot.frame`` (a dark run 4 to 10 pixels thick) and fitted with a straight line
  (``_frame_line``; the scan's sides tilt by up to 0.005, a few pixels across a plot).
  A pixel is referred to the frame: its row to the left axis's row at the same height
  in the frame (``_left_row``), its column to the bottom axis along the frame's sides
  (``_bottom_column``).
* **the wavelength axis is piecewise linear**: the author expanded it at the
  resonances and compressed it elsewhere (p. 841). It is calibrated segment by
  segment (``Plot.x_segments``), each a run of ticks at one step: the labelled ticks
  and, in the compressed stretches, the unlabelled minor ticks between two labels,
  which divide the interval evenly (values by that division). The frame's left and
  right lines are ticks where a label stands at the corner. Each segment is a straight
  line through its ticks, checked (``Axis.check``, ``X_TOLERANCE``); a value is read
  only inside a segment. Stretches spanned by one division between two labels (two
  ticks, so no check) are not segments, and are left out (comments in ``PLOTS``).
* **the cross-section axis** (left) and **the ratio axis** (right) are logarithmic,
  ticked at every decade and at 2-9 times it; each decade is calibrated by itself from
  its ticks (``Y_TOLERANCE``), as the decades differ in length by up to 3 % in one
  figure. Ticks hidden by a curve touching the frame are not used (``Plot``).
* every tick written in ``Plot`` is found again in the image by the script (a dark run
  4-14 pixels wide crossing a probe line 7 pixels inside the frame and at least 10
  long) within ``TICK_FIND`` pixels.
* **the cross section** is traced column by column from seed pixels on the solid line
  (``Plot.sigma_seeds``), following the dark run nearest to the previous position
  within ``STEP`` pixels, below the R = 1.5 line (``RATIO_MARGIN``). It is sampled at
  regular wavelengths (``Plot.steps``, per segment); at each, the traced column whose
  wavelength is nearest, and the centre of the line's run there. A sample is taken
  only where that run is at most ``MAX_RUN`` pixels tall (the line is 8-12 thick): the
  steep flanks and the peaks of the resonances are left out, as are the stretches
  listed in ``Plot.sigma_left_out``.
* **the ratio** is traced the same way (``Plot.ratio_seeds``) in the band below the
  R = 1.5 line, whose position is fitted like a frame line; the crosses and the joining
  line are one curve, the dark pixels of a column within ``JOIN`` pixels of each other
  one cluster, whose centre is read (a cross is symmetric about its centre). A sample
  is taken where the cluster is at most ``MAX_CLUSTER`` pixels tall and does not touch
  the R = 1.5 line (where the curve runs into the line it cannot be told from it).

Each value's ``printed`` cell gives its pixel (column, row) in the page image. One
pixel is 0.7-0.9 % in sigma and R (265-335 pixels per decade); the line is 8-12 pixels
wide. ``checks`` compares the reading with the numbers printed for these curves: the
step of Fig. 4 'at 1843 A' (p. 844), the drawn 'R = 1.5' lines, and McCown et al.'s
quotation of the Xe cross section at 193 nm (mccown_1982); ``extract`` fails if one
is missed (the record's notes give the results).
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, write

DIRECTORY = "McGuire_1981_TwoThreePhoton"
DARK = 128
# How far, pixels, a written tick may be from the one found in the image.
TICK_FIND = 1.5
# The largest residual, pixels, of a tick on its segment's straight line (wavelength)
# or decade (log axes).
X_TOLERANCE = 4.0
Y_TOLERANCE = 4.5
# The probe line for ticks, pixels inside the frame line's centre, and the shortest
# tick.
PROBE = 7
MIN_TICK = 10
# Tracing: the largest jump, pixels, from one column's position to the next; the
# tallest run of the cross-section line read as a value; the gap, pixels, that joins
# the dark runs of the ratio's crosses and line into one cluster; the tallest cluster.
STEP = 12
MAX_RUN = 18
NEIGHBOURS = 3
JOIN = 4
MAX_CLUSTER = 45
# A sample is read at the traced column whose wavelength is within this many columns
# width of the one wanted.
HALF_COLUMN = 0.6
# Dark runs shorter than this, pixels, are specks of the scan, not ink of a line.
SPECK = 3
# The most columns a trace skips (no run near the last one) before it stops.
MAX_GAP = 15
# The cross section is traced below the R = 1.5 line by this many pixels.
RATIO_MARGIN = 14
# The ratio band: from this many pixels below the R = 1.5 line's centre (its edge)
# down ``RATIO_BAND`` pixels.
RATIO_TOP = 6
RATIO_BAND = 420


@dataclass(frozen=True)
class Plot:
    """One figure: where it is, its ticks, its seeds and what is left out."""

    figure: str
    gas: str
    page: int
    frame: tuple[float, float, float, float]
    """Approximate centres of the frame's top and bottom rows, left and right
    columns."""
    ratio_line: float
    """Approximate centre row of the drawn R = 1.5 line."""
    x_segments: tuple[tuple[tuple[float, float], ...], ...]
    """The wavelength axis: per segment, the bottom ticks (column, angstrom)."""
    left_ticks: tuple[tuple[float, float], ...]
    """The cross-section axis: ticks (row, cm^4 s) on the left side."""
    right_ticks: tuple[tuple[float, float], ...]
    """The ratio axis: ticks (row, R) on the right side."""
    steps: tuple[float, ...]
    """The sampling step of each wavelength segment, angstrom."""
    sigma_seeds: tuple[tuple[int, int], ...]
    ratio_seeds: tuple[tuple[int, int], ...]
    sigma_left_out: tuple[tuple[float, float], ...] = ()
    """Wavelength ranges (angstrom) of the cross section not read (see the record)."""
    ratio_left_out: tuple[tuple[float, float], ...] = ()


PLOTS: dict[str, Plot] = {
    "fig1": Plot(
        figure="Fig. 1",
        gas="Ne",
        page=8,
        frame=(670, 2521, 685, 4189),
        ratio_line=880,
        # The frame's left line is 592 (the label at the corner), its right line 1150.
        x_segments=(
            (
                (683.0, 592),
                (793.0, 593),
                (899.0, 594),
                (1003.5, 595),
                (1112.5, 596),
                (1218.5, 597),
                (1327.5, 598),
                (1435.5, 599),
                (1541.0, 600),
                (1648.5, 601),
                (1754.5, 602),
                (1862.5, 603),
                (1968.5, 604),
            ),
            # 604 to 610: five minor ticks, 1 angstrom apart.
            (
                (1968.5, 604),
                (1985.5, 605),
                (2007.0, 606),
                (2029.0, 607),
                (2049.5, 608),
                (2071.0, 609),
                (2093.0, 610),
            ),
            (
                (2093.0, 610),
                (2197.0, 615),
                (2305.0, 620),
                (2412.0, 625),
                (2516.0, 630),
                (2623.0, 635),
                (2728.0, 640),
            ),
            # 640 to 700: five minor ticks, 10 angstrom apart.
            (
                (2728.0, 640),
                (2748.5, 650),
                (2769.5, 660),
                (2789.0, 670),
                (2812.5, 680),
                (2833.5, 690),
                (2854.0, 700),
            ),
            # 700 to 730: two minor ticks, 10 angstrom apart.
            ((2854.0, 700), (2874.5, 710), (2894.0, 720), (2919.5, 730)),
            (
                (2919.5, 730),
                (3020.0, 735),
                (3124.5, 740),
                (3231.0, 745),
                (3335.0, 750),
            ),
            (
                (3335.0, 750),
                (3441.5, 800),
                (3551.5, 850),
                (3658.5, 900),
                (3762.0, 950),
                (3871.0, 1000),
                (3975.0, 1050),
                (4083.5, 1100),
                (4188.5, 1150),
            ),
        ),
        # The ticks at 4 and 3e-50 are hidden where the curve meets the frame.
        left_ticks=(
            (682.0, 9e-46),
            (695.0, 8e-46),
            (709.5, 7e-46),
            (727.5, 6e-46),
            (749.5, 5e-46),
            (773.0, 4e-46),
            (805.5, 3e-46),
            (855.5, 2e-46),
            (931.5, 1e-46),
            (944.5, 9e-47),
            (960.0, 8e-47),
            (974.0, 7e-47),
            (992.0, 6e-47),
            (1012.0, 5e-47),
            (1038.0, 4e-47),
            (1070.0, 3e-47),
            (1117.0, 2e-47),
            (1197.0, 1e-47),
            (1212.5, 9e-48),
            (1225.0, 8e-48),
            (1241.0, 7e-48),
            (1258.5, 6e-48),
            (1277.5, 5e-48),
            (1304.0, 4e-48),
            (1337.0, 3e-48),
            (1384.0, 2e-48),
            (1462.0, 1e-48),
            (1477.5, 9e-49),
            (1491.5, 8e-49),
            (1506.0, 7e-49),
            (1524.0, 6e-49),
            (1543.5, 5e-49),
            (1567.0, 4e-49),
            (1602.5, 3e-49),
            (1651.0, 2e-49),
            (1727.5, 1e-49),
            (1742.0, 9e-50),
            (1756.0, 8e-50),
            (1771.5, 7e-50),
            (1789.0, 6e-50),
            (1810.0, 5e-50),
            (1914.5, 2e-50),
            (1994.0, 1e-50),
            (2009.0, 9e-51),
            (2022.5, 8e-51),
            (2039.0, 7e-51),
            (2056.0, 6e-51),
            (2075.5, 5e-51),
            (2101.0, 4e-51),
            (2135.0, 3e-51),
            (2182.0, 2e-51),
            (2258.5, 1e-51),
            (2273.0, 9e-52),
            (2287.5, 8e-52),
            (2303.0, 7e-52),
            (2321.0, 6e-52),
            (2340.0, 5e-52),
            (2362.0, 4e-52),
            (2398.0, 3e-52),
            (2444.0, 2e-52),
        ),
        right_ticks=(
            (682.0, 9e0),
            (696.0, 8e0),
            (711.5, 7e0),
            (727.0, 6e0),
            (746.0, 5e0),
            (773.0, 4e0),
            (806.0, 3e0),
            (852.5, 2e0),
            (931.0, 1e0),
            (944.0, 9e-1),
            (957.0, 8e-1),
            (971.5, 7e-1),
            (989.5, 6e-1),
            (1012.0, 5e-1),
            (1036.0, 4e-1),
            (1067.5, 3e-1),
            (1115.5, 2e-1),
            (1194.5, 1e-1),
            (1207.5, 9e-2),
            (1221.5, 8e-2),
            (1237.5, 7e-2),
            (1255.0, 6e-2),
            (1275.0, 5e-2),
            (1301.0, 4e-2),
            (1334.0, 3e-2),
            (1380.0, 2e-2),
            (1462.5, 1e-2),
            (1475.5, 9e-3),
            (1489.0, 8e-3),
            (1503.0, 7e-3),
            (1520.5, 6e-3),
            (1542.0, 5e-3),
            (1566.5, 4e-3),
            (1600.5, 3e-3),
            (1646.0, 2e-3),
            (1725.0, 1e-3),
        ),
        steps=(0.2, 1, 1, 10, 10, 1, 10),
        sigma_seeds=(
            (1300, 2263),
            (2330, 2220),
            (2460, 2027),
            (2580, 2204),
            (3600, 2237),
        ),
        ratio_seeds=((884, 918), (1990, 918), (2580, 908)),
        # The dashed enlargement of the 3d'' resonance (0.02 angstrom per division)
        # crosses and meets the solid line from 611 to 617 angstrom, where the solid
        # line also breaks.
        sigma_left_out=((611.0, 617.0),),
    ),
    "fig2": Plot(
        figure="Fig. 2",
        gas="Ar",
        page=8,
        frame=(3505, 5550, 663, 4136),
        ratio_line=3738,
        # The frame's left line is 825 (the label at the corner), its right line 1400;
        # 910-950, 1000-1040 and 1070-1100 have minor ticks 10 angstrom apart.
        x_segments=(
            (
                (662.5, 825),
                (783.0, 830),
                (900.5, 835),
                (1020.0, 840),
                (1139.0, 845),
                (1257.5, 850),
                (1374.0, 855),
                (1493.0, 860),
                (1612.5, 865),
                (1731.0, 870),
                (1849.5, 875),
                (1969.5, 880),
                (2089.5, 885),
                (2206.0, 890),
                (2328.0, 895),
                (2446.0, 900),
                (2566.0, 905),
                (2683.5, 910),
            ),
            (
                (2683.5, 910),
                (2706.5, 920),
                (2728.0, 930),
                (2753.5, 940),
                (2775.0, 950),
            ),
            # 950 to 1000: one division, left out.
            (
                (2896.0, 1000),
                (2919.0, 1010),
                (2942.0, 1020),
                (2967.5, 1030),
                (2994.5, 1040),
            ),
            (
                (2994.5, 1040),
                (3111.5, 1050),
                (3232.5, 1060),
                (3349.0, 1070),
            ),
            (
                (3349.0, 1070),
                (3374.0, 1080),
                (3397.0, 1090),
                (3421.5, 1100),
            ),
            (
                (3421.5, 1100),
                (3539.0, 1150),
                (3656.5, 1200),
                (3774.5, 1250),
                (3892.0, 1300),
                (4011.5, 1350),
                (4132.0, 1400),
            ),
        ),
        # Ticks hidden where a curve meets the frame are not listed.
        left_ticks=(
            (3519.0, 9e-44),
            (3533.5, 8e-44),
            (3551.5, 7e-44),
            (3572.0, 6e-44),
            (3594.5, 5e-44),
            (3621.5, 4e-44),
            (3660.0, 3e-44),
            (3712.0, 2e-44),
            (3799.5, 1e-44),
            (3815.0, 9e-45),
            (3830.0, 8e-45),
            (3847.0, 7e-45),
            (3867.5, 6e-45),
            (3889.5, 5e-45),
            (3917.5, 4e-45),
            (3953.0, 3e-45),
            (4004.5, 2e-45),
            (4094.0, 1e-45),
            (4109.0, 9e-46),
            (4124.0, 8e-46),
            (4141.5, 7e-46),
            (4161.0, 6e-46),
            (4185.0, 5e-46),
            (4212.5, 4e-46),
            (4246.0, 3e-46),
            (4303.0, 2e-46),
            (4390.5, 1e-46),
            (4407.0, 9e-47),
            (4421.5, 8e-47),
            (4439.0, 7e-47),
            (4456.5, 6e-47),
            (4479.0, 5e-47),
            (4507.5, 4e-47),
            (4544.0, 3e-47),
            (4599.0, 2e-47),
            (4685.5, 1e-47),
            (4698.5, 9e-48),
            (4716.5, 8e-48),
            (4732.5, 7e-48),
            (4752.0, 6e-48),
            (4776.0, 5e-48),
            (4802.5, 4e-48),
            (4838.0, 3e-48),
            (4890.5, 2e-48),
            (4974.0, 1e-48),
            (4989.0, 9e-49),
            (5003.5, 8e-49),
            (5020.5, 7e-49),
            (5041.0, 6e-49),
            (5062.0, 5e-49),
            (5091.5, 4e-49),
            (5127.0, 3e-49),
            (5175.0, 2e-49),
            (5262.0, 1e-49),
            (5277.0, 9e-50),
            (5292.0, 8e-50),
            (5309.0, 7e-50),
            (5327.0, 6e-50),
            (5348.0, 5e-50),
            (5377.0, 4e-50),
            (5414.5, 3e-50),
            (5464.0, 2e-50),
        ),
        right_ticks=(
            (3518.0, 9e0),
            (3533.0, 8e0),
            (3550.5, 7e0),
            (3569.5, 6e0),
            (3591.5, 5e0),
            (3621.5, 4e0),
            (3659.0, 3e0),
            (3712.0, 2e0),
            (3815.0, 9e-1),
            (3828.5, 8e-1),
            (3845.5, 7e-1),
            (3867.5, 6e-1),
            (3888.5, 5e-1),
            (3917.0, 4e-1),
            (3954.0, 3e-1),
            (4005.5, 2e-1),
            (4093.5, 1e-1),
            (4108.0, 9e-2),
            (4123.5, 8e-2),
            (4140.5, 7e-2),
            (4159.0, 6e-2),
            (4185.0, 5e-2),
            (4211.0, 4e-2),
            (4249.0, 3e-2),
            (4301.0, 2e-2),
            (4387.5, 1e-2),
        ),
        steps=(1, 10, 10, 2, 10, 10),
        sigma_seeds=(
            (800, 5428),
            (1200, 5387),
            (1780, 5350),
            (1930, 5221),
            (2200, 5423),
            (2500, 5342),
            (3600, 5471),
            (3900, 5491),
        ),
        ratio_seeds=(
            (900, 3788),
            (1600, 3776),
            (2000, 3775),
            (2900, 3756),
            (3150, 3788),
            (3300, 3795),
        ),
    ),
    "fig3": Plot(
        figure="Fig. 3",
        gas="Kr",
        page=9,
        frame=(700, 2882, 681, 4144),
        ratio_line=946,
        # The frame's left line is 945 (the label at the corner; the right line is
        # unlabelled); 945-950 is ticked every angstrom, 1050-1100, 1170-1200 and
        # 1200-1230 have minor ticks 10 angstrom apart. Labels on a second row (1150,
        # 1200, 1250) stand under their own ticks.
        x_segments=(
            (
                (683.5, 945),
                (733.5, 946),
                (786.0, 947),
                (834.5, 948),
                (885.5, 949),
                (934.0, 950),
            ),
            (
                (934.0, 950),
                (1060.0, 960),
                (1186.0, 970),
                (1311.5, 980),
                (1438.5, 990),
                (1563.5, 1000),
                (1689.0, 1010),
                (1811.5, 1020),
                (1939.5, 1030),
                (2062.5, 1040),
                (2187.0, 1050),
            ),
            (
                (2187.0, 1050),
                (2210.5, 1060),
                (2236.0, 1070),
                (2262.5, 1080),
                (2289.0, 1090),
                (2311.5, 1100),
            ),
            # 1100 to 1150 and 1150 to 1160: one division each, left out.
            (
                (2460.0, 1160),
                (2581.5, 1165),
                (2709.5, 1170),
            ),
            (
                (2709.5, 1170),
                (2734.5, 1180),
                (2762.5, 1190),
                (2787.5, 1200),
            ),
            (
                (2787.5, 1200),
                (2812.5, 1210),
                (2836.5, 1220),
                (2859.5, 1230),
            ),
            (
                (2859.5, 1230),
                (2986.0, 1235),
                (3111.5, 1240),
            ),
            # 1240 to 1250: one division, left out.
            (
                (3139.0, 1250),
                (3265.5, 1300),
                (3390.5, 1350),
                (3513.5, 1400),
                (3638.0, 1450),
                (3766.0, 1500),
                (3893.5, 1550),
                (4017.5, 1600),
            ),
        ),
        # Ticks hidden where a curve meets the frame are not listed.
        left_ticks=(
            (714.5, 9e-44),
            (733.0, 8e-44),
            (749.0, 7e-44),
            (770.0, 6e-44),
            (793.5, 5e-44),
            (823.5, 4e-44),
            (864.0, 3e-44),
            (918.5, 2e-44),
            (1041.5, 8e-45),
            (1060.5, 7e-45),
            (1080.5, 6e-45),
            (1104.5, 5e-45),
            (1134.0, 4e-45),
            (1171.0, 3e-45),
            (1230.0, 2e-45),
            (1321.0, 1e-45),
            (1338.5, 9e-46),
            (1355.5, 8e-46),
            (1372.0, 7e-46),
            (1392.0, 6e-46),
            (1416.5, 5e-46),
            (1449.0, 4e-46),
            (1485.5, 3e-46),
            (1539.5, 2e-46),
            (1631.5, 1e-46),
            (1648.0, 9e-47),
            (1664.0, 8e-47),
            (1684.0, 7e-47),
            (1702.0, 6e-47),
            (1727.5, 5e-47),
            (1757.5, 4e-47),
            (1795.0, 3e-47),
            (1852.5, 2e-47),
            (1946.5, 1e-47),
            (1962.5, 9e-48),
            (1979.0, 8e-48),
            (1996.5, 7e-48),
            (2016.0, 6e-48),
            (2039.5, 5e-48),
            (2070.0, 4e-48),
            (2110.0, 3e-48),
            (2329.0, 6e-49),
            (2351.5, 5e-49),
            (2379.0, 4e-49),
            (2419.5, 3e-49),
            (2473.0, 2e-49),
            (2568.0, 1e-49),
            (2582.0, 9e-50),
            (2597.5, 8e-50),
            (2618.0, 7e-50),
            (2641.0, 6e-50),
            (2662.5, 5e-50),
            (2692.5, 4e-50),
            (2732.0, 3e-50),
            (2790.5, 2e-50),
        ),
        right_ticks=(
            (721.0, 9e0),
            (738.5, 8e0),
            (756.0, 7e0),
            (776.0, 6e0),
            (800.0, 5e0),
            (831.5, 4e0),
            (870.0, 3e0),
            (920.5, 2e0),
            (1020.5, 1e0),
            (1034.5, 9e-1),
            (1047.5, 8e-1),
            (1065.0, 7e-1),
            (1085.0, 6e-1),
            (1112.0, 5e-1),
            (1141.0, 4e-1),
            (1178.5, 3e-1),
            (1231.5, 2e-1),
            (1327.5, 1e-1),
            (1344.5, 9e-2),
            (1359.5, 8e-2),
            (1377.0, 7e-2),
            (1398.5, 6e-2),
            (1421.0, 5e-2),
            (1451.5, 4e-2),
            (1491.0, 3e-2),
            (1541.5, 2e-2),
            (1637.5, 1e-2),
        ),
        steps=(0.2, 2, 10, 1, 10, 10, 1, 10),
        sigma_seeds=(
            (1200, 2516),
            (1500, 2636),
            (2100, 2479),
            (2300, 2581),
            (2800, 2607),
            (3400, 2651),
            (3900, 2685),
        ),
        ratio_seeds=(
            (900, 970),
            (2100, 963),
            (2300, 962),
            (2450, 986),
            (2650, 995),
            (2800, 986),
            (3000, 1014),
            (3400, 976),
            (3900, 1002),
        ),
    ),
    "fig4": Plot(
        figure="Fig. 4",
        gas="Xe",
        page=9,
        frame=(3276, 5600, 1064, 3827),
        ratio_line=3553,
        # The frame's left line is 1100 (the label at the corner); 1300-1340 is labelled
        # every 10 angstrom on two rows at the 20-angstrom scale of 1100-1300 (one
        # segment), 1340-1400 has minor ticks 20 angstrom apart; 2000 is the last tick
        # before the frame's right line.
        x_segments=(
            (
                (1060.5, 1100),
                (1191.5, 1120),
                (1325.0, 1140),
                (1460.5, 1160),
                (1593.0, 1180),
                (1729.0, 1200),
                (1863.0, 1220),
                (1998.0, 1240),
                (2132.5, 1260),
                (2265.5, 1280),
                (2400.0, 1300),
                (2469.5, 1310),
                (2535.0, 1320),
                (2604.5, 1330),
                (2667.5, 1340),
            ),
            (
                (2667.5, 1340),
                (2695.0, 1360),
                (2723.5, 1380),
                (2751.5, 1400),
            ),
            # 1400 to 1460: one division, left out.
            (
                (2835.0, 1460),
                (2970.5, 1470),
                (3102.5, 1480),
            ),
            # 1480 to 1500: one division, left out.
            (
                (3133.5, 1500),
                (3265.0, 1600),
                (3399.0, 1700),
                (3532.0, 1800),
                (3667.5, 1900),
                (3801.0, 2000),
            ),
        ),
        # Ticks hidden where a curve meets the frame are not listed.
        left_ticks=(
            (3304.0, 8e-44),
            (3324.0, 7e-44),
            (3345.5, 6e-44),
            (3370.0, 5e-44),
            (3403.5, 4e-44),
            (3444.5, 3e-44),
            (3502.5, 2e-44),
            (3603.5, 1e-44),
            (3619.0, 9e-45),
            (3636.0, 8e-45),
            (3657.5, 7e-45),
            (3678.5, 6e-45),
            (3704.5, 5e-45),
            (3737.5, 4e-45),
            (3780.5, 3e-45),
            (3839.5, 2e-45),
            (3937.5, 1e-45),
            (3951.0, 9e-46),
            (3972.0, 8e-46),
            (3989.5, 7e-46),
            (4012.0, 6e-46),
            (4039.0, 5e-46),
            (4070.5, 4e-46),
            (4114.0, 3e-46),
            (4171.5, 2e-46),
            (4271.5, 1e-46),
            (4286.5, 9e-47),
            (4372.5, 5e-47),
            (4405.0, 4e-47),
            (4447.5, 3e-47),
            (4506.0, 2e-47),
            (4603.5, 1e-47),
            (4620.0, 9e-48),
            (4640.5, 8e-48),
            (4656.5, 7e-48),
            (4680.5, 6e-48),
            (4705.5, 5e-48),
            (4738.5, 4e-48),
            (4782.5, 3e-48),
            (4839.0, 2e-48),
            (4938.5, 1e-48),
            (4954.5, 9e-49),
            (4976.5, 8e-49),
            (4993.5, 7e-49),
            (5016.0, 6e-49),
            (5043.0, 5e-49),
            (5073.0, 4e-49),
            (5112.5, 3e-49),
            (5175.0, 2e-49),
            (5271.0, 1e-49),
            (5287.0, 9e-50),
            (5344.5, 6e-50),
            (5370.0, 5e-50),
            (5402.0, 4e-50),
            (5443.0, 3e-50),
            (5498.5, 2e-50),
        ),
        right_ticks=(
            (3297.5, 9e0),
            (3315.5, 8e0),
            (3334.5, 7e0),
            (3355.0, 6e0),
            (3381.5, 5e0),
            (3414.0, 4e0),
            (3456.0, 3e0),
            (3514.5, 2e0),
            (3614.5, 1e0),
            (3632.5, 9e-1),
            (3648.5, 8e-1),
            (3668.5, 7e-1),
            (3689.5, 6e-1),
            (3718.0, 5e-1),
            (3749.5, 4e-1),
            (3793.0, 3e-1),
            (3849.5, 2e-1),
            (3950.0, 1e-1),
            (3964.0, 9e-2),
            (3984.0, 8e-2),
            (4002.0, 7e-2),
            (4025.0, 6e-2),
            (4050.5, 5e-2),
            (4080.0, 4e-2),
            (4124.5, 3e-2),
            (4181.5, 2e-2),
            (4281.5, 1e-2),
        ),
        steps=(2, 20, 1, 10),
        sigma_seeds=(
            (1400, 5280),
            (1800, 5140),
            (2300, 5085),
            (2450, 5131),
            (2600, 5173),
            (2800, 5204),
            (3300, 5240),
            (3500, 5255),
        ),
        ratio_seeds=(
            (1400, 3566),
            (1550, 3588),
            (2300, 3567),
            (2800, 3568),
            (3000, 3606),
            (3300, 3570),
            (3500, 3583),
        ),
    ),
}

FILES = {
    (key, kind): f"{plot.gas}/mcguire_1981__{key}_{kind}_{plot.gas}.csv"
    for key, plot in PLOTS.items()
    for kind in ("sigma2", "ratio")
}


def _spans(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


@dataclass(frozen=True)
class Line:
    """A straight frame line: the row (horizontal) or column (vertical) a + b t."""

    a: float
    b: float
    worst: float
    """The largest residual of the fitted samples, pixels."""

    def __call__(self, t: float) -> float:
        """Return the line's row or column at column or row t."""
        return self.a + self.b * t


def _frame_line(
    dark: np.ndarray,
    approx: float,
    low: int,
    high: int,
    horizontal: bool,
    thickest: int = 10,
) -> Line:
    """Fit a frame line near ``approx`` between ``low`` and ``high`` (rows or columns).

    Every 20 pixels along it, the one dark run 4 to ``thickest`` pixels thick within 14
    pixels of
    ``approx``; a straight line through their centres, refitted without samples more
    than 2 pixels off (where a tick or a curve joins the line).
    """
    centre = round(approx)
    along, across = [], []
    for t in range(low, high, 20):
        if horizontal:
            strip = dark[centre - 14 : centre + 14, t]
        else:
            strip = dark[t, centre - 14 : centre + 14]
        runs = [(a, b) for a, b in _spans(strip) if 4 <= b - a + 1 <= thickest]
        if len(runs) == 1:
            along.append(t)
            across.append(centre - 14 + (runs[0][0] + runs[0][1]) / 2)
    t_values, values = np.array(along, dtype=float), np.array(across)
    keep = np.ones(len(values), dtype=bool)
    slope, intercept = 0.0, 0.0
    for _ in range(5):
        slope, intercept = np.polyfit(t_values[keep], values[keep], 1)
        keep = np.abs(values - (intercept + slope * t_values)) < 2.0
    if keep.sum() < 0.8 * len(values) or len(values) < 20:
        msg = f"a frame line near {approx}: {keep.sum()} of {len(values)} samples fit"
        raise ValueError(msg)
    residual = values[keep] - (intercept + slope * t_values[keep])
    return Line(float(intercept), float(slope), float(np.abs(residual).max()))


@dataclass(frozen=True)
class Frame:
    """The four sides of a plot's frame."""

    top: Line
    bottom: Line
    left: Line
    right: Line


def frame(dark: np.ndarray, plot: Plot) -> Frame:
    """Return a plot's frame, fitted to its four lines."""
    top, bottom, left, right = plot.frame
    first_row, last_row = round(top) + 30, round(bottom) - 30
    first_column, last_column = round(left) + 30, round(right) - 30
    return Frame(
        _frame_line(dark, top, first_column, last_column, horizontal=True),
        _frame_line(dark, bottom, first_column, last_column, horizontal=True),
        _frame_line(dark, left, first_row, last_row, horizontal=False),
        _frame_line(dark, right, first_row, last_row, horizontal=False),
    )


def _left_row(box: Frame, column: float, row: float) -> float:
    """Return the row on the left axis at the same height in the frame as a pixel."""
    top, bottom = box.top(column), box.bottom(column)
    edge = box.left(row)
    return box.top(edge) + (row - top) / (bottom - top) * (
        box.bottom(edge) - box.top(edge)
    )


def _bottom_column(box: Frame, column: float, row: float) -> float:
    """Return the bottom axis's column reached from a pixel along the frame's sides.

    The sides' slopes (column per row) are interpolated across the frame.
    """
    left, right = box.left(row), box.right(row)
    u = (column - left) / (right - left)
    slope = (1 - u) * box.left.b + u * box.right.b
    return column + slope * (box.bottom(column) - row)


def _ticks_on(
    dark: np.ndarray, line: Line, low: int, high: int, side: str
) -> list[float]:
    """Return the centres of the ticks on a frame line between ``low`` and ``high``.

    A tick is a dark run 4-14 pixels wide crossing the probe line ``PROBE`` pixels
    inside the frame line, and at least ``MIN_TICK`` pixels long from the frame line.
    """
    inward = {"bottom": -1, "top": 1, "left": 1, "right": -1}[side]
    horizontal = side in {"bottom", "top"}

    def pixel(t: int, depth: int) -> bool:
        across = round(line(t) + inward * depth)
        return bool(dark[across, t] if horizontal else dark[t, across])

    probe = np.array([pixel(t, PROBE) for t in range(low, high)])
    found = []
    for first, last in _spans(probe):
        if not 3 <= last - first <= 13:
            continue
        depth = PROBE
        while depth < 300 and any(
            pixel(t, depth) for t in range(low + first, low + last + 1)
        ):
            depth += 1
        if depth >= MIN_TICK:
            found.append(low + (first + last) / 2)
    return found


def _check_ticks(dark: np.ndarray, box: Frame, plot: Plot) -> None:
    """Raise unless every tick written in ``plot`` is found in the image."""
    top, bottom, left, right = (round(v) for v in plot.frame)
    row = box.bottom(box.left(plot.frame[1]))
    corners = (box.left(row), box.right(box.bottom(box.right(plot.frame[1]))))
    columns = sorted({c for segment in plot.x_segments for c, _ in segment})
    for column in columns:
        # A tick at a corner is the frame's side line, checked against its fit.
        if min(abs(column - corner) for corner in corners) < 3 * TICK_FIND and (
            min(abs(column - corner) for corner in corners) > TICK_FIND
        ):
            msg = f"{plot.figure}: the corner tick {column} is not on the frame"
            raise ValueError(msg)
    written = {
        "bottom": [
            c for c in columns if min(abs(c - corner) for corner in corners) > TICK_FIND
        ],
        "left": [r for r, _ in plot.left_ticks],
        "right": [r for r, _ in plot.right_ticks],
    }
    lines = {"bottom": box.bottom, "left": box.left, "right": box.right}
    for side, positions in written.items():
        if side == "bottom":
            found = _ticks_on(dark, lines[side], left + 10, right - 9, side)
        else:
            found = _ticks_on(dark, lines[side], top + 10, bottom - 9, side)
        for position in positions:
            if not found or min(abs(f - position) for f in found) > TICK_FIND:
                msg = f"{plot.figure}: no {side} tick at {position}"
                raise ValueError(msg)


@dataclass(frozen=True)
class Scale:
    """A piecewise axis: per segment, a straight line through its ticks."""

    segments: tuple[Axis, ...]
    ranges: tuple[tuple[float, float], ...]
    """Each segment's range of positions."""

    def __call__(self, position: float) -> float | None:
        """Return the value at a position, or None outside every segment."""
        for axis, (low, high) in zip(self.segments, self.ranges, strict=True):
            if low <= position <= high:
                return float(axis(position))
        return None

    def worst(self) -> list[float]:
        """Return each segment's largest tick residual, pixels."""
        out = []
        for axis in self.segments:
            y = np.log10(axis.values) if axis.log else np.asarray(axis.values)
            slope, intercept = np.polyfit(axis.positions, y, 1)
            out.append(float(np.max(np.abs((y - intercept) / slope - axis.positions))))
        return out


def x_scale(box: Frame, plot: Plot) -> Scale:
    """Return the wavelength axis of a plot, in bottom-axis columns."""
    segments, ranges = [], []
    for ticks in plot.x_segments:
        if len(ticks) < 3:
            msg = f"{plot.figure}: a wavelength segment of {len(ticks)} ticks"
            raise ValueError(msg)
        positions = tuple(
            _bottom_column(box, c, box.bottom(c) - PROBE) for c, _ in ticks
        )
        axis = Axis(positions, tuple(float(v) for _, v in ticks), tolerance=X_TOLERANCE)
        axis.check()
        segments.append(axis)
        ranges.append((positions[0], positions[-1]))
    del box
    return Scale(tuple(segments), tuple(ranges))


def log_scale(box: Frame, ticks: tuple[tuple[float, float], ...], side: str) -> Scale:
    """Return a logarithmic axis calibrated decade by decade, in left-axis rows."""
    line = box.left if side == "left" else box.right
    decades: dict[int, list[tuple[float, float]]] = {}
    for row, value in ticks:
        position = _left_row(box, line(row), row)
        exponent = math.floor(math.log10(value) + 1e-9)
        decades.setdefault(exponent, []).append((position, value))
        if value == 10.0**exponent:  # a decade tick belongs to the decade below too
            decades.setdefault(exponent - 1, []).append((position, value))
    segments, ranges = [], []
    for exponent in sorted(decades):
        points = sorted(decades[exponent])
        if len(points) < 3:
            continue
        axis = Axis(
            tuple(p for p, _ in points),
            tuple(v for _, v in points),
            log=True,
            tolerance=Y_TOLERANCE,
        )
        axis.check()
        slope, intercept = np.polyfit(axis.positions, np.log10(axis.values), 1)
        ends = sorted(
            ((exponent - intercept) / slope, (exponent + 1 - intercept) / slope)
        )
        segments.append(axis)
        ranges.append((float(ends[0]), float(ends[1])))
    return Scale(tuple(segments), tuple(ranges))


def _ratio_line(dark: np.ndarray, box: Frame, plot: Plot) -> Line:
    """Return the drawn R = 1.5 line, fitted like a frame line."""
    left, right = round(box.left(plot.ratio_line)), round(box.right(plot.ratio_line))
    return _frame_line(
        dark, plot.ratio_line, left + 30, right - 30, horizontal=True, thickest=13
    )


Trace = dict[int, tuple[int, int]]
"""A traced curve: per column, the first and last row of its run (or cluster)."""


def _clusters(column: np.ndarray, offset: int, join: int) -> list[tuple[int, int]]:
    """Return the dark runs of a column, those within ``join`` pixels merged.

    Runs shorter than ``SPECK`` pixels (specks of the scan) are dropped.
    """
    out: list[tuple[int, int]] = []
    for first, last in _spans(column):
        if last - first + 1 < SPECK:
            continue
        if out and first - out[-1][1] <= join:
            out[-1] = (out[-1][0], last)
        else:
            out.append((first, last))
    return [(offset + first, offset + last) for first, last in out]


def trace(
    dark: np.ndarray,
    seed: tuple[int, int],
    rows: tuple[np.ndarray, np.ndarray],
    columns: tuple[int, int],
    join: int = 0,
) -> Trace:
    """Trace a curve column by column, both ways from a seed pixel.

    In each column, between the rows ``rows`` gives for it (first, last; arrays
    indexed from ``columns[0]``), the dark run (or the cluster of runs within ``join``
    pixels) nearest to the run taken in the previous column, if it overlaps it or is
    within ``STEP`` pixels of it (ties: the longer overlap, then the nearer centre).
    A line is 8-connected, so on a steep stretch the run of each column touches the
    next one's. A column with no such run is skipped; after ``MAX_GAP`` skipped
    columns the trace stops.
    """
    first_column, last_column = columns
    found: Trace = {}
    for direction in (1, -1):
        previous = (seed[1], seed[1])
        c, skipped = seed[0], 0
        while first_column <= c <= last_column and skipped <= MAX_GAP:
            low = int(rows[0][c - first_column])
            high = int(rows[1][c - first_column])
            runs = _clusters(dark[low:high, c], low, join)
            centre = (previous[0] + previous[1]) / 2
            near = [
                (
                    max(0, a - previous[1], previous[0] - b),
                    -(min(b, previous[1]) - max(a, previous[0])),
                    abs((a + b) / 2 - centre),
                    (a, b),
                )
                for a, b in runs
            ]
            near = [item for item in near if item[0] <= STEP]
            if near:
                previous = min(near)[3]
                found[c] = previous
                skipped = 0
            else:
                skipped += 1
            c += direction
    return found


def _limits(
    box: Frame, columns: tuple[int, int], top: Line | None, below: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return per column the first and last row searched.

    Inside the frame (clear of its ticks), and below ``top`` (the R = 1.5 line) by
    ``below`` pixels.
    """
    cs = [float(c) for c in range(columns[0], columns[1] + 1)]
    first = np.array([box.top(c) + 2 * MIN_TICK + 15 for c in cs])
    if top is not None:
        first = np.array([top(c) + below for c in cs])
    last = np.array([box.bottom(c) - 2 * MIN_TICK - 15 for c in cs])
    return np.round(first).astype(int), np.round(last).astype(int)


def _sample(
    plot: Plot,
    box: Frame,
    traced: Trace,
    value: Scale,
    usable: Callable[[int, tuple[int, int]], bool],
    left_out: tuple[tuple[float, float], ...],
) -> list[Point]:
    """Sample a traced curve at the plot's wavelengths.

    For each wanted wavelength, the traced column whose wavelength (at the curve's
    row there) is nearest, if within half a column; the value at the centre of its
    run, if the run is usable and the value inside the axis's calibrated decades.
    """
    wavelength = x_scale(box, plot)
    columns = np.array(sorted(traced))
    centres = np.array([(traced[c][0] + traced[c][1]) / 2 for c in columns])
    lam = np.array(
        [
            np.nan if (w := wavelength(_bottom_column(box, c, r))) is None else w
            for c, r in zip(columns, centres, strict=True)
        ]
    )
    points: dict[str, Point] = {}
    for ticks, step in zip(plot.x_segments, plot.steps, strict=True):
        low, high = ticks[0][1], ticks[-1][1]
        count = round((high - low) / step)
        per_column = (high - low) / abs(
            _bottom_column(box, ticks[-1][0], box.bottom(ticks[-1][0]))
            - _bottom_column(box, ticks[0][0], box.bottom(ticks[0][0]))
        )
        for k in range(count + 1):
            wanted = round(low + k * step, 6)
            if any(a <= wanted <= b for a, b in left_out):
                continue
            row_id = f"{wanted:g}"
            if row_id in points:
                continue
            index = int(np.nanargmin(np.abs(lam - wanted)))
            if not abs(lam[index] - wanted) <= HALF_COLUMN * per_column:
                continue
            column = int(columns[index])
            run = traced[column]
            if not usable(column, run):
                continue
            centre = float(centres[index])
            read = value(_left_row(box, column, centre))
            if read is None:
                continue
            points[row_id] = Point(
                row_id,
                f"{lam[index]:.2f}",
                f"{read:.3e}",
                "",
                f"pixel ({column}, {centre:.1f})",
            )
    return sorted(points.values(), key=lambda p: float(p.x))


@dataclass(frozen=True)
class Reading:
    """A plot read from its image: the frame, scales, R = 1.5 line and traces."""

    plot: Plot
    box: Frame
    wavelength: Scale
    sigma: Scale
    ratio: Scale
    line: Line
    sigma_trace: Trace
    ratio_trace: Trace


def reading(key: str) -> Reading:
    """Return a plot's calibration and traced curves."""
    plot = PLOTS[key]
    dark = image(DIRECTORY, plot.page) < DARK
    box = frame(dark, plot)
    _check_ticks(dark, box, plot)
    line = _ratio_line(dark, box, plot)
    columns = (
        round(box.left(plot.frame[1])) + 2 * MIN_TICK,
        round(box.right(plot.frame[1])) - 2 * MIN_TICK,
    )
    sigma_trace: Trace = {}
    for seed in plot.sigma_seeds:
        sigma_trace |= trace(
            dark, seed, _limits(box, columns, line, RATIO_MARGIN), columns
        )
    ratio_trace: Trace = {}
    limits = _limits(box, columns, line, RATIO_TOP)
    limits = (limits[0], np.minimum(limits[1], limits[0] + RATIO_BAND))
    for seed in plot.ratio_seeds:
        ratio_trace |= trace(dark, seed, limits, columns, join=JOIN)
    return Reading(
        plot,
        box,
        x_scale(box, plot),
        log_scale(box, plot.left_ticks, "left"),
        log_scale(box, plot.right_ticks, "right"),
        line,
        sigma_trace,
        ratio_trace,
    )


def readings(key: str) -> dict[str, list[Point]]:
    """Return a plot's cross section and ratio, read from the image."""
    read = reading(key)
    plot, box, line = read.plot, read.box, read.line

    def narrow(column: int, run: tuple[int, int]) -> bool:
        # The run and its traced neighbours ``NEIGHBOURS`` columns either side are no
        # taller than MAX_RUN: not a steep flank, and not the tip of a resonance,
        # whose apex column alone can be short.
        around = [read.sigma_trace.get(column + k) for k in (-NEIGHBOURS, NEIGHBOURS)]
        return all(r is not None and r[1] - r[0] + 1 <= MAX_RUN for r in (run, *around))

    def clear(column: int, run: tuple[int, int]) -> bool:
        # Not touching the R = 1.5 line: the band's first row is the line's edge.
        touching = run[0] <= line(column) + RATIO_TOP + 0.5
        return run[1] - run[0] + 1 <= MAX_CLUSTER and not touching

    return {
        FILES[key, "sigma2"]: _sample(
            plot, box, read.sigma_trace, read.sigma, narrow, plot.sigma_left_out
        ),
        FILES[key, "ratio"]: _sample(
            plot, box, read.ratio_trace, read.ratio, clear, plot.ratio_left_out
        ),
    }


# The step in Fig. 4 where the 2P1/2 channel closes, 'at 1843 A' (p. 844), and how far,
# angstrom, the drawn step may be from it (two columns of Fig. 4's 100-angstrom
# divisions, 0.75 angstrom each, and the bottom ticks' residual).
STEP_PRINTED = 1843.0
STEP_MARGIN = 3.0
# How far the drawn R = 1.5 line may read from 1.5 on the right axis before the
# calibration is taken to have failed (it reads 1.50 to 1.71; the record's notes).
RATIO_LINE_MARGIN = 0.15
# McCown et al. (1982), Table I and p. 3323: McGuire's coefficient for Xe at 193 nm,
# 'sigma (circular polarization) is found to be 8 x 10^-32 cm4 W-1' (alpha, with
# sigma_2 = alpha hbar omega; mccown_1982), and its one significant figure.
MCCOWN_ALPHA = 8e-32
MCCOWN_MARGIN = 0.10
HC_EV_ANGSTROM = 12398.42  # hc, eV angstrom (CODATA 2018: 12398.419843...)
EV = 1.602176634e-19  # J


def _step(read: Reading) -> float:
    """Return the wavelength of the drawn step of Fig. 4 between 1820 and 1870 A.

    The first traced column, going to longer wavelengths, whose run is taller than
    ``MAX_RUN`` (the line drops there by 30 pixels within a column or two).
    """
    for column in sorted(read.sigma_trace):
        first, last = read.sigma_trace[column]
        centre = (first + last) / 2
        at = read.wavelength(_bottom_column(read.box, column, centre))
        if at is not None and 1820 <= at <= 1870 and last - first + 1 > MAX_RUN:
            return at
    msg = "Fig. 4: no step between 1820 and 1870 angstrom"
    raise ValueError(msg)


def checks(files: dict[str, list[Point]] | None = None) -> dict[str, object]:
    """Return the checks of the reading against the numbers printed for these curves.

    Raises if one fails its margin: the step of Fig. 4 (``STEP_MARGIN``), the drawn
    R = 1.5 lines (``RATIO_LINE_MARGIN``), McCown et al.'s quotation of the Xe
    cross section at 193 nm (``MCCOWN_MARGIN``).
    """
    files = extract(check=False) if files is None else files
    out: dict[str, object] = {}
    failures = []
    xenon = reading("fig4")
    step = _step(xenon)
    out["Xe step (angstrom)"] = (step, STEP_PRINTED)
    if abs(step - STEP_PRINTED) > STEP_MARGIN:
        failures.append(f"Fig. 4's step at {step:.1f} angstrom, printed 1843")
    for key in PLOTS:
        read = xenon if key == "fig4" else reading(key)
        box = read.box
        ends = (box.left(box.bottom(0)) + 100, box.right(box.bottom(0)) - 100)
        values = [
            read.ratio(_left_row(box, c, read.line(c))) for c in np.linspace(*ends, 9)
        ]
        lows = [v for v in values if v is not None]
        out[f"{read.plot.figure} R = 1.5 line reads"] = (min(lows), max(lows))
        if max(abs(v / 1.5 - 1) for v in lows) > RATIO_LINE_MARGIN:
            failures.append(f"{read.plot.figure}: the R = 1.5 line reads {lows}")
    # McCown et al.: 8e-32 cm4/W at 193 nm, times hbar omega (1930 angstrom).
    points = {p.row: p for p in files[FILES["fig4", "sigma2"]]}
    read_1930 = float(points["1930"].value)
    quoted = MCCOWN_ALPHA * HC_EV_ANGSTROM / float(points["1930"].x) * EV
    out["Xe at 1930 angstrom (cm4 s): read, McCown's quotation"] = (read_1930, quoted)
    if abs(read_1930 / quoted - 1) > MCCOWN_MARGIN:
        failures.append(
            f"Xe at 1930 angstrom reads {read_1930:.3e}, quoted {quoted:.3e}"
        )
    if failures:
        raise ValueError("; ".join(failures))
    return out


def extract(check: bool = True) -> dict[str, list[Point]]:
    """Return the cross sections and ratios read from Figs. 1-4, checked."""
    files: dict[str, list[Point]] = {}
    for key in PLOTS:
        files |= readings(key)
    if check:
        checks(files)
    return files


def main() -> None:
    """Write the data files and print the checks."""
    files = extract(check=False)
    for label, value in checks(files).items():
        print(f"{label}: {value}")
    write(files)


if __name__ == "__main__":
    main()

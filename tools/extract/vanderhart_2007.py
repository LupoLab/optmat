# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract van der Hart (2007): Ar at 390 nm, time-dependent R-matrix rates.

**Table I** (p. 4, ``tbl_2c7127524ab5``): the ionisation rates of Ar (s^-1) at 0.48,
0.75 and 1.50e14 W/cm^2, with one Ar+ target state (3s2 3p5 2P) and with two (adding
3s 3p6 2S). Read from the library's cells, which print 'a x 10^b' as ``a``, a
backspace character and ``10b`` (``9.79\x081010`` for 9.79 x 10^10, checked against
the table image): the script checks that shape and reads the exponent after the 10.

**Fig. 4** (p. 3, library figure ``fig_3c19189b753a``): the same rates at 12 intensities
(solid circles) and the R-matrix Floquet rates of van der Hart 2006 (the solid line).
The figure is a vector drawing, read from the PDF (data rule 8, TC-699, TC-700):

* the axes: the x ticks, 13 at 0.25e14 W/cm^2 steps, are the vertices of the frame
  path (52 segments); the y ticks, 1e8 to 1e14 s^-1 in two-decade steps, are the
  vertices of the 16-segment tick path. Every tick is checked against the fitted
  scale.
* the time-dependent rates: the 13 filled circles (a 13-segment path, drawn twice, as
  fill and outline, 2.8 and 3.4 pt wide), less the legend's (at x = 456.95 pt). Their
  centres. The three intensities of Table I are left out of the dataset (their numbers
  are printed) and serve as the check of the reading: figure over table 1.013, 1.057,
  1.003.
* the R-matrix Floquet curve: the 1221-vertex polyline, every vertex.

The ADK curve (the dash-dot line) is not taken.
"""

import re

from _figure import Axis, paths, unique
from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "VanderHart_2007_TimeDependentMultielectron"
FIGURE_PAGE = 3
REGION = (330.0, 140.0, 540.0, 290.0)
LEGEND_X = 456.949
TABLE_INTENSITIES = (0.48, 0.75, 1.50)
X_TICKS = tuple(0.25 * k for k in range(13))
# The frame's bottom and top edges, pt (the 52-segment frame path's extreme y).
FRAME_BOTTOM, FRAME_TOP = 153.093, 287.920
Y_TICKS = (1e8, 1e10, 1e12, 1e14)


def _axes(drawn: list) -> tuple[Axis, Axis]:
    """Return the intensity (1e14 W/cm^2) and rate (s^-1) axes of Fig. 4."""
    frame = [d for d in drawn if len(d.vertices) == 52]
    ticks = [d for d in drawn if len(d.vertices) == 16]
    if len(frame) != 1 or len(ticks) != 1:
        msg = "Fig. 4: the frame or the tick path is not where expected"
        raise ValueError(msg)
    xs = sorted({round(float(x), 3) for x in frame[0].vertices[:, 0]})
    ys = sorted({round(float(y), 3) for y in ticks[0].vertices[:, 1]})
    if len(xs) != len(X_TICKS) or len(ys) != len(Y_TICKS):
        msg = f"Fig. 4: {len(xs)} x ticks and {len(ys)} y ticks"
        raise ValueError(msg)
    x_axis = Axis(tuple(xs), X_TICKS)
    y_axis = Axis(tuple(ys), Y_TICKS, log=True)
    x_axis.check()
    y_axis.check()
    return x_axis, y_axis


def _figure_points() -> tuple[list[Point], list[Point], list[tuple[float, float]]]:
    """Return Fig. 4's time-dependent points, its Floquet curve and the checks."""
    drawn = paths(DIRECTORY, FIGURE_PAGE, REGION)
    x_axis, y_axis = _axes(drawn)
    markers = unique([d.centre for d in drawn if len(d.vertices) == 13 and d.size < 4])
    markers = [m for m in markers if abs(m[0] - LEGEND_X) > 0.01]
    if len(markers) != 12:
        msg = f"Fig. 4: {len(markers)} markers"
        raise ValueError(msg)
    rates, checks = [], []
    for x, y in markers:
        intensity = float(x_axis(x))
        rate = float(y_axis(y))
        near = [t for t in TABLE_INTENSITIES if abs(t - intensity) < 1e-3]
        if near:
            checks.append((near[0], rate))
            continue
        printed = f"marker at ({x:.3f}, {y:.3f}) pt"
        rates.append(
            Point(f"{intensity:.3f}", f"{intensity:.3f}e14", f"{rate:.4e}", "", printed)
        )
    curves = [d for d in drawn if len(d.vertices) == 1221]
    if len(curves) != 1:
        msg = "Fig. 4: the Floquet curve is not where expected"
        raise ValueError(msg)
    floquet = []
    # The frame clips the curve: keep the vertices it shows.
    frame = x_axis.positions[0], x_axis.positions[-1], FRAME_BOTTOM, FRAME_TOP
    for k, (x, y) in enumerate(curves[0].vertices, start=1):
        if not (frame[0] <= x <= frame[1] and frame[2] <= y <= frame[3]):
            continue
        intensity, rate = float(x_axis(x)), float(y_axis(y))
        floquet.append(
            Point(
                f"v{k:04d}",
                f"{intensity:.5f}e14",
                f"{rate:.4e}",
                "",
                f"vertex {k} at ({x:.3f}, {y:.3f}) pt",
            )
        )
    return rates, floquet, checks


def figure_check() -> list[tuple[float, float, float]]:
    """Return Table I's intensities, its one-state rate and Fig. 4's reading of it."""
    _, _, checks = _figure_points()
    table = {float(p.row): float(p.value) for p in _table()["one"]}
    return [(t, table[t], rate) for t, rate in sorted(checks)]


_POWER = re.compile(r"^(\d+\.\d+)\x0810(\d+)$")


def _table() -> dict[str, list[Point]]:
    """Return Table I's rates with one and with two Ar+ target states."""
    grid = json_table(DIRECTORY, "tbl_2c7127524ab5")
    if grid[:2] != [
        ["Intensity", "One-target", "Two-target"],
        ["1014 W/cm2", "s-1", "s-1"],
    ]:
        msg = f"unexpected header of Table I: {grid[:2]!r}"
        raise ValueError(msg)
    out: dict[str, list[Point]] = {"one": [], "two": []}
    for cells in grid[2:]:
        intensity = number(cells[0])
        printed = " | ".join(cells).replace("\x08", " x ")
        for name, column in (("one", 1), ("two", 2)):
            out[name].append(
                Point(intensity, f"{intensity}e14", _rate(cells[column]), "", printed)
            )
    return out


def _rate(cell: str) -> str:
    """Return a rate cell, 'a', backspace, '10b' for a x 10^b, as 'aeb'."""
    match = _POWER.match(cell)
    if match is None:
        msg = f"unexpected rate cell of Table I: {cell!r}"
        raise ValueError(msg)
    return f"{number(match.group(1))}e{match.group(2)}"


# The figure's circles against Table I: how far the reading may stray (it is 1.3, 5.6
# and 0.3 %).
AGREEMENT = 0.06


def extract() -> dict[str, list[Point]]:
    """Return Table I's rates and Fig. 4's time-dependent points and Floquet curve."""
    table = _table()
    for intensity, printed, read in figure_check():
        if abs(read / printed - 1) > AGREEMENT:
            msg = f"Fig. 4 reads {read:.3e} at {intensity}e14, Table I {printed:.3e}"
            raise ValueError(msg)
    rates, floquet, _ = _figure_points()
    return {
        "Ar/vanderhart_2007__table1_one_state_Ar.csv": table["one"],
        "Ar/vanderhart_2007__table1_two_states_Ar.csv": table["two"],
        "Ar/vanderhart_2007__fig4_td_rmatrix_Ar.csv": rates,
        "Ar/vanderhart_2007__fig4_floquet_Ar.csv": floquet,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

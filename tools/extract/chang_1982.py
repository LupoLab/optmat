# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chang and Kim (1982): photoionisation of the excited p states of Ne-Xe.

The library copy is a scan of Phys. Rev. A 26, 2728 (1982) (PDF page = journal page -
2727), each page one embedded image (4938 x 6686 pixels) with a text layer.

* **Table II** (p. 2731, ``tbl_c2bfbdfde088``, JSON only): the calculated (length
  form, p. 2732) cross sections of the first excited p state of each gas, in Mb, at
  308 nm (Ne, Xe) or 248 nm (Ar, Kr), beside the values of Duzy and Hyman ('taken
  from Ref. 7'). Each is one file of one point; the printed row is kept in
  ``printed``. Checked against the table image.
* **Table I** (p. 2730, ``tbl_a1b5c6b30139``): the calculated effective quantum
  numbers n* of the excited mp orbitals, from the orbital eigenvalues of the
  effective Hamiltonian (p. 2731). Not data here: the binding energy of a state is
  Ry / n*^2, which turns the figures' photoelectron energy into the photon energy.
* **Figs. 1-4** (pp. 2729-2730): the velocity-form cross sections (Mb, log axis,
  1e-4 to 1e2) of the 3p-6p (Ne), 4p-7p (Ar), 5p-8p (Kr) and 6p-9p (Xe) states
  against the photoelectron energy (eV, log axis, 0.1 to the frame at about 400 eV),
  four curves each. Raster scans: read from the page image at native resolution.

The figures are read as follows (:func:`curve`):

* the axes: the frame's sides are hand-drafted and the scan is turned by up to 0.3
  degree, so each side has its own logarithmic :class:`Axis`, calibrated on its
  decade ticks (written below: the pixel row or column of each labelled decade, read
  from the image), and a point's value is interpolated between the two sides by its
  position across the frame. Every decade tick is checked to be a dark run of its
  side's scan line (``SCAN`` pixels inside the frame) within ``FOUND`` pixels, every
  fit within ``TOLERANCE``, and every minor tick (2-9 x 10^k) is checked to be a dark
  run within ``MINOR`` pixels of the place its two decades give it (``_check_ticks``).
  The right frame line (about 400 eV) is not labelled and is not used; in Fig. 1 a
  short tick at about 356 eV, between the 300 eV tick and the frame, is on no value
  and is ignored.
* the curves: in a column, the dark runs between ``MARGIN`` pixels inside the top and
  bottom frame lines, runs closer than ``JOIN`` pixels joined (the scan breaks lines)
  and specks shorter than ``SPECK`` left out; where at least four runs are found, the
  lowest four are the four curves, in the order of their labels (a figure's own
  label, 'Kr' say, lies above them all). A value is the centre of the run (for a
  steep line, the centre of its vertical chord). The trace starts at the figure's
  ``seed`` energy and goes both ways one column at a time; it stops where fewer than
  four runs are found (the curves touch or cross) or the curve's row jumps by more
  than ``STEP`` pixels, and at the figure's ``last`` energy on the right (before the
  curves' end labels). The range read is therefore the range in which the four curves
  are resolved; below it, near threshold, they cross (Ne) or merge.
* the curves written (``TRACED``): the first excited p state of each gas and Kr 6p,
  one point every ``EVERY`` columns of the trace (about 0.06 decade); ``printed``
  holds the pixel and the photoelectron energy read.

The photon energy of a point is its photoelectron energy plus the binding energy of
the state from Table I, Ry / n*^2 (Ry = 13.605693 eV, R_infinity hc, CODATA 2018): the
calculation's own threshold, the conversion the record states. One pixel is 0.0028
decade in either axis (0.65 %); the lines are 6-10 pixels wide, the decade ticks lie
within 5 pixels of their fits and the minor ticks within 10 of their places.
"""

import math
import re
from dataclasses import dataclass
from functools import cache

import numpy as np
from _figure import Axis, image
from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Chang_1982_PhotoionizationExcitedStates"
TABLE_I = "tbl_a1b5c6b30139"
TABLE_II = "tbl_c2bfbdfde088"
RYDBERG_EV = 13.605693122994
"""R_infinity hc in eV (CODATA 2018); the paper's energies are in Rydberg units."""
HC_EV_NM = 1239.841984
"""hc in eV nm (CODATA 2018)."""
GASES = {"Ne": "3p", "Ar": "4p", "Kr": "5p", "Xe": "6p"}
"""Each gas's first excited p state, Table II's rows and the curves taken."""
TRACED = (*GASES.items(), ("Kr", "6p"))
"""The curves written: each gas's first excited p state, and Kr 6p (the text's)."""
STATE = re.compile(r"^(\d)p\^\{5\}(\d)p$")
"""Table II's state cells, ``4p^{5}5p``."""

DARK = 128
SCAN = 12
"""How far inside a frame line, pixels, the ticks are looked for."""
JOIN_TICKS = 4
"""Runs of a tick scan line closer than this, pixels, are one tick (drawn double)."""
JOIN = 6
"""Runs of a column closer than this, pixels, are one line (the scan breaks lines)."""
FOUND = 3.0
"""The farthest, pixels, a written decade tick may be from the centre of a dark run."""
TOLERANCE = 5.0
"""The largest residual, pixels, of a decade tick on its side's fit (hand-drafted)."""
MINOR = 10.0
"""The farthest, pixels, a minor tick predicted by the fit may be from a dark run."""
MARGIN = 45
"""Rows this close to the top or bottom frame line are not searched (the ticks)."""
EDGE = 60
"""Columns this close to the left or right frame line are not searched."""
SPECK = 4
"""Runs shorter than this, pixels, are specks of the scan, not lines."""
STEP = 25
"""The largest jump, pixels, of a curve's row from one column to the next."""
EVERY = 20
"""A point is written every this many columns of the trace."""
CURVES = 4
"""Curves per figure."""


@dataclass(frozen=True)
class Frame:
    """A figure: page, frame lines (approximate), decade ticks per side, figure id."""

    figure: str
    page: int
    columns: tuple[int, int]
    """The left and right frame lines' columns (at mid-height)."""
    rows: tuple[int, int]
    """The top and bottom frame lines' rows (at mid-width)."""
    left: tuple[float, ...]
    """Rows of the decades 1e2, 1e1, ..., 1e-4 Mb on the left side."""
    right: tuple[float, ...]
    top: tuple[float, ...]
    """Columns of 0.1, 1, 10 and 100 eV along the top side."""
    bottom: tuple[float, ...]
    states: tuple[str, ...]
    """The curves from the top one down, as labelled (at high energy)."""
    seed: float
    """The photoelectron energy (eV) the trace starts at."""
    last: float
    """The photoelectron energy (eV) the trace stops at on the right."""


Y_DECADES = (1e2, 1e1, 1.0, 1e-1, 1e-2, 1e-3, 1e-4)
X_DECADES = (0.1, 1.0, 10.0, 100.0)
MINOR_Y = tuple(m * 10.0**k for k in range(-4, 2) for m in range(2, 10))
"""The minor ticks of the sigma axes, 2e-4 to 9e1 Mb."""
MINOR_X = tuple(m * 10.0**k for k in range(-1, 3) for m in range(2, 10))
"""The minor ticks of the energy axes, 0.2 to 300 eV (the frame is at about 400)."""

# Decade ticks read from the page images (centres of the dark runs SCAN pixels inside
# each side); the frame lines are the 1e2 and 1e-4 Mb ticks and the 0.1 eV tick.
FIGURES = {
    "Ne": Frame(
        figure="fig_b8e5ffd05efa",
        page=2,
        columns=(924, 2217),
        rows=(826, 2958),
        left=(826.0, 1180.5, 1531.5, 1892.5, 2249.0, 2601.5, 2958.0),
        right=(825.5, 1184.5, 1534.5, 1897.5, 2251.0, 2607.0, 2957.0),
        top=(924.0, 1284.5, 1641.0, 2002.0),
        bottom=(924.5, 1285.0, 1642.5, 1999.5),
        states=("3p", "4p", "5p", "6p"),
        seed=20.0,
        last=58.0,
    ),
    "Ar": Frame(
        figure="fig_d4adf2172464",
        page=2,
        columns=(3002, 4290),
        rows=(831, 2962),
        left=(831.5, 1184.5, 1539.0, 1896.0, 2252.0, 2606.5, 2960.0),
        right=(836.5, 1192.5, 1546.5, 1903.5, 2260.0, 2616.0, 2965.0),
        top=(3005.5, 3362.5, 3719.5, 4084.5),
        bottom=(2995.5, 3351.0, 3708.0, 4068.0),
        states=("4p", "5p", "6p", "7p"),
        seed=5.0,
        last=55.0,
    ),
    "Kr": Frame(
        figure="fig_0aeec97192fc",
        page=3,
        columns=(916, 2201),
        rows=(845, 2974),
        left=(847.0, 1199.0, 1555.0, 1909.0, 2266.0, 2619.5, 2975.0),
        right=(842.0, 1194.5, 1549.5, 1903.5, 2260.5, 2615.0, 2973.0),
        top=(916.5, 1267.5, 1622.0, 1981.5),
        bottom=(921.5, 1275.5, 1633.5, 1989.0),
        states=("5p", "6p", "7p", "8p"),
        seed=10.0,
        last=58.0,
    ),
    "Xe": Frame(
        figure="fig_26a28978918d",
        page=3,
        columns=(912, 2203),
        rows=(3518, 5649),
        left=(3521.5, 3881.5, 4236.0, 4592.0, 4947.0, 5301.0, 5649.0),
        right=(3514.5, 3874.5, 4229.0, 4588.5, 4942.0, 5298.0, 5647.5),
        top=(910.5, 1264.5, 1624.5, 1984.5),
        bottom=(920.5, 1278.0, 1636.5, 1996.0),
        states=("6p", "7p", "8p", "9p"),
        seed=10.0,
        last=58.0,
    ),
}


def _spans(line: np.ndarray, join: int) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True, joining near runs."""
    index = np.flatnonzero(line)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > join) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _axis(positions: tuple[float, ...], values: tuple[float, ...]) -> Axis:
    """Return a logarithmic axis calibrated on a side's decade ticks."""
    return Axis(positions, values, log=True, tolerance=TOLERANCE)


@dataclass(frozen=True)
class Calibrated:
    """A figure's four side axes and its frame, for reading a pixel."""

    frame: Frame
    left: Axis
    right: Axis
    top: Axis
    bottom: Axis

    def value(self, column: float, row: float) -> tuple[float, float]:
        """Return (photoelectron energy in eV, sigma in Mb) at a pixel."""
        left, right = self.frame.columns
        top, bottom = self.frame.rows
        across = (column - left) / (right - left)
        down = (row - top) / (bottom - top)
        log_e = (1 - down) * math.log10(float(self.top(column))) + down * math.log10(
            float(self.bottom(column))
        )
        log_s = (1 - across) * math.log10(float(self.left(row))) + across * math.log10(
            float(self.right(row))
        )
        return 10**log_e, 10**log_s

    def column(self, energy: float) -> int:
        """Return the column of a photoelectron energy at mid-height."""
        top = np.log10(float(self.top(0.0))), np.log10(float(self.top(1.0)))
        bottom = np.log10(float(self.bottom(0.0))), np.log10(float(self.bottom(1.0)))
        intercept = (top[0] + bottom[0]) / 2
        slope = (top[1] - top[0] + bottom[1] - bottom[0]) / 2
        return round((math.log10(energy) - intercept) / slope)


def _check_ticks(dark: np.ndarray, frame: Frame) -> None:
    """Raise unless the decade and minor ticks are dark runs where they should be."""
    left, right = frame.columns
    top, bottom = frame.rows
    sides = {
        "left": (frame.left, Y_DECADES, dark[:, left + SCAN]),
        "right": (frame.right, Y_DECADES, dark[:, right - SCAN]),
        "top": (frame.top, X_DECADES, dark[top + SCAN, :]),
        "bottom": (frame.bottom, X_DECADES, dark[bottom - SCAN, :]),
    }
    for side, (ticks, values, line) in sides.items():
        centres = np.array([(a + b) / 2 for a, b in _spans(line, JOIN_TICKS)])
        axis = _axis(ticks, values)
        axis.check()
        inner = ticks[1:] if side in {"top", "bottom"} else ticks[1:-1]
        for position in inner:
            if np.min(np.abs(centres - position)) > FOUND:
                msg = f"{frame.figure}: no {side} decade tick at {position}"
                raise ValueError(msg)
        order = np.argsort(np.log10(values))
        logs = np.log10(values)[order]
        positions = np.asarray(ticks)[order]
        for wanted in MINOR_Y if side in {"left", "right"} else MINOR_X:
            if wanted < 10 ** logs[0] or wanted > 10 ** logs[-1] * 3:
                continue
            if wanted > 10 ** logs[-1]:  # beyond the last decade: its spacing
                spacing = (positions[-1] - positions[-2]) / (logs[-1] - logs[-2])
                guess = positions[-1] + spacing * (math.log10(wanted) - logs[-1])
            else:  # between the two decades either side, as drafted
                guess = float(np.interp(math.log10(wanted), logs, positions))
            if np.min(np.abs(centres - guess)) > MINOR:
                msg = f"{frame.figure}: no {side} minor tick near {guess:.0f}"
                raise ValueError(msg)


@cache
def calibrate(gas: str) -> tuple[Calibrated, np.ndarray]:
    """Return a figure's calibration (its ticks checked) and its dark pixels."""
    frame = FIGURES[gas]
    dark = image(DIRECTORY, frame.page) < DARK
    _check_ticks(dark, frame)
    calibrated = Calibrated(
        frame,
        _axis(frame.left, Y_DECADES),
        _axis(frame.right, Y_DECADES),
        _axis(frame.top, X_DECADES),
        _axis(frame.bottom, X_DECADES),
    )
    return calibrated, dark


def _runs(dark: np.ndarray, frame: Frame, column: int) -> list[tuple[int, int]]:
    """Return the dark runs (first and last row) of a column inside the margins.

    Specks (runs shorter than ``SPECK``) are left out.
    """
    top, bottom = frame.rows
    first = top + MARGIN
    return [
        (first + a, first + b)
        for a, b in _spans(dark[first : bottom - MARGIN, column], JOIN)
        if b - a + 1 >= SPECK
    ]


@cache
def curve(gas: str, state: str) -> list[tuple[int, float, float, float]]:
    """Return a state's curve: (column, row, photoelectron eV, sigma Mb) per column.

    The ``CURVES`` lowest runs of a column are the curves in label order; the trace
    runs from the seed both ways while there are that many runs and the row moves by
    at most ``STEP`` pixels per column.
    """
    calibrated, dark = calibrate(gas)
    frame = calibrated.frame
    rank = frame.states.index(state) - CURVES
    left, right = frame.columns
    stop = min(calibrated.column(frame.last), right - EDGE)
    found: dict[int, tuple[int, int]] = {}
    seed = calibrated.column(frame.seed)
    for direction in (1, -1):
        column, previous = seed, None
        while left + EDGE <= column <= stop:
            runs = _runs(dark, frame, column)
            if len(runs) < CURVES:
                break
            a, b = runs[rank]
            if previous is not None and abs((a + b) / 2 - previous) > STEP:
                break
            found[column] = (a, b)
            previous = (a + b) / 2
            column += direction
    out = []
    for column in sorted(found):
        a, b = found[column]
        energy, sigma = calibrated.value(column, (a + b) / 2)
        out.append((column, (a + b) / 2, energy, sigma))
    return out


def _at(trace: list[tuple[int, float, float, float]], energy: float) -> float:
    """Return a curve's sigma at a photoelectron energy (log-log interpolation)."""
    log_e = np.log10([p[2] for p in trace])
    log_s = np.log10([p[3] for p in trace])
    if not log_e.min() <= math.log10(energy) <= log_e.max():
        msg = f"{energy:.3g} eV is outside the curve read"
        raise ValueError(msg)
    return float(10 ** np.interp(math.log10(energy), log_e, log_s))


def n_star() -> dict[tuple[str, str], str]:
    """Return Table I's effective quantum numbers by (gas, state), as printed."""
    rows = table(DIRECTORY, TABLE_I)
    gases = rows[0][1:]
    out = {}
    for cells in rows[1:]:
        for gas, cell in zip(gases, cells[1:], strict=True):
            if cell:
                out[(gas, cells[0])] = cell
    if len(out) != 16 or list(gases) != list(GASES):
        msg = "Table I: cells moved"
        raise ValueError(msg)
    return out


def binding(gas: str, state: str) -> float:
    """Return a state's binding energy (eV) from Table I's n*: Ry / n*^2."""
    return RYDBERG_EV / float(n_star()[(gas, state)]) ** 2


def table_ii() -> dict[str, tuple[str, str, str, list[str]]]:
    """Return Table II by gas: wavelength (nm), present work, Duzy and Hyman, row."""
    rows = json_table(DIRECTORY, TABLE_II)
    header = [rows[0][2], rows[0][4], *rows[1][2:]]
    if header != [
        r"\lambda = 308 \text{ nm}",
        r"\lambda = 248 \text{ nm}",
        "Present work",
        "Duzy and Hyman",
        "Present work",
        "Duzy and Hyman",
    ]:
        msg = f"Table II: headers moved {header}"
        raise ValueError(msg)
    out = {}
    for cells in rows[2:]:
        match = STATE.match(cells[1])
        if cells[0] not in GASES or match is None or f"{match[2]}p" != GASES[cells[0]]:
            msg = f"Table II: an unexpected row {cells}"
            raise ValueError(msg)
        filled = [i for i in (2, 4) if cells[i]]
        if len(filled) != 1 or not cells[filled[0] + 1]:
            msg = f"Table II: {cells[0]} has no single wavelength"
            raise ValueError(msg)
        i = filled[0]
        wavelength = "308" if i == 2 else "248"
        out[cells[0]] = (wavelength, cells[i], cells[i + 1], cells)
    return out


TEXT_KR_6P = (193.5, 7.4e-19)
"""Kr 6p (J = 0), velocity form, 'at lambda = 193.5 nm' (p. 2731), nm and cm^2."""
KROLL = {248.0: 3.7, 193.0: 1.8}
"""kroll_1990's Table II, 'ref. 14' (this paper) for Xe 2p5, 1e-18 cm^2, by nm."""
KR_6P_RATIO = 1.10
"""Fig. 3's 6p curve at 193.5 nm over the text's 7.4e-19 cm^2: 10 % (the register)."""
GEOHEGAN = 8.0e-19
"""geohegan_1986's 'Chang and Kim' value for Kr 6p[3/2]2, cm^2 (p. 273)."""
LENGTH_VELOCITY = 0.10
"""'the length and velocity results generally agree with each other to 10% or less
near the ionization threshold (to about 5 eV above)' for Ne, Ar and Kr (p. 2730)."""


def checks() -> dict[str, object]:
    """Return the figures read against the printed numbers.

    Raises unless the Kr 6p curve gives ``KR_6P_RATIO`` times the text's 7.4e-19 cm^2
    at 193.5 nm within 3 %, and unless the Ar and Kr curves (velocity) give Table
    II's length values within 10 % plus 5 % for the reading (the Ne value lies where
    the curves cross and is not read; Xe's forms differ more).
    """
    out: dict[str, object] = {}
    printed = table_ii()
    for gas, state in GASES.items():
        trace = curve(gas, state)
        wavelength, value, _, _ = printed[gas]
        electron = HC_EV_NM / float(wavelength) - binding(gas, state)
        out[f"{gas} {state}: range read (photoelectron eV)"] = (
            round(trace[0][2], 3),
            round(trace[-1][2], 1),
        )
        if electron < trace[0][2]:  # Ne 3p at 308 nm: where the curves cross
            out[f"{gas} {state} at {wavelength} nm"] = f"{electron:.2f} eV, not read"
            continue
        read = _at(trace, electron)
        ratio = read / float(value)
        out[f"{gas} {state} at {wavelength} nm: Fig. (velocity) / Table II"] = round(
            ratio, 3
        )
        if gas != "Xe" and abs(ratio - 1) > LENGTH_VELOCITY + 0.05:
            msg = f"{gas}: the figure is {ratio:.2f} x Table II"
            raise ValueError(msg)
    wavelength, value = TEXT_KR_6P
    electron = HC_EV_NM / wavelength - binding("Kr", "6p")
    read = _at(curve("Kr", "6p"), electron) * 1e-18
    out["Kr 6p at 193.5 nm: Fig. / text 7.4e-19"] = round(read / value, 3)
    out["Kr 6p at 193.5 nm: geohegan_1986's 8.0e-19 / Fig."] = round(GEOHEGAN / read, 3)
    if abs(read / value - KR_6P_RATIO) > 0.03:
        msg = f"Kr 6p: the figure gives {read:.3g} cm2, the text 7.4e-19"
        raise ValueError(msg)
    xe = curve("Xe", "6p")
    for nm, quoted in KROLL.items():
        sigma = _at(xe, HC_EV_NM / nm - binding("Xe", "6p"))
        out[f"Xe 6p at {nm:.0f} nm (Mb), Fig.; kroll_1990 quotes {quoted}"] = round(
            sigma, 2
        )
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table II's values and the curves of ``TRACED``."""
    checks()
    files: dict[str, list[Point]] = {}
    for gas, (wavelength, present, duzy, cells) in table_ii().items():
        state = GASES[gas]
        row = " | ".join(cells)
        files[f"{gas}/chang_1982__table2_{state}.csv"] = [
            Point(wavelength, wavelength, present, "", row)
        ]
        files[f"{gas}/chang_1982__table2_duzy_hyman_{state}.csv"] = [
            Point(wavelength, wavelength, duzy, "", row)
        ]
    for gas, state in TRACED:
        threshold = binding(gas, state)
        trace = curve(gas, state)
        points = []
        for column, row, energy, sigma in trace[::EVERY]:
            points.append(
                Point(
                    f"{column}",
                    f"{energy + threshold:.4f}",
                    f"{sigma:.3e}",
                    "",
                    f"pixel ({column}, {row:.1f}); photoelectron {energy:.4g} eV",
                )
            )
        files[f"{gas}/chang_1982__fig_{state}.csv"] = points
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

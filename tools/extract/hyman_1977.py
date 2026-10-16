# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hyman (1977): photoionisation cross sections of excited Ar and Kr at 2486 A.

The library copy is AIP's scan of Appl. Phys. Lett. 31, 14 (1977) behind a cover page
(PDF p. 2 = journal p. 14, PDF p. 3 = p. 15), each page one embedded image (300 dpi)
with a text layer.

* **Table I** (p. 15, ``tbl_ae5d166a6610``): the cross sections of Ar(4s), Ar(4p),
  Kr(5s) and Kr(5p) at 2486 A, in cm^2, one value each. The library's cells are an OCR
  of the scan (``4. 3x 10-18``); each is read by ``CELL`` as mantissa and power of ten,
  and was checked against the table image and the rendered page. The text repeats two
  of them ('~4.5 x 10^-18 cm2' in the abstract; 'the calculated Ar (4p) value of
  4.3 x 10^-18 cm2', p. 15).

The figures are not data: Figs. 2 (Ar, p. 14) and 3 (Kr, p. 15) draw the cross
sections against the photoelectron's excess energy epsilon (Ry), and the paper prints
no ionisation energy of the averaged levels (Fig. 1 draws the unaveraged levels on an
energy axis and prints none), so no photon energy can be given to them.
:func:`checks` reads the p-state curves from the page images at their native
resolution to locate Table I on them and to test dakka_2018's quotation:

* the axes: the tick marks (inward, both sides of the frame) are written below,
  (pixel, value), and each is checked to be a dark run of the image within ``FOUND``
  pixels. The figures are hand-drafted: the x ticks at the top and bottom of a frame
  differ by up to 4.5 pixels, the y ticks on its left and right by up to 3, so each
  side has its own :class:`Axis` (every tick within ``TOLERANCE`` of its fit) and a
  point's value is interpolated between the two sides by its position across the
  frame.
* the curve: traced column by column from a seed (the topmost dark run at
  ``SEED_EPSILON``, where the p curve is the top one), following the run nearest the
  previous position within ``STEP`` pixels (the curve label lies above the curve and is
  not reached); a value is the centre of the run in a column.
* Table I on the curve: the excess energy at which the curve has Table I's value, by
  log-linear interpolation between columns, and with the photon energy at 2486 A in
  the paper's Rydberg ('1 Ry = 13.6 eV', p. 14), the implied binding energy of the
  level. For Kr(5p) the curve is then read 0.0573 Ry higher, at the photon energy of
  215 nm, where dakka_2018 quotes 3.3 x 10^-18 cm^2 'Calculated from' this paper.

One pixel is 0.002 Ry in epsilon and 1.2 % in sigma; the lines are 3-4 pixels wide.
"""

import math
import re
from dataclasses import dataclass

import numpy as np
from _figure import Axis, image
from _library import Point, table, write

DIRECTORY = "Hyman_1977_PhotoionizationCrossSections"
TABLE_I = "tbl_ae5d166a6610"
STATES = {
    "Ar(4s)": ("Ar", "4s"),
    "Ar(4p)": ("Ar", "4p"),
    "Kr(5s)": ("Kr", "5s"),
    "Kr(5p)": ("Kr", "5p"),
}
"""Table I's rows, in order: material and level."""
CELL = re.compile(r"^(\d)\.\s?(\d)\s?x\s?10-(\d{2})$")
"""An OCR'd cell, ``4. 3x 10-18``: mantissa digits and the (negative) power of ten."""
WAVELENGTH_A = "2486"
"""Table I's wavelength, 'lambda = 2486 A' (caption)."""
RYDBERG_EV = 13.6
"""The paper's conversion, '1 Ry = 13.6 eV' (p. 14)."""
HC_EV_NM = 1239.841984
"""hc in eV nm (CODATA 2018), for the photon energies of 2486 A and 215 nm."""
DAKKA_NM = 215.0
"""dakka_2018's wavelength for its 'Literature' sigma_ei ('215 nm')."""
DAKKA_VALUE = 3.3e-18
"""dakka_2018's 'Calculated from' this paper (Table I, 'Literature'), cm^2."""

DARK = 128
FOUND = 1.5
"""The farthest, pixels, a written tick may be from the centre of a dark run."""
TOLERANCE = 3.0
"""The largest residual, pixels, of a tick on its side's fitted scale (hand-drafted)."""
STEP = 8
MAX_WIDTH = 8
"""The widest run, pixels, a value is read from (the line is 3-4 pixels wide)."""
SEED_EPSILON = 0.5


@dataclass(frozen=True)
class Frame:
    """A figure's tick marks, (pixel, value), on each side, and its p curve."""

    page: int
    left: tuple[tuple[float, float], ...]
    right: tuple[tuple[float, float], ...]
    top: tuple[tuple[float, float], ...]
    bottom: tuple[tuple[float, float], ...]
    columns: tuple[int, int]
    """The left and right frame lines' columns (for the y ticks' scan lines)."""
    rows: tuple[int, int]
    """The top and bottom frame lines' rows (for the x ticks' scan lines)."""
    state: str


def _y(rows: tuple[float, ...]) -> tuple[tuple[float, float], ...]:
    """Pair a side's 21 tick rows with 1e-17, 8, 6, 4, 2e-18, 1e-18, ... 1e-21 cm^2."""
    values = []
    for decade in range(-17, -21, -1):
        values += [10.0**decade] + [m * 10.0 ** (decade - 1) for m in (8, 6, 4, 2)]
    values.append(1e-21)
    return tuple(zip(rows, values, strict=True))


def _x(columns: tuple[float, ...]) -> tuple[tuple[float, float], ...]:
    """Pair a side's six tick columns with 0, 0.2, ..., 1.0 Ry."""
    return tuple(zip(columns, (0.0, 0.2, 0.4, 0.6, 0.8, 1.0), strict=True))


# Fig. 2 (Ar, p. 14) and Fig. 3 (Kr, p. 15): tick centres read from the image, the
# frame lines standing for the end ticks (1e-17 and 1e-21 cm^2; 0 and 1.0 Ry).
FIGURES = {
    "Ar": Frame(
        page=2,
        left=_y(
            (
                2206.0,
                2223.0,
                2247.5,
                2281.5,
                2339.0,
                2396.5,
                2414.0,
                2441.5,
                2475.0,
                2533.0,
                2591.0,
                2609.5,
                2632.5,
                2667.0,
                2723.5,
                2780.0,
                2799.5,
                2823.0,
                2855.0,
                2911.5,
                2970.5,
            )
        ),
        right=_y(
            (
                2207.0,
                2226.0,
                2250.0,
                2284.0,
                2341.5,
                2398.5,
                2418.0,
                2441.5,
                2476.0,
                2533.0,
                2592.0,
                2611.0,
                2633.5,
                2666.0,
                2724.0,
                2781.5,
                2799.0,
                2823.0,
                2857.0,
                2914.0,
                2971.0,
            )
        ),
        top=_x((1412.5, 1509.5, 1606.0, 1702.5, 1799.0, 1897.0)),
        bottom=_x((1409.0, 1507.5, 1603.0, 1700.0, 1796.0, 1893.5)),
        columns=(1411, 1895),
        rows=(2206, 2971),
        state="Ar(4p)",
    ),
    "Kr": Frame(
        page=3,
        left=_y(
            (
                230.5,
                249.5,
                274.0,
                307.0,
                367.5,
                426.0,
                446.5,
                469.0,
                503.0,
                562.5,
                618.5,
                639.0,
                662.5,
                696.0,
                754.5,
                812.5,
                831.5,
                856.5,
                890.0,
                948.0,
                1007.5,
            )
        ),
        right=_y(
            (
                231.0,
                249.5,
                274.5,
                307.5,
                368.0,
                426.0,
                446.5,
                469.0,
                503.0,
                561.5,
                618.5,
                638.5,
                662.5,
                695.5,
                754.5,
                812.0,
                831.0,
                856.0,
                891.0,
                949.5,
                1007.0,
            )
        ),
        top=_x((303.0, 399.0, 495.5, 594.0, 691.5, 788.0)),
        bottom=_x((303.0, 400.5, 493.5, 590.5, 687.0, 786.5)),
        columns=(303, 788),
        rows=(230, 1008),
        state="Kr(5p)",
    ),
}
SCAN_OFFSET = 6
"""How far inside a frame line, pixels, the ticks are looked for."""


def _spans(mask: np.ndarray) -> list[tuple[int, int]]:
    """Return the first and last index of each run of True."""
    index = np.flatnonzero(mask)
    if index.size == 0:
        return []
    groups = np.split(index, np.flatnonzero(np.diff(index) > 1) + 1)
    return [(int(group[0]), int(group[-1])) for group in groups]


def _check_ticks(dark: np.ndarray, frame: Frame) -> None:
    """Raise unless every written tick is a dark run of its side's scan line."""
    left, right = frame.columns
    top, bottom = frame.rows
    sides = {
        "left": (frame.left, dark[:, left + SCAN_OFFSET]),
        "right": (frame.right, dark[:, right - SCAN_OFFSET]),
        "top": (frame.top, dark[top + SCAN_OFFSET, :]),
        "bottom": (frame.bottom, dark[bottom - SCAN_OFFSET, :]),
    }
    for side, (ticks, line) in sides.items():
        centres = [(a + b) / 2 for a, b in _spans(line)]
        for position, value in ticks[1:-1]:
            if min(abs(c - position) for c in centres) > FOUND:
                msg = f"{frame.state}: no {side} tick at {position} ({value})"
                raise ValueError(msg)
    for ticks, log in ((frame.left, True), (frame.right, True)):
        _axis(ticks, log).check()
    for ticks in (frame.top, frame.bottom):
        _axis(ticks, False).check()


def _axis(ticks: tuple[tuple[float, float], ...], log: bool) -> Axis:
    """Return the axis calibrated on one side's ticks."""
    return Axis(tuple(p for p, _ in ticks), tuple(v for _, v in ticks), log, TOLERANCE)


def _value(frame: Frame, column: float, row: float) -> tuple[float, float]:
    """Return (epsilon in Ry, sigma in cm^2) at a pixel, between the frame's sides."""
    left, right = frame.columns
    top, bottom = frame.rows
    across = (column - left) / (right - left)
    down = (row - top) / (bottom - top)
    epsilon = (1 - down) * float(_axis(frame.top, False)(column)) + down * float(
        _axis(frame.bottom, False)(column)
    )
    log_sigma = (1 - across) * math.log10(
        float(_axis(frame.left, True)(row))
    ) + across * math.log10(float(_axis(frame.right, True)(row)))
    return epsilon, 10**log_sigma


def trace(gas: str) -> list[tuple[int, float, float, float]]:
    """Return the p curve of a figure: (column, row, epsilon, sigma) per column.

    Traced from the topmost dark run at ``SEED_EPSILON`` both ways, the run nearest
    the previous position within ``STEP`` pixels; columns where the run is wider than
    ``MAX_WIDTH`` are left out.
    """
    frame = FIGURES[gas]
    dark = image(DIRECTORY, frame.page) < DARK
    _check_ticks(dark, frame)
    left, right = frame.columns
    top, bottom = frame.rows
    inside = (left + 2 * SCAN_OFFSET + 2, right - 2 * SCAN_OFFSET - 2)
    seed = round(left + SEED_EPSILON * (right - left))
    runs = _spans(dark[top + 3 : bottom - 2, seed])
    first = top + 3 + (runs[0][0] + runs[0][1]) / 2
    found: dict[int, tuple[int, int]] = {}
    for direction in (1, -1):
        current, column = first, seed
        while inside[0] <= column <= inside[1]:
            spans = [
                (top + 3 + a, top + 3 + b)
                for a, b in _spans(dark[top + 3 : bottom - 2, column])
            ]
            near = [
                (min(max(current, a), b), (a, b))
                for a, b in spans
                if abs(min(max(current, a), b) - current) <= STEP
            ]
            if near:
                current, span = min(near, key=lambda item: abs(item[0] - current))
                found[column] = span
            column += direction
    out = []
    for column in sorted(found):
        a, b = found[column]
        if b - a + 1 <= MAX_WIDTH:
            row = (a + b) / 2
            epsilon, sigma = _value(frame, column, row)
            out.append((column, row, epsilon, sigma))
    return out


def _at(curve: list[tuple[int, float, float, float]], epsilon: float) -> float:
    """Return the curve's sigma at an excess energy (log-linear interpolation)."""
    eps = np.array([p[2] for p in curve])
    log_sigma = np.log10([p[3] for p in curve])
    order = np.argsort(eps)
    return float(10 ** np.interp(epsilon, eps[order], log_sigma[order]))


def _where(curve: list[tuple[int, float, float, float]], sigma: float) -> float:
    """Return the excess energy at which the (falling) curve has a cross section."""
    eps = np.array([p[2] for p in curve])
    log_sigma = np.log10([p[3] for p in curve])
    order = np.argsort(log_sigma)
    if not log_sigma.min() < math.log10(sigma) < log_sigma.max():
        msg = f"{sigma:.2e} cm2 is not on the curve"
        raise ValueError(msg)
    return float(np.interp(math.log10(sigma), log_sigma[order], eps[order]))


def rows() -> list[tuple[str, str, list[str]]]:
    """Return Table I's rows: the state, its value (cm^2) as text, and the cells."""
    out = []
    for cells in table(DIRECTORY, TABLE_I)[1:]:
        match = CELL.match(cells[1].strip())
        if cells[0] not in STATES or match is None:
            msg = f"Table I: an unexpected row {cells}"
            raise ValueError(msg)
        a, b, power = match.groups()
        out.append((cells[0], f"{a}.{b}e-{power}", cells))
    if [state for state, _, _ in out] != list(STATES):
        msg = "Table I: rows moved"
        raise ValueError(msg)
    return out


def checks() -> dict[str, object]:
    """Return Table I located on the figures, and the 215 nm reading for dakka_2018.

    Raises if a figure's ticks or curve are not where they should be, or if the text's
    repetitions of Table I differ from it.
    """
    values = {state: float(value) for state, value, _ in rows()}
    if values["Ar(4p)"] != 4.3e-18 or values["Kr(5p)"] != 4.5e-18:
        msg = "Table I differs from the text's 4.3e-18 (Ar 4p) and ~4.5e-18"
        raise ValueError(msg)
    photon_ry = HC_EV_NM / (float(WAVELENGTH_A) / 10) / RYDBERG_EV
    out: dict[str, object] = {"photon energy at 2486 A (Ry, 13.6 eV)": photon_ry}
    for gas, frame in FIGURES.items():
        curve = trace(gas)
        epsilon = _where(curve, values[frame.state])
        out[f"{frame.state}: Table I's value on Fig. at epsilon (Ry)"] = epsilon
        out[f"{frame.state}: implied binding energy (eV)"] = (
            photon_ry - epsilon
        ) * RYDBERG_EV
        if gas == "Kr":
            shift = (HC_EV_NM / DAKKA_NM / RYDBERG_EV) - photon_ry
            sigma = _at(curve, epsilon + shift)
            out["Kr(5p) at 215 nm (cm2), read"] = sigma
            out["Kr(5p) at 215 nm: dakka_2018's 3.3e-18 / read"] = DAKKA_VALUE / sigma
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's four values, one file per state."""
    checks()
    files: dict[str, list[Point]] = {}
    for state, value, cells in rows():
        gas, level = STATES[state]
        files[f"{gas}/hyman_1977__table1_{level}.csv"] = [
            Point(WAVELENGTH_A, WAVELENGTH_A, value, "", " | ".join(cells))
        ]
    return files


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hartmann et al. (2022): RMT's xenon ground-state survival at 760 nm.

The paper prints no number of it; Fig. 5(a) draws it ('Xe_sim'). The authors' data
archive (Zenodo 10.5281/zenodo.6683951, the tagged repository ``ATAS_xenon`` v1.1) is
attached to the paper in the library byte for byte, as ``ATAS_xenon-v1.1.zip``; the
script reads it there and checks its sha256 (:data:`ARCHIVE_SHA256`). From it:

* ``data/rmt/gs_pop.csv``: the ground-state population at the end of each RMT run,
  t = 2999.81 au, for seven peak intensities, as the label (1e14 W/cm2) and the value
  printed by RMT (``archer_files/get_pop_GS.sh`` takes the last line of
  ``pop_GS.Xe_ATAS00000200`` in the run with the XUV 5.75 fs after the NIR; the
  README: 'the population agrees to ~6s.f. across the time delays').
* ``archer_files/input.conf``: the RMT input of these runs. Checked: the ground
  state's energy is set to -0.444 au below the ion's (``Adjust_GS_Energy``,
  ``GS_Energy_desired``), the field read from a file, ``Final_T`` 3000 au, ``ml_max``
  0 (linear polarisation along z).
* ``data/rmt/field<I>.csv``: the field RMT wrote (``EField.Xe*``, by
  ``bin/makedf.py``), one column per XUV delay, every 1 au. The NIR is separated from
  the XUV (2.05 au) by a low-pass filter below :data:`CUTOFF` and fitted to E0
  exp(-((t - t0)/w)^2) cos(omega (t - t0) + phi) in every column; the fit must give,
  in every column and at every intensity, E0 = 0.05336 sqrt(I / 1e14 W/cm2) (the
  convention of ``archer_files/gaussian.f90``) within :data:`E0_TOLERANCE`, omega =
  :data:`OMEGA` (760 nm), w = :data:`WIDTH` au and phi = 0 (a cosine carrier) within
  their tolerances. The fit leaves at most 2e-10 au away from the XUV: the NIR is
  analytic, omega = 0.06 au (759.4 nm, the paper's '760 nm'), w = 186 au, t0 = 900 au
  plus the delay. The archived ``gaussian.f90`` is not the generator of these fields:
  it sets omega = 0.06798 au (670 nm) and w = 155 au. The paper's 'Gaussian envelope
  (FWHM 4.5 fs)' is w = 186 au = 4.50 fs, the 1/e half-width of the field; the
  intensity's full width at half maximum is w sqrt(2 ln 2) = 219 au = 5.30 fs.

One data file, ``Xe/hartmann_2022__rmt_ground_survival_Xe.csv``: the row id is the
intensity label as printed, ``x`` the intensity in W/cm2, the value the population as
printed; ``printed`` keeps the file's line.
"""

import csv
import hashlib
import io
import math
import re
import warnings
import zipfile

import numpy as np
from _library import Point, attachment, number, write
from scipy.optimize import OptimizeWarning, curve_fit

PAPER = "Hartmann_2022_CoreResonanceLine"
ARCHIVE = "ATAS_xenon-v1.1.zip"
ARCHIVE_SHA256 = "baac753dcf137b589ac624694e2b544aa4a7210a024aae3b8e2bfa322c7c5ebb"
ROOT = "lhutcheson-ATAS_xenon-90b11b2"
FILE = "Xe/hartmann_2022__rmt_ground_survival_Xe.csv"
INTENSITIES = ("1.1", "1.3", "1.5", "1.7", "1.9", "2.2", "2.5")
# input.conf: the settings this record relies on, as printed.
SETTINGS = {
    "Adjust_GS_Energy": ".true.",
    "GS_Energy_desired": "-0.444",
    "use_field_from_file": ".true.",
    "Final_T": "3000.0",
    "ml_max": "0",
}
FIELD_PER_ROOT_INTENSITY = 0.05336  # au per sqrt(1e14 W/cm2), gaussian.f90
OMEGA = 0.06  # au, 759.4 nm
WIDTH = 186.0  # au, the field's 1/e half-width
CENTRE = 900.0  # au, the NIR's peak in the column without delay
CUTOFF = 0.5  # au, the low-pass filter separating the NIR from the XUV
E0_TOLERANCE = 1e-6  # relative
OMEGA_TOLERANCE = 1e-6  # relative
WIDTH_TOLERANCE = 1e-3  # au
PHASE_TOLERANCE = 1e-4  # rad
END_TIME = "0.29998100000000E+004"  # au, the time get_pop_GS.sh replaces
_LINE = re.compile(r"^\s*(\d\.\d),\s*(\S+)\s*$")


def _archive() -> zipfile.ZipFile:
    """Return the authors' archive, after checking its checksum."""
    path = attachment(PAPER, ARCHIVE)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA256:
        msg = f"{ARCHIVE}: sha256 {digest}, not {ARCHIVE_SHA256}"
        raise ValueError(msg)
    return zipfile.ZipFile(path)


def _member(archive: zipfile.ZipFile, name: str) -> str:
    return archive.read(f"{ROOT}/{name}").decode("utf-8")


def check_input(archive: zipfile.ZipFile) -> None:
    """Raise unless input.conf states :data:`SETTINGS`."""
    text = _member(archive, "archer_files/input.conf")
    found = dict(re.findall(r"^\s*(\w+)\s*=\s*(\S+)", text, flags=re.MULTILINE))
    wrong = {k: found.get(k) for k, v in SETTINGS.items() if found.get(k) != v}
    if wrong:
        msg = f"input.conf: {wrong}, not {SETTINGS}"
        raise ValueError(msg)
    script = _member(archive, "archer_files/get_pop_GS.sh")
    if "delay_+_5.75/data/pop_GS" not in script or END_TIME not in script:
        msg = "get_pop_GS.sh no longer reads the last line of the +5.75 fs run"
        raise ValueError(msg)


def _envelope(
    t: np.ndarray, e0: float, t0: float, w: float, omega: float, phi: float
) -> np.ndarray:
    return e0 * np.exp(-(((t - t0) / w) ** 2)) * np.cos(omega * (t - t0) + phi)


def fit_field(time: np.ndarray, field: np.ndarray) -> tuple[float, ...]:
    """Return the NIR's (E0, t0, w, omega, phi) in one column of a field file."""
    step = time[1] - time[0]
    frequencies = np.fft.rfftfreq(len(field), step) * 2 * math.pi
    spectrum = np.fft.rfft(field)
    spectrum[frequencies > CUTOFF] = 0
    nir = np.fft.irfft(spectrum, len(field))
    peak = int(np.argmax(np.abs(nir)))
    guess = (nir[peak], time[peak], WIDTH, OMEGA, 0.0)
    near = np.abs(time - time[peak]) < 3 * WIDTH
    with warnings.catch_warnings():  # an exact fit leaves no covariance to estimate
        warnings.simplefilter("ignore", OptimizeWarning)
        found, _ = curve_fit(_envelope, time[near], nir[near], p0=guess)
    e0, t0, w, omega, phi = found
    if e0 < 0:
        e0, phi = -e0, phi + math.pi
    phi = (phi + math.pi) % (2 * math.pi) - math.pi
    return e0, t0, abs(w), omega, phi


def check_fields(archive: zipfile.ZipFile) -> dict[str, tuple[float, ...]]:
    """Raise unless every column of every field file is the NIR of the docstring.

    Returns, per intensity, the fit of the column whose NIR peaks at :data:`CENTRE`.
    """
    fits = {}
    for label in INTENSITIES:
        text = _member(archive, f"data/rmt/field{label}.csv")
        rows = list(csv.reader(io.StringIO(text)))
        header, data = rows[0], np.array([[float(c) for c in r] for r in rows[1:]])
        time = data[:, 0]
        expected = FIELD_PER_ROOT_INTENSITY * math.sqrt(float(label))
        for column, name in enumerate(header[1:], start=1):
            e0, t0, w, omega, phi = fit_field(time, data[:, column])
            problems = [
                abs(e0 / expected - 1) > E0_TOLERANCE,
                abs(omega / OMEGA - 1) > OMEGA_TOLERANCE,
                abs(w - WIDTH) > WIDTH_TOLERANCE,
                abs(phi) > PHASE_TOLERANCE,
            ]
            if any(problems):
                msg = (
                    f"field{label}.csv column {name}: E0 {e0:.6f} ({expected:.6f}), "
                    f"omega {omega:.6f}, w {w:.2f}, phi {phi:.4f}"
                )
                raise ValueError(msg)
            if abs(t0 - CENTRE) < 1:
                fits[label] = (e0, t0, w, omega, phi)
        if label not in fits:
            msg = f"field{label}.csv: no column peaks at {CENTRE} au"
            raise ValueError(msg)
    return fits


def survival(archive: zipfile.ZipFile) -> list[tuple[str, str, str]]:
    """Return (label, value, line) for each row of the RMT gs_pop.csv."""
    lines = _member(archive, "data/rmt/gs_pop.csv").splitlines()
    if lines[0].strip() != "Intensity,Pop_GS":
        msg = f"gs_pop.csv: header {lines[0]!r}"
        raise ValueError(msg)
    rows = []
    for line in lines[1:]:
        match = _LINE.match(line)
        if match is None:
            msg = f"gs_pop.csv: unexpected line {line!r}"
            raise ValueError(msg)
        rows.append((match[1], number(match[2]), line.strip()))
    if tuple(label for label, _, _ in rows) != INTENSITIES:
        msg = f"gs_pop.csv: intensities {[r[0] for r in rows]}"
        raise ValueError(msg)
    values = [float(value) for _, value, _ in rows]
    if values != sorted(values, reverse=True) or not 0 < values[-1] < values[0] < 1:
        msg = "gs_pop.csv: the survival does not fall with intensity within (0, 1)"
        raise ValueError(msg)
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the ground-state survival against intensity, after the checks."""
    with _archive() as archive:
        check_input(archive)
        check_fields(archive)
        rows = survival(archive)
    points = [Point(row, f"{row}e14", value, "", line) for row, value, line in rows]
    return {FILE: points}


def main() -> None:
    """Write the data file and print the fitted pulses."""
    with _archive() as archive:
        for label, (e0, t0, w, omega, phi) in check_fields(archive).items():
            fwhm = w * math.sqrt(2 * math.log(2)) * 0.0241888
            print(
                f"{label}e14: E0 {e0:.6f} au, t0 {t0:.2f} au, w {w:.2f} au "
                f"(intensity FWHM {fwhm:.3f} fs), omega {omega:.6f} au "
                f"({45.5634 / omega:.1f} nm), phi {phi:+.4f}"
            )
    write(extract())


if __name__ == "__main__":
    main()

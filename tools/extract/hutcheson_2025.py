# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Check Hutcheson et al. (2025)'s record against the authors' data: the NIR's carrier.

No value is taken from this source (its record says why). The script checks what the
record states about the RMT populations of Fig. 4, from the authors' archive (Zenodo
10.5281/zenodo.15432565, the tagged repository ``dipole_response`` 1.0), attached to the
paper in the library byte for byte: ``RMTpopulations.csv`` (time, and the ground-state,
selected bound-state and r > 20 au populations at eight NIR intensities) and the whole
archive ``dipole_response-1.0.zip``.

* The populations oscillate at twice the NIR's carrier frequency (the paper, p. 2: 'the
  ground-state behavior consists of ... oscillations at twice the driving NIR
  frequency'). The peak of the spectrum of d(ground)/dt between :data:`BAND` lies within
  :data:`AGREEMENT` of 2 x 0.06798 au (670 nm) at every intensity, and at least
  :data:`APART` from 2 x 0.05995 au (760 nm, the paper's wavelength).
* The archive's ``bin/plot_rmt_populations.py`` draws these populations with
  ``pulse(0.8, IR_freq=0.06798)``, and ``bin/helper_functions.py``'s ``pulse`` defaults
  to ``FWHM=186`` (au, the field's 1/e half-width) and E0 = 0.05336 sqrt(I).
"""

import csv
import math
import zipfile

import numpy as np
from _library import Point, attachment

PAPER = "Hutcheson_2025_PhaseEvolutionStrong"
POPULATIONS = "RMTpopulations.csv"
ARCHIVE = "dipole_response-1.0.zip"
ROOT = "lhutcheson-dipole_response-9d46dd1"
INTENSITIES = ("0.8", "1.1", "1.3", "1.5", "1.7", "1.9", "2.2", "2.5")
PLOTTED = 0.06798  # au, the carrier of plot_rmt_populations.py
STATED = 45.5634 / 760  # au, the paper's 760 nm
WINDOW = (400.0, 1400.0)  # au, the NIR (peak at 900 au)
BAND = (0.09, 0.2)  # au
AGREEMENT = 0.02  # relative, of 2 omega
APART = 0.08  # relative, of 2 omega
PADDING = 8
SNIPPETS = {
    "bin/plot_rmt_populations.py": ("pulse(0.8, IR_freq=0.06798)",),
    "bin/helper_functions.py": ("FWHM=186", "E0 = 0.05336*np.sqrt(intensity)"),
}


def _populations() -> tuple[list[str], np.ndarray]:
    with attachment(PAPER, POPULATIONS).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    header = rows[0][1:]
    data = np.array([[float(cell) for cell in row[1:]] for row in rows[1:]])
    return header, data


def carriers() -> dict[str, float]:
    """Return omega (au), half the 2-omega peak of d(ground)/dt, per intensity."""
    header, data = _populations()
    time = data[:, header.index("Time")]
    step = time[1] - time[0]
    near = (time > WINDOW[0]) & (time < WINDOW[1])
    out = {}
    for label in INTENSITIES:
        ground = data[:, header.index(f"Ground {label}")][near]
        rate = np.gradient(ground, step) * np.hanning(near.sum())
        frequencies = np.fft.rfftfreq(PADDING * len(rate), step) * 2 * math.pi
        spectrum = np.abs(np.fft.rfft(rate, PADDING * len(rate)))
        band = (frequencies > BAND[0]) & (frequencies < BAND[1])
        out[label] = float(frequencies[band][np.argmax(spectrum[band])]) / 2
    return out


def check_archive() -> None:
    """Raise unless the archive's plotting code states the pulse of the docstring."""
    with zipfile.ZipFile(attachment(PAPER, ARCHIVE)) as archive:
        for name, snippets in SNIPPETS.items():
            text = archive.read(f"{ROOT}/{name}").decode("utf-8")
            missing = [s for s in snippets if s not in text]
            if missing:
                msg = f"{name}: {missing} not found"
                raise ValueError(msg)


def extract() -> dict[str, list[Point]]:
    """Check the record's statements; there are no data files to write."""
    check_archive()
    for label, omega in carriers().items():
        if abs(omega / PLOTTED - 1) > AGREEMENT or abs(omega / STATED - 1) < APART:
            msg = f"{label}e14: the populations oscillate at 2 x {omega:.4f} au"
            raise ValueError(msg)
    return {}


def main() -> None:
    """Run the checks and print the carriers and the final ground-state populations."""
    extract()
    header, data = _populations()
    for label, omega in carriers().items():
        ground = data[-1, header.index(f"Ground {label}")]
        print(
            f"{label}e14: omega {omega:.4f} au ({45.5634 / omega:.0f} nm), "
            f"ground {ground:.5f} at {data[-1, header.index('Time')]:.2f} au"
        )


if __name__ == "__main__":
    main()

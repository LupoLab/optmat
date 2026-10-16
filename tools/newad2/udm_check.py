# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Check optmat's universal dispersion model against the authors' own program.

Franta et al. publish ``udm.cxx`` with their newAD2 parameter files
(http://newad.physics.muni.cz/). It is not committed here: it bundles GPL-3 code
(TC-320).
Download it, then::

    NEWAD2_UDM_SOURCE=/path/to/udm.cxx uv run python tools/newad2/udm_check.py

The script compiles a copy that prints 17 significant digits (the program prints 7),
runs it on every newAD2 parameter file attached to Franta et al.'s papers in the paper
library (``$PAPER_LIBRARY``) and reports the largest
relative difference in epsilon above and below 1 eV. Above 1 eV the two agree to about
1e-10. Below 1 eV the program approximates the higher-energy term, and on a machine
whose ``long double`` is only double precision (Apple silicon) its closed forms lose
digits; there optmat's implementation is the exact one (TC-321). A C++ compiler and
Boost's headers are needed.
"""

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from optmat.solid.udm import from_file as udm_from_file

# The papers the newAD2 files are attached to, in the paper library.
PAPERS = (
    "Franta_2016_OpticalCharacterizationSiO",
    "Franta_2015_DispersionModelOptical",
    "Franta_2017_TemperatureDependentDispersion",
)
LIBRARY = Path(os.environ.get("PAPER_LIBRARY", "PAPER_LIBRARY-is-not-set"))


def _compile(source: Path, work: Path) -> Path:
    text = source.read_text(encoding="utf-8")
    precise, count = re.subn(
        r'printf\("%e\\t%e\\t%e\\n",val,real\(oc\),imag\(oc\)\);',
        'printf("%.17e\\\\t%.17e\\\\t%.17e\\\\n",val,real(oc),imag(oc));',
        text,
    )
    if count != 1:
        sys.exit("udm.cxx: the output statement was not found; is this the right file?")
    patched = work / "udm.cxx"
    patched.write_text(precise, encoding="utf-8")
    binary = work / "udm"
    include = os.environ.get("BOOST_INCLUDE", "/opt/homebrew/include")
    compiler = os.environ.get("CXX", "c++")
    command = [
        compiler,
        "-O2",
        "-std=c++17",
        f"-I{include}",
        "-o",
        str(binary),
        str(patched),
    ]
    subprocess.run(command, check=True)
    return binary


def main() -> None:
    """Compile the program and compare it with every parameter file."""
    source = os.environ.get("NEWAD2_UDM_SOURCE")
    if not source:
        sys.exit("set NEWAD2_UDM_SOURCE to the authors' udm.cxx")
    with tempfile.TemporaryDirectory() as work:
        binary = _compile(Path(source), Path(work))
        files = sorted(p for d in PAPERS for p in (LIBRARY / d).glob("data/*/*.par"))
        for par in files:
            output = subprocess.run(
                [str(binary), str(par)], check=True, capture_output=True, text=True
            ).stdout
            table = np.loadtxt(output.splitlines())
            energy = table[:, 0]
            theirs = (table[:, 1] + 1j * table[:, 2]) ** 2
            ours = udm_from_file(par).permittivity(energy)
            finite = np.isfinite(theirs)
            relative = np.abs(ours - theirs) / np.abs(theirs)
            high = finite & (energy > 1)
            low = finite & (energy <= 1)
            if not finite.any():
                print(f"{par.stem:24s} the program returns NaN on this machine")
                continue
            above, below = relative[high].max(), relative[low].max()
            print(
                f"{par.stem:24s} max relative difference in epsilon: "
                f"above 1 eV {above:.1e}, below {below:.1e}"
            )


if __name__ == "__main__":
    main()

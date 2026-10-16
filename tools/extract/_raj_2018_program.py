# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

# /// script
# requires-python = "==3.11.*"
# dependencies = ["numpy==1.24.4", "scipy==1.10.1"]
# ///
r"""Run Raj, Hamaguchi & Witek's ``rovibME`` Python module, unmodified, on a job list.

    uv run --script tools/extract/_raj_2018_program.py < jobs.json > results.jsonl

Called by ``raj_2018.py --program``; not an extraction script (the leading underscore
keeps it out of the library test). The authors' module (``rovibME.py`` in the
``python-module`` directory of their H2-PolarizabilityMatrixElements package, J. Chem.
Phys. 148, 104308 (2018), Ref. 39 and supplementary material) is imported from where
it lies and never copied or edited. It needs an environment of its time: it passes
``(n, 1)`` arrays where NumPy 2 requires scalars, and calls
``scipy.integrate.quadrature``, removed in SciPy 1.15. The inline metadata above pins
Python 3.11, NumPy 1.24.4 and SciPy 1.10.1, for which the module runs without a warning
(the one exception, an invalid ``\g`` escape in a docstring, is a compile-time
DeprecationWarning of Python 3.11, ignored at import).

``rovibME.compute`` returns nothing: it prints ``<mol> < v=.. J=.. | isotropy | v=..
J=.. >  =  <value> a.u. (Integration err: <err>)`` with the value rounded to seven
decimals. Each job's printed line is captured and parsed; the line is passed back
whole, for the data files' ``printed`` column.

Input (stdin, JSON): ``{"module": <python-module directory>, "jobs": [[molecule,
bra_v, bra_J, ket_v, ket_J, wavelength, unit], ...]}`` with ``wavelength`` a number or
``"static"`` and ``unit`` one of the module's ("nm", "A", "H"). Output (stdout): one
JSON header line with the versions and the module's SHA-256, then one line per job:
``[job, value, integration_error, printed_line]``.
"""

import contextlib
import hashlib
import io
import json
import os
import platform
import re
import sys
import warnings
from pathlib import Path

LINE = re.compile(
    r"^(?P<line>\S+ < v=\d+ J=\d+ \| isotropy \| v=\d+ J=\d+ >\s+=\s+"
    r"(?P<value>[-+0-9.eE]+) a\.u\. \(Integration err: (?P<err>[-+0-9.eE]+)\))\s*$",
    re.MULTILINE,
)


def main() -> None:
    """Run every job and print the results as JSON lines."""
    request = json.load(sys.stdin)
    module = Path(request["module"])
    sys.dont_write_bytecode = True  # never write a __pycache__ into the authors' tree
    os.chdir(module)  # the module reads ./data and ./wavefunctions relative to the cwd
    sys.path.insert(0, str(module))
    with (
        contextlib.redirect_stdout(io.StringIO()),
        warnings.catch_warnings(),
    ):
        warnings.simplefilter("ignore", DeprecationWarning)  # the docstring escape
        import numpy  # noqa: PLC0415
        import rovibME  # pyright: ignore[reportMissingImports]  # noqa: PLC0415
        import scipy  # noqa: PLC0415
    source = (module / "rovibME.py").read_bytes()
    print(
        json.dumps(
            {
                "python": platform.python_version(),
                "numpy": numpy.__version__,
                "scipy": scipy.__version__,
                "rovibME_sha256": hashlib.sha256(source).hexdigest(),
            }
        )
    )
    for job in request["jobs"]:
        molecule, bra_v, bra_j, ket_v, ket_j, wavelength, unit = job
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer), warnings.catch_warnings():
            # A warning would mean a changed numeric path: make it an error.
            warnings.simplefilter("error")
            rovibME.compute(
                molecule, bra_v, bra_j, ket_v, ket_j, wavelength, unit, "mp"
            )
        found = LINE.findall(buffer.getvalue())
        if len(found) != 1:
            msg = f"unexpected output for {job}: {buffer.getvalue()!r}"
            raise SystemExit(msg)
        line, value, err = found[0]
        print(json.dumps([job, value, err, line]))


if __name__ == "__main__":
    main()

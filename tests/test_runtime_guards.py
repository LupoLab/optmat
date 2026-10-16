# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Public validation also works in ordinary Python, without the test import hook."""

import dataclasses
import subprocess
import sys
from pathlib import Path
from typing import cast

import numpy as np
import pytest


def _materials():
    from optmat.exceptions import OptmatError  # noqa: PLC0415
    from optmat.materials import materials  # noqa: PLC0415

    argon = materials()["Ar"]
    assert dataclasses.replace(argon, name="synthetic argon").kind == "gas"
    with pytest.raises(OptmatError, match="unknown kind 'missing'"):
        dataclasses.replace(argon, kind="missing")
    with pytest.raises(OptmatError, match="unknown optics 'missing'"):
        dataclasses.replace(argon, optics="missing")


def _distributions():
    from optmat.gas.euv import OscillatorDistribution  # noqa: PLC0415
    from optmat.gas.gases import Gas  # noqa: PLC0415

    valid = OscillatorDistribution(
        Gas.He,
        "synthetic",
        np.array([1.0]),
        np.array([0.5]),
        np.array([2.0, 3.0]),
        np.array([0.3, 0.3]),
    )
    assert valid.threshold == 2.0
    # Without validation, NumPy can silently broadcast one line's strength over two
    # energies. Check both independent table pairs, not just a failing calculation.
    with pytest.raises(ValueError, match="tables of unequal length"):
        dataclasses.replace(valid, line_energies=np.array([1.0, 1.5]))
    with pytest.raises(ValueError, match="tables of unequal length"):
        dataclasses.replace(valid, continuum_density=np.array([0.3]))


def _coulomb():
    from optmat.gas.ionisation_rates import (  # noqa: PLC0415
        CoulombFactor,
        coulomb_factor,
    )

    # kappa=1 and n*=1: the static factor is (2/F)^2 = 1600.
    assert coulomb_factor(np.array(0.05), 0.05, 0.5) == pytest.approx(1600, rel=2e-15)
    with pytest.raises(ValueError, match="unknown Coulomb factor 'missing'"):
        coulomb_factor(np.array(0.05), 0.05, 0.5, cast(CoulombFactor, "missing"))


def _raman_components():
    from optmat.exceptions import DataError  # noqa: PLC0415
    from optmat.gas.raman_response import Component, raman_response  # noqa: PLC0415
    from optmat.solid.raman_response import (  # noqa: PLC0415
        Component as SolidComponent,
    )
    from optmat.solid.raman_response import (  # noqa: PLC0415
        silica_raman_response,
    )

    gas = raman_response("H2", 300.0)
    np.testing.assert_array_equal(gas.component_amplitudes("parallel"), gas.amplitude)
    with pytest.raises(DataError, match="unknown component 'missing'"):
        gas.component_amplitudes(cast(Component, "missing"))
    solid = silica_raman_response("lin_2006")
    assert solid.response(0.0, "parallel") == 0.0
    with pytest.raises(DataError, match="unknown component 'missing'"):
        solid.response(0.0, cast(SolidComponent, "missing"))


CHECKS = {
    "materials": _materials,
    "distributions": _distributions,
    "coulomb": _coulomb,
    "raman_components": _raman_components,
}


@pytest.mark.parametrize("case", CHECKS)
def test_public_guards_in_normal_python(case, tmp_path):
    # Run this file as a normal script: the suite itself keeps its jaxtyping hook.
    # -W error preserves the suite's warning policy in the child process.
    result = subprocess.run(
        [sys.executable, "-W", "error", str(Path(__file__).resolve()), case],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


if __name__ == "__main__":
    CHECKS[sys.argv[1]]()

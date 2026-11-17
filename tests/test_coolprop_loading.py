# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""CoolProp loaded without superancillaries: the same values, quietly (TC-1096).

optmat loads CoolProp without the superancillaries it would build for every fluid. The
densities and saturation pressures optmat reads must agree with those of a CoolProp
loaded with them, over every gas's phase domain from its triple point to 600 K and 200
bar, and CoolProp's notice of the switch must not reach standard output. Each
loading needs a process of its own: CoolProp reads the switch once, when it loads.
"""

import json
import os
import subprocess
import sys

import pytest

from optmat.gas import eos

_SWITCH = "COOLPROP_DISABLE_SUPERANCILLARIES_ENTIRELY"

_SAMPLE = """
import json
import numpy as np
from optmat.conditions import State
from optmat.gas.eos import equation_of_state
from optmat.gas.gases import Gas, info

values = {}
for gas in (g for g in Gas if info(g).fluid is not None):
    e = equation_of_state(gas)
    low, high = e.temperature_range
    for t in np.geomspace(low * (1 + 1e-9), min(high, 600.0), 12):
        bottom, top = e.phase_domain(float(t))
        # The same states in both runs: pressures on a fixed grid, clear of the top,
        # which moves with the saturation pressure.
        for p in np.geomspace(1e2, 2e7, 40):
            if bottom <= p < top * (1 - 1e-6):
                state = State(float(p), float(t))
                values[f"{gas.value} {state!r}"] = e.number_density(state)
        saturation = e.saturation_pressure(float(t))
        if saturation is not None:
            values[f"{gas.value} {float(t)!r} K saturation"] = saturation
print(json.dumps(values))
"""


def _run(code: str, *, superancillaries: bool) -> str:
    """Return the standard output of code run in a fresh interpreter.

    - ``superancillaries``: whether CoolProp loads them; it is imported first, before
      optmat would load it without them, to keep them.
    """
    environment = {k: v for k, v in os.environ.items() if k != _SWITCH}
    first = "import CoolProp.CoolProp\n" if superancillaries else ""
    result = subprocess.run(
        [sys.executable, "-c", first + code],
        capture_output=True,
        text=True,
        check=True,
        env=environment,
    )
    return result.stdout


def test_the_values_read_agree_with_and_without_superancillaries():
    without = json.loads(_run(_SAMPLE, superancillaries=False))
    with_them = json.loads(_run(_SAMPLE, superancillaries=True))
    assert without.keys() == with_them.keys()
    assert sum(key.endswith("saturation") for key in without) > 50
    for key, value in with_them.items():
        # Measured 2026-11-16: every density bitwise the same; the saturation pressures
        # within 3e-11, but oxygen's at its triple point (146 Pa) within 2.0e-9, the
        # tolerance of CoolProp's saturation solve, which they then come from.
        bound = 5e-9 if key.endswith("saturation") else 1e-12
        assert without[key] == pytest.approx(value, rel=bound, abs=0.0), key


_PROBE = """
from optmat.gas.eos import _abstract_state
state = _abstract_state("Helium")
try:
    state.update_QT_pure_superanc(1.0, 4.0)
except ValueError:
    print("without")
else:
    print("with")
"""


def test_coolprop_is_loaded_without_superancillaries_and_says_nothing():
    # Standard output holds the probe's word alone: CoolProp's notice is not there.
    assert _run(_PROBE, superancillaries=False).split() == ["without"]


def test_child_processes_do_not_inherit_the_switch():
    # A child would take an inherited switch for the user's and print CoolProp's notice.
    code = _PROBE + "import os; print(os.environ.get('" + _SWITCH + "'))\n"
    assert _run(code, superancillaries=False).split() == ["without", "None"]


_EARLIER_OUTPUT = """
import ctypes
from optmat.gas import eos
libc = ctypes.CDLL(None)
libc.printf(b"before\\n")
with eos._quiet_stdout():
    libc.printf(b"hidden\\n")
libc.printf(b"after\\n")
libc.fflush(None)
"""


@pytest.mark.skipif(os.name != "posix", reason="libc's printf is reached through POSIX")
def test_native_output_written_before_the_block_is_kept():
    # Buffered in C (standard output is a pipe) until the block flushes it first.
    assert _run(_EARLIER_OUTPUT, superancillaries=False).split() == ["before", "after"]


def test_a_coolprop_already_loaded_keeps_its_superancillaries():
    assert _run(_PROBE, superancillaries=True).split() == ["with"]


def test_a_switch_set_by_the_user_is_kept_and_coolprop_says_so():
    environment = {**os.environ, _SWITCH: "yes"}
    code = _PROBE + "import os; print(os.environ['" + _SWITCH + "'])\n"
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
        env=environment,
    )
    # The user's value stands, and CoolProp's own notice of it is not hidden.
    assert result.stdout.split()[:2] == ["without", "yes"]
    assert "superancillaries have been disabled" in result.stdout


def test_without_a_standard_output_descriptor_the_block_still_runs(monkeypatch):
    def refused(fd: int) -> int:
        raise OSError(fd)

    monkeypatch.setattr(eos.os, "dup", refused)
    ran = []
    with eos._quiet_stdout():
        ran.append(True)
    assert ran == [True]

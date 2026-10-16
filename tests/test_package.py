# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import pkgutil
import subprocess
import sys

import pytest

import optmat


def test_submodule_list_matches_the_modules_on_disk():
    on_disk = {
        module.name
        for module in pkgutil.iter_modules(optmat.__path__)
        if not module.name.startswith("_")
    }
    assert set(optmat._SUBMODULES) == on_disk


def test_submodules_load_lazily_and_are_listed():
    assert optmat.gas.gases.Gas.Ar == "Ar"
    assert "conditions" in dir(optmat)


def test_unknown_attribute_raises():
    with pytest.raises(AttributeError, match="no attribute"):
        _ = optmat.no_such_module


def test_import_does_not_load_coolprop():
    code = "import sys, optmat; print('CoolProp' in sys.modules)"
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "False"


def test_the_propagation_surface_defers_the_slow_scipy_modules():
    # Importing the modules a propagator uses loads neither scipy.signal (which brings
    # scipy.stats) nor scipy.ndimage, scipy.interpolate or scipy.optimize; O3's measured
    # absorption loads them on first use (TC-810).
    code = """
import sys
import numpy as np
import optmat.gas.air, optmat.gas.eos, optmat.gas.ionisation_rates
import optmat.gas.models.registry, optmat.gas.molecular_ionisation_rates
import optmat.gas.nonlinear, optmat.gas.raman_response, optmat.gas.transport
import optmat.solid.bundle, optmat.solid.models.registry
slow = ["scipy.signal", "scipy.stats", "scipy.ndimage", "scipy.interpolate",
        "scipy.optimize"]
print([m for m in slow if m in sys.modules])
alpha = optmat.gas.models.registry.model("fit/O3").complex_polarizability(255e-9)
print(bool(np.isfinite(alpha) and alpha.imag > 0))
print([m for m in slow if m not in sys.modules])
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.split("\n")[:3] == ["[]", "True", "[]"]


def test_the_core_imports_no_medium():
    # The core (records, vocabulary, materials, states, constants) serves every medium
    # and imports none of their subpackages (TC-317, TC-318).
    import ast  # noqa: PLC0415
    from pathlib import Path  # noqa: PLC0415

    package = Path(__file__).parents[1] / "src" / "optmat"
    media = {path.name for path in package.iterdir() if (path / "__init__.py").exists()}
    assert "gas" in media
    for path in package.glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom) and node.module:
                names = {node.module}
            elif isinstance(node, ast.Import):
                names = {alias.name for alias in node.names}
            else:
                continue
            for name in names:
                parts = name.split(".")
                assert not (
                    parts[0] == "optmat" and len(parts) > 1 and parts[1] in media
                ), f"{path.name} imports {name}"


@pytest.mark.parametrize("subpackage", ["gas", "solid"])
def test_a_subpackages_list_matches_its_modules_on_disk(subpackage):
    package = getattr(optmat, subpackage)
    on_disk = {
        module.name
        for module in pkgutil.iter_modules(package.__path__)
        if not module.name.startswith("_")
    }
    assert set(package._SUBMODULES) == on_disk

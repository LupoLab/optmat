# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The public surface is declared, and its changes are deliberate (OEP 1, §6.1).

``tests/data/public_api.json`` is written by ``tools/public_api.py``. When the surface
changes on purpose, re-run the tool and add a changelog line in the same commit.
"""

import difflib
import importlib
import importlib.util
import inspect
import json
import sys
from pathlib import Path

import pytest

# The recorded signatures are the annotations as written, which only Python 3.14's
# annotationlib returns (the tool imports it); the surface does not depend on the
# version.
if sys.version_info < (3, 14):
    pytest.skip(
        "the snapshot is Python 3.14's annotation text", allow_module_level=True
    )

ROOT = Path(__file__).parents[1]
RECORDED = (ROOT / "tests" / "data" / "public_api.json").read_text(encoding="utf-8")


def _tool():
    path = ROOT / "tools" / "public_api.py"
    spec = importlib.util.spec_from_file_location("public_api", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("public_api", module)
    spec.loader.exec_module(module)
    return module


def test_the_surface_matches_its_snapshot():
    now = json.dumps(_tool().snapshot(), indent=1, sort_keys=True) + "\n"
    if now != RECORDED:
        diff = "".join(
            difflib.unified_diff(
                RECORDED.splitlines(keepends=True),
                now.splitlines(keepends=True),
                "recorded",
                "now",
                n=2,
            )
        )
        pytest.fail(
            "the public surface changed; if that is intended, re-run "
            f"tools/public_api.py and add a changelog line:\n{diff}"
        )


@pytest.mark.parametrize("name", _tool().SURFACE_MODULES)
def test_every_public_definition_is_declared(name):
    module = importlib.import_module(name)
    defined = {
        n
        for n, value in vars(module).items()
        if not n.startswith("_")
        and (inspect.isclass(value) or inspect.isfunction(value))
        and value.__module__ == name
    }
    assert defined <= set(module.__all__), sorted(defined - set(module.__all__))
    assert all(hasattr(module, n) for n in module.__all__)

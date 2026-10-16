# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The ionisation page's tables are what the curated parameters give today.

``docs/ionisation.md`` holds three tables written by ``tools/ionisation_tables.py``: the
parameters of each gas, the tunnelling rate against the ab initio rates, and solprop's
values against the curated ones. When a parameter or its source changes, run the script
and paste its output over the tables (and update the page's text where a finding
changes).
"""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "ionisation.md").read_text()


def _tables():
    path = ROOT / "tools" / "ionisation_tables.py"
    spec = importlib.util.spec_from_file_location("ionisation_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("ionisation_tables", module)
    spec.loader.exec_module(module)
    return module


def _block(name: str) -> str:
    text = (
        (ROOT / "tests/data/historical_documentation_tables.md").read_text()
        if name in ("ionisation-solprop",)
        else PAGE
    )
    return text.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


@pytest.mark.parametrize(
    ("name", "function"),
    [
        ("ionisation-parameters", "parameter_table"),
        ("ionisation-validation", "validation_table"),
        ("ionisation-solprop", "solprop_table"),
    ],
)
def test_the_page_table_is_current(name, function):
    assert _block(name) == getattr(_tables(), function)()

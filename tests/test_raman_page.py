# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The Raman page's tables are what the curated parameters give today.

``docs/raman.md`` holds three tables written by ``tools/raman_tables.py``. When a
parameter or its source changes, run the script and paste its output over the tables
(and update the page's text where a finding changes).
"""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "raman.md").read_text()


def _tables():
    path = ROOT / "tools" / "raman_tables.py"
    spec = importlib.util.spec_from_file_location("raman_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("raman_tables", module)
    spec.loader.exec_module(module)
    return module


def _block(name: str) -> str:
    location = (
        "tests/data/historical_documentation_tables.md"
        if name == "raman-solprop"
        else "docs/user-guide/raman-validation.md"
    )
    page = (ROOT / location).read_text()
    # Some response tables describe the interface and stay with its usage chapter.
    if f"<!-- {name}:start -->" not in page:
        page = PAGE
    return page.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


@pytest.mark.parametrize(
    ("name", "function"),
    [
        ("raman-parameters", "parameter_table"),
        ("raman-gain", "gain_table"),
        ("raman-rotational", "rotational_table"),
        ("raman-solprop", "solprop_table"),
    ],
)
def test_the_page_table_is_current(name, function):
    assert _block(name) == getattr(_tables(), function)()

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The molecular-rates page's tables are what the models give today.

``docs/molecular_ionisation_rates.md`` holds the tables written by
``tools/molecular_rate_tables.py``. When a model or a reference changes, run the
script and paste its output over the tables (and update the page's text where a
finding changes).
"""

import importlib.util
import sys
import warnings
from pathlib import Path

import pytest

from optmat.exceptions import ModelRangeWarning

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "molecular_ionisation_rates.md").read_text(encoding="utf-8")
PAGE += (ROOT / "docs" / "user-guide" / "molecular-ionisation-validation.md").read_text(
    encoding="utf-8"
)
NAMES = (
    "molecular-rate-constants",
    "molecular-rate-models",
    "molecular-rate-forms",
    "molecular-rate-static",
    "molecular-rate-vibration",
    "molecular-rate-wahlstrand",
    "molecular-rate-solprop",
)


def _tables():
    path = ROOT / "tools" / "molecular_rate_tables.py"
    spec = importlib.util.spec_from_file_location("molecular_rate_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("molecular_rate_tables", module)
    spec.loader.exec_module(module)
    return module


def _block(name: str) -> str:
    text = (
        (ROOT / "tests/data/historical_documentation_tables.md").read_text(
            encoding="utf-8"
        )
        if name in ("molecular-rate-solprop",)
        else PAGE
    )
    return text.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


def test_the_script_writes_every_table_of_the_page():
    assert tuple(_tables().TABLES) == NAMES


@pytest.mark.parametrize("name", NAMES)
def test_the_page_table_is_current(name):
    with warnings.catch_warnings():
        # The forms are shown above the barrier-suppression field, where they have
        # no correction, as the script's main() does.
        warnings.simplefilter("ignore", ModelRangeWarning)
        assert _block(name) == _tables().TABLES[name]()

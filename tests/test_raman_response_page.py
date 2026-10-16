# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The page *The Raman response*'s tables are what the responses give today.

``docs/raman_response.md`` holds tables written by ``tools/raman_response_tables.py``.
When a response, a parameter or its source changes, run the script and paste its
output over the tables (and update the page's text where a finding changes).
"""

import importlib.util
import sys
import warnings
from pathlib import Path

import pytest

from optmat.exceptions import ModelRangeWarning

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "raman_response.md").read_text()


def _tables():
    path = ROOT / "tools" / "raman_response_tables.py"
    spec = importlib.util.spec_from_file_location("raman_response_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("raman_response_tables", module)
    spec.loader.exec_module(module)
    return module


def _block(name: str) -> str:
    location = (
        "tests/data/historical_documentation_tables.md"
        if name == "response-solprop"
        else "docs/user-guide/raman-validation.md"
    )
    page = (ROOT / location).read_text()
    # Some response tables describe the interface and stay with its usage chapter.
    if f"<!-- {name}:start -->" not in page:
        page = PAGE
    return page.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


@pytest.mark.parametrize("name", list(_tables().TABLES))
def test_the_page_table_is_current(name):
    with warnings.catch_warnings():
        # The validators' and the rank study's conditions lie outside some ranges.
        warnings.simplefilter("ignore", ModelRangeWarning)
        table = _tables().TABLES[name]()
    assert _block(name) == table

# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The ionisation-rates page's tables are what the models give today.

``docs/ionisation_rates.md`` holds the tables written by
``tools/ionisation_rate_tables.py``. When a model, a constant or a reference changes,
run the script and paste its output over the tables (and update the page's text where
a finding changes).
"""

import importlib.util
import sys
import warnings
from pathlib import Path

import pytest

from optmat.exceptions import ModelRangeWarning

ROOT = Path(__file__).parents[1]
PAGE = (ROOT / "docs" / "ionisation_rates.md").read_text(encoding="utf-8")
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-theory.md").read_text(
    encoding="utf-8"
)
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-static.md").read_text(
    encoding="utf-8"
)
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-laser.md").read_text(
    encoding="utf-8"
)
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-one-electron.md").read_text(
    encoding="utf-8"
)
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-measured.md").read_text(
    encoding="utf-8"
)
PAGE += (ROOT / "docs" / "user-guide" / "ionisation-accuracy.md").read_text(
    encoding="utf-8"
)
NAMES = (
    "ionisation-rate-models",
    "ionisation-rate-electrons",
    "ionisation-rate-weak-field",
    "ionisation-rate-atoms",
    "ionisation-rate-static",
    "ionisation-rate-helium",
    "ionisation-rate-helium-two-electron",
    "ionisation-rate-helium-matched",
    "ionisation-rate-xenon-rmt",
    "ionisation-rate-xenon-measured",
    "ionisation-rate-rmatrix",
    "ionisation-rate-tdse",
    "ionisation-rate-hydrogen",
    "ionisation-rate-tdse42",
    "ionisation-rate-saturation",
    "ionisation-rate-wahlstrand",
    "ionisation-rate-wallace",
    "ionisation-rate-wallace-points",
    "ionisation-rate-accuracy",
    "ionisation-rate-legacy",
    "ionisation-rate-legacy-tests",
    "ionisation-rate-phase",
    "ionisation-rate-luna",
    "ionisation-rate-luna-tests",
)


def _tables():
    path = ROOT / "tools" / "ionisation_rate_tables.py"
    spec = importlib.util.spec_from_file_location("ionisation_rate_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("ionisation_rate_tables", module)
    spec.loader.exec_module(module)
    return module


def _block(name: str) -> str:
    text = (
        (ROOT / "tests/data/historical_documentation_tables.md").read_text(
            encoding="utf-8"
        )
        if name in ("ionisation-rate-legacy", "ionisation-rate-legacy-tests")
        else PAGE
    )
    return text.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


def test_the_script_writes_every_table_of_the_page():
    assert tuple(_tables().TABLES) == NAMES


@pytest.mark.parametrize("name", NAMES)
def test_the_page_table_is_current(name):
    with warnings.catch_warnings():
        # The literature's forms are compared beyond the range of their missing
        # static correction, as the script's main() does.
        warnings.simplefilter("ignore", ModelRangeWarning)
        assert _block(name) == _tables().TABLES[name]()

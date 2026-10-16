# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The Python snippets of the README and the manual's pages run, page by page.

A page's ``python`` fenced blocks share one namespace, as a reader would run them.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
PAGES = (
    "README.md",
    "docs/introduction/installation.md",
    "docs/quickstart.md",
    "docs/identity.md",
    "docs/models.md",
    "docs/conventions.md",
    "docs/transport.md",
    "docs/mixtures.md",
    "docs/propagation.md",
    "docs/nonlinear.md",
    "docs/user-guide/nonlinear-support.md",
    "docs/developers/fitting.md",
    "docs/fitting.md",
    "docs/raman.md",
    "docs/raman_response.md",
    "docs/ionisation.md",
    "docs/ionisation_rates.md",
    "docs/user-guide/ionisation-theory.md",
    "docs/two_photon.md",
    "docs/user-guide/two-photon-excitation.md",
    "docs/molecular_ionisation_rates.md",
    "docs/silica.md",
    "docs/crystals.md",
    "docs/metals.md",
)
BLOCK = re.compile(r"^```python\n(.*?)^```", re.MULTILINE | re.DOTALL)


@pytest.mark.parametrize("page", PAGES)
def test_the_page_snippets_run(page):
    namespace: dict[str, object] = {}
    blocks = BLOCK.findall((ROOT / page).read_text(encoding="utf-8"))
    assert blocks, page
    for number, block in enumerate(blocks, 1):
        exec(compile(block, f"{page}, block {number}", "exec"), namespace)

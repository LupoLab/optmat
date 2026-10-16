# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Test session configuration.

Shape and dtype annotations are checked at runtime during tests only: the jaxtyping
import hook instruments every ``optmat`` module imported after it with
beartype, configured with the PEP 484 numeric tower as in solprop.
"""

from jaxtyping import install_import_hook

JAXTYPING_HOOK = install_import_hook(
    "optmat",
    "beartype.beartype(conf=beartype.BeartypeConf(is_pep484_tower=True))",
)

import optmat  # noqa: E402, F401 - imported after the hook so it is checked

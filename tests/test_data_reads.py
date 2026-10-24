# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""How many data files the entry points read, in a fresh process.

A propagation code builds a few models at start-up; each must read its own records and
data files, not the catalogue (some 8,300 files, seconds on a cold disk or under an
on-access scanner). The bounds are about twice what each reads today, so a model that
grows a little passes and one that loads the catalogue fails.
"""

import subprocess
import sys

import pytest

COUNT = """
import sys
opened = 0
def hook(event, args):
    global opened
    path = args[0] if event == "open" and isinstance(args[0], str) else ""
    if "/optmat/data/" in path.replace("\\\\", "/"):
        opened += 1
sys.addaudithook(hook)
{code}
print(opened)
"""


@pytest.mark.parametrize(
    ("code", "bound"),
    [
        (
            "from optmat.gas.molecular_ionisation_rates import molecular_rate_model\n"
            "molecular_rate_model('N2', 800e-9)",
            80,
        ),
        (
            "from optmat.gas.ionisation_rates import ionisation_rate_model\n"
            "ionisation_rate_model('He', 800e-9)",
            80,
        ),
        (
            "from optmat.gas import ionisation_rates, molecular_ionisation_rates\n"
            "ionisation_rates.ionisation_rate_coverage()\n"
            "molecular_ionisation_rates.molecular_rate_coverage()",
            5,
        ),
        (
            "from optmat.gas.models.registry import default_model\ndefault_model('He')",
            60,
        ),
        (
            "from optmat.gas.raman_response import raman_response\n"
            "raman_response('N2', 293.15)",
            60,
        ),
        ("from optmat.gas.air import air_model\nair_model()", 50),
        (
            "from optmat.solid.models.registry import default_model\n"
            "default_model('SiO2'); default_model('N-BK7')",
            20,
        ),
        ("from optmat.solid.bundle import silica_bundle\nsilica_bundle()", 20),
        ("from optmat.datasets import sources\nsources()", 800),
    ],
    ids=[
        "molecular_rate",
        "atomic_rate",
        "rate_coverage",
        "gas_default",
        "raman",
        "air",
        "solid_defaults",
        "silica_bundle",
        "every_record",
    ],
)
def test_an_entry_point_reads_only_its_own_data_files(code, bound):
    result = subprocess.run(
        [sys.executable, "-c", COUNT.format(code=code)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert int(result.stdout.split()[-1]) <= bound

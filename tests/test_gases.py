# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

import pytest

from optmat.exceptions import OptmatError
from optmat.gas.gases import Gas, as_gas, info


def test_every_gas_has_info():
    for gas in Gas:
        assert info(gas).gas is gas


@pytest.mark.parametrize(
    ("gas", "electrons"),
    [
        ("He", 2),
        ("Ne", 10),
        ("Ar", 18),
        ("Kr", 36),
        ("Xe", 54),
        ("N2", 14),
        ("SF6", 70),
    ],
)
def test_electron_counts(gas, electrons):
    assert info(gas).electrons == electrons


def test_only_noble_gases_are_atomic():
    atomic = {gas for gas in Gas if info(gas).atomic}
    assert atomic == {Gas.He, Gas.Ne, Gas.Ar, Gas.Kr, Gas.Xe}


def test_as_gas_accepts_names_and_members():
    assert as_gas("Kr") is Gas.Kr
    assert as_gas(Gas.Kr) is Gas.Kr


def test_as_gas_rejects_unknown_species():
    with pytest.raises(OptmatError, match="unknown gas 'Og'"):
        as_gas("Og")

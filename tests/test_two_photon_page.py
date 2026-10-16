# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The manual page *Two-photon ionisation and absorption*: its tables and claims.

``tools/two_photon_tables.py`` writes the page's tables from the curated records (the
reconciliation of the two-photon curation); these tests check that the page holds its
output (regenerate with ``tools/two_photon_tables.py --update-docs``) and the numbers
the page's text states.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]


def _tool():
    path = ROOT / "tools" / "two_photon_tables.py"
    spec = importlib.util.spec_from_file_location("two_photon_tables", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("two_photon_tables", module)
    spec.loader.exec_module(module)
    return module


T = _tool()


def _block(name: str) -> str:
    text = "\n".join(path.read_text(encoding="utf-8") for path in T.PAGES)
    return text.split(f"<!-- {name}:start -->\n")[1].split(f"\n<!-- {name}:end -->")[0]


@pytest.mark.parametrize("name", sorted(T.tables()))
def test_the_page_holds_the_tool_output(name):
    assert _block(name) == T.tables()[name]


def test_the_thresholds_are_nists_halved():
    # Argon: 15.7596119 eV and the Ar+ 2P1/2 level 1431.5831 cm-1 (nist_asd_2026),
    # halved; 1e-9 relative is the arithmetic, not a physical uncertainty.
    t = T.thresholds("Ar")
    assert t["lower"] == pytest.approx(15.7596119 / 2, rel=1e-12, abs=0.0)
    assert t["upper"] == pytest.approx(
        (15.7596119 + 1431.5831 * T.HC_EV_CM) / 2, rel=1e-9, abs=0.0
    )
    assert t["upper"] == pytest.approx(7.9686, abs=5e-5)
    assert T.thresholds("He")["upper"] == T.thresholds("He")["lower"]


def test_mccowns_coefficient_converts_by_the_photon_energy():
    # McCown et al.'s Eq. 9: dn/dt = alpha [Xe] I^2 / hbar omega, so sigma = alpha
    # hbar omega: 4e-32 cm4 W-1 at 193 nm is 4.1e-50 cm4 s (the value Lagutin et al.
    # draw for it; mccown_1982). 1e-3: the printed 193 nm against the line centre.
    (value,) = T.cm4_s(T.dataset("mccown_1982/alpha_Xe"))
    assert value == pytest.approx(
        4e-32 * 1239.841984 / 193 * 1.602176634e-19, rel=1e-3, abs=0.0
    )


def test_argon_calculations_agree_within_the_stated_margin():
    # The page: Petrov et al. with core polarisation within 10 % of McKenna and van der
    # Hart from 8 to 11 eV (both gauges).
    _, curves = T.COMPARISONS["Ar"]
    for energy in (8.0, 9.0, 10.0, 11.0):
        base = T.curve_at(curves[0], energy)
        for curve in curves[1:3]:
            assert T.curve_at(curve, energy) / base == pytest.approx(1.0, abs=0.10)


def test_values_are_not_extrapolated():
    assert T.at("lagutin_2017/fig8_cipfcp_linear", 6.42) is None
    assert T.at("mccown_1982/alpha_Xe", 7.0) is None

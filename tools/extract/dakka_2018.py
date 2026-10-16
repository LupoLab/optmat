# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Dakka et al. (2018): Kr 2p6 two-photon and photoionisation cross sections.

The library holds the submitted version (arXiv:1805.05669v1, 15 May 2018), behind a
one-page repository cover (PDF p. 1). Table I (``tbl_599c680c31e3``, PDF p. 5; image
checked) prints, for 'This work', the branching ratio of the 2p6 level to the metastable
1s5 (75 +- 4 %, not a quantity of the vocabulary), the two-photon cross section
sigma_0^(2) of the ground-2p6 transition (5.2 +- 2.1, 10^-35 cm^4) and the
photoionisation cross section of the 2p6 level sigma_ei (3.8 +- 1.7, 10^-18 cm^2), both
fitted with the authors' density-matrix model; and a 'Literature' row (not extracted:
the branching ratio's average of refs. 26, 27, 31-39; Kroll and Bischel's xenon value;
a value 'calculated from' Hyman 1977).

:func:`checks` compares the 'Literature' sigma_0^(2), 2.4 +- 0.8, with kroll_1990's
2p6 value: it is the as-measured product sigma_tot G2(0), not the G2-corrected 1.7.
"""

from _library import Point, number, table, write

from optmat.datasets import DATA_DIRECTORY, load_source

DIRECTORY = "Dakka_2018_LaserBasedMetastable"
TABLE_I = "tbl_599c680c31e3"
SIGMA_0, SIGMA_EI = 2, 3
"""Table I's columns of sigma_0^(2) and sigma_ei."""
IONISING_NM = "215"
"""The laser that both excites and ionises the 2p6 level, 'a 215 nm laser' (p. 2),
'near 215 nm' (abstract): the wavelength as printed."""


def _split(cell: str) -> tuple[str, str]:
    """Return a ``value +- sigma`` cell's two numbers, as printed."""
    value, _, sigma = cell.partition("±")
    return number(value), number(sigma)


def checks() -> dict[str, object]:
    """Return the comparison of Table I's 'Literature' sigma_0^(2) with kroll_1990.

    Raises if the 'Literature' value is not Kroll and Bischel's as-measured 2p6 value.
    """
    rows = table(DIRECTORY, TABLE_I)
    literature = rows[2]
    printed = number(literature[SIGMA_0].split("±")[0])
    record = load_source(DATA_DIRECTORY / "sources" / "kroll_1990.toml")
    kroll = {d.id: d for d in record.datasets}
    measured = kroll["sigma_tot_g2"].values[kroll["sigma_tot_g2"].rows.index("2p6")]
    corrected = kroll["sigma_tot"].values[kroll["sigma_tot"].rows.index("2p6")]
    out: dict[str, object] = {
        "Literature sigma_0 (Xe 2p6, ref. 41)": printed,
        "kroll_1990 2p6: sigma_tot G2(0), sigma_tot": (measured, corrected),
    }
    if printed != measured:
        msg = f"Table I's literature value {printed} is not kroll_1990's {measured}"
        raise ValueError(msg)
    value, sigma = _split(rows[1][SIGMA_0])
    out["sigma_0 relative uncertainty"] = float(sigma) / float(value)
    return out


def extract() -> dict[str, list[Point]]:
    """Return Table I's two cross sections of this work."""
    rows = table(DIRECTORY, TABLE_I)
    this_work = rows[1]
    if this_work[0] != "This work":
        msg = "Table I's first row is not 'This work'"
        raise ValueError(msg)
    printed = " | ".join(this_work)
    sigma_0 = _split(this_work[SIGMA_0])
    sigma_ei = _split(this_work[SIGMA_EI])
    return {
        "Kr/dakka_2018__sigma0_2p6.csv": [Point("2p6", "2p6", *sigma_0, printed)],
        "Kr/dakka_2018__photoionization_2p6.csv": [
            Point(IONISING_NM, IONISING_NM, *sigma_ei, printed)
        ],
    }


def main() -> None:
    """Write the data files and print the checks."""
    write(extract())
    for label, result in checks().items():
        print(f"check {label}: {result}")


if __name__ == "__main__":
    main()

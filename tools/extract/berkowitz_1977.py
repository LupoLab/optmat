# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Berkowitz & Eland (1977): fragment-ion branching ratios of N2O.

The ionisation energies and appearance energies of this source are text values, inline
in the record; this script writes only the two branching-ratio tables (TC-965).

* **Table IV** (``tbl_130104e0aa54``, p. 2745, image checked): 'Branching ratios at
  selected wavelengths in the 740-840 A region', the ion intensities N2O+ : O+ :
  (15NO+ + 14NO+) with N2O+ = 100, at three wavelengths each printed with a "~" and
  "(peak)" or "(valley)". The wavelength is the number after the "~"; the cell stays
  whole in ``printed``.
* **Table VII** (``tbl_be3196b2f333``, p. 2748, image checked): 'Branching ratios in
  the 600-750 A region', N2O+ : NO+ : N2+ : N+ : O+, and 15NO+/14NO+. Footnote a:
  'Figures marked with an asterisk (*) refer to amplitudes of features, i.e., peak
  heights or step heights. A dash (-) indicates that a feature was too weak or to [sic]
  badly overlapped to measure.' So the table holds two kinds of ratio, written to
  separate files: the ion intensities at a wavelength (the first three rows, no
  asterisk: ``t7_yield_*``) and the amplitudes of the autoionising resonances over that
  of N2O+ (the next five rows, every cell with an asterisk: ``t7_feature_*``). The
  asterisk is removed from the value and kept in ``printed``; a dash is no value and
  writes no row.

  The last three rows (the B 2Pi regions 705-685 and 685-640 A and C 2Sigma+ (0,0,0) at
  618 A) have a dash for N2O+, so their numbers are on no stated reference; they are not
  written (the record's notes give them). The 15NO+/14NO+ column is a ratio between the
  two NO+ of the 14N15NO isotopologue, not a channel of N2O; it is not written either.

  The description of the 666.7 A row is extracted as "nso, n^{*}=3"; the image prints
  "ns sigma, n* = 3" (:data:`ROWS`).
"""

from _library import Point, number, table, write

DIRECTORY = "Berkowitz_1977_PhotoionizationN2OMechanisms"
TABLE_4 = "tbl_130104e0aa54"
TABLE_7 = "tbl_be3196b2f333"
SOURCE = "berkowitz_1977"

TABLE_4_ROWS = {
    "~798.3(peak)": ("798.3_peak", "798.3"),
    "~792(valley)": ("792_valley", "792"),
    "~782(peak)": ("782_peak", "782"),
}
"""Table IV's first cell: the row id and the wavelength, angstrom, as printed."""
TABLE_4_COLUMNS = {"O^{+}": "t4_Oplus", "^{15}NO^{+} + ^{14}NO^{+}": "t4_NOplus"}

ROWS = {
    ("600", "flat"): ("600_flat", "yield"),
    ("660", "valley"): ("660_valley", "yield"),
    ("668.3", "nd\\pi peak"): ("668.3_ndpi_peak", "yield"),
    ("668.3", "nd\\pi, n^{*}=2.95"): ("668.3_ndpi_n2.95", "feature"),
    ("644.3", "nd\\pi, n^{*}=3.97"): ("644.3_ndpi_n3.97", "feature"),
    ("633.9", "nd\\pi, n^{*}=4.95"): ("633.9_ndpi_n4.95", "feature"),
    ("648.9", "nd\\sigma, n^{*}=3.7"): ("648.9_ndsigma_n3.7", "feature"),
    # The image prints "ns sigma, n* = 3"; the extraction reads the sigma as "o".
    ("666.7", "nso, n^{*}=3"): ("666.7_nssigma_n3", "feature"),
}
"""Table VII's wavelength and description cells: the row id and the kind of ratio."""
NOT_WRITTEN = {"705-685", "685-640", "618"}
"""Table VII's rows with no N2O+ entry (the B and C states): no stated reference."""
TABLE_7_COLUMNS = {
    "NO^{+}": "NOplus",
    "N_{2}^{+}": "N2plus",
    "N^{+}": "Nplus",
    "O^{+}": "Oplus",
}
DASHES = {"-", "\N{MINUS SIGN}", "\N{EN DASH}"}


def _table_4() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_4)
    header = [cell.strip() for cell in grid[0]]
    if header[1] != "N_{2}O^{+}":
        msg = f"Table IV header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        f"N2O/{SOURCE}__{dataset}.csv": [] for dataset in TABLE_4_COLUMNS.values()
    }
    for cells in grid[1:]:
        row, x = TABLE_4_ROWS[cells[0].strip()]
        if number(cells[1]) != "100":
            msg = f"Table IV row {cells!r} is not on N2O+ = 100"
            raise ValueError(msg)
        printed = " | ".join(cells)
        for column, dataset in TABLE_4_COLUMNS.items():
            value = number(cells[header.index(column)])
            files[f"N2O/{SOURCE}__{dataset}.csv"].append(
                Point(row, x, value, "", printed)
            )
    return files


def _table_7() -> dict[str, list[Point]]:
    grid = table(DIRECTORY, TABLE_7)
    header = [cell.strip() for cell in grid[0]]
    if header[2] != "N_{2}O^{+}":
        msg = f"Table VII header {header!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        wavelength = cells[0].strip()
        if wavelength in NOT_WRITTEN:
            if cells[2].strip() not in DASHES:
                msg = f"Table VII row {cells!r} has an N2O+ entry"
                raise ValueError(msg)
            continue
        row, kind = ROWS[wavelength, cells[1].strip()]
        marked = kind == "feature"
        parent = cells[2].strip()
        if parent.endswith("*") != marked or number(parent.removesuffix("*")) != "100":
            msg = f"Table VII row {cells!r}: N2O+ is not {'100*' if marked else '100'}"
            raise ValueError(msg)
        printed = " | ".join(cells)
        for column, ion in TABLE_7_COLUMNS.items():
            cell = cells[header.index(column)].strip()
            if cell in DASHES:
                continue
            if cell.endswith("*") != marked:
                msg = f"Table VII row {row}, {column}: {cell!r} breaks its asterisks"
                raise ValueError(msg)
            value = number(cell.removesuffix("*"))
            relative = f"N2O/{SOURCE}__t7_{kind}_{ion}.csv"
            files.setdefault(relative, []).append(
                Point(row, wavelength, value, "", printed)
            )
    return files


def extract() -> dict[str, list[Point]]:
    """Return the branching ratios of Tables IV and VII, one file per ion and kind."""
    return _table_4() | _table_7()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

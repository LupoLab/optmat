# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Hagemann, Gudat & Kunz (1974), DESY SR-74/7: Tables 2-6 and 9.

The report's Tables 1-9 are line-printer output of the Kramers-Kronig analysis
("Numerical values of the optical functions from the Kramers-Kronig-analysis of
absorption data, N-EFF according to Eq. 6", PDF page 36), one table per substance, each
over two or three pages, every page headed "OPTICAL CONSTANTS OF <substance>" and
labelled in the margin "Table N, page M". The pages and their order, from the images
(the library's cells carry no captions; its titles are misread, "CL" for CU, "AL" for
AU, "AI" for BI):

- Table 2, Al: PDF pp. 40-42, 1.00E-03 to 1.20E+05 eV;
- Table 3, Cu: pp. 43-45, 2.00E-03 to 5.00E+04 eV;
- Table 4, Ag: pp. 46-48, 5.00E-03 to 5.00E+05 eV;
- Table 5, Au: pp. 49-51, 1.00E-03 to 1.50E+05 eV;
- Table 6, Au "FROM REFLECTIVITY-DATA FITTED TO CANFIELD ET AL (10)": pp. 52-54, 1.50 to
  350 eV;
- Table 9, Al2O3 (evaporated, amorphous): pp. 60-61, 6.00 to 1.60E+03 eV.

(Table 1, Mg: pp. 37-39; Table 7, Bi: pp. 55-57; Table 8, C: pp. 58-59; not extracted.)

Every page prints nine columns: ENERGY (eV), EPSILON1, EPSILON2, N, K, ABSCOEFF (mu =
2 omega k / c in cm-1, checked below), REFLECT% (normal-incidence reflectance in
percent), ENERGYLOSS (-Im 1/epsilon) and N-EFF (Eq. 6). The last page of each metal's
table (from 65 eV for Cu to 122 eV for Al) prints "1.-EPSL1" and "1.-N" (1 - epsilon_1
and 1 - n, three decimals) instead of EPSILON1 and N: the n and epsilon_1 datasets end
there, "1.-N" becomes the ``delta`` dataset (refractive_index_decrement, delta = 1 - n,
as printed) and "1.-EPSL1", which has no quantity, stays in ``printed``. Both are
checked with the rest of the row (n = 1 - delta, epsilon_1 = 1 - printed value).
Numbers are printed "5.19E-03" and "1.17E 05": ``printed`` keeps them so, and the
stored number writes the positive exponent "E+05". Energies are the row ids.

The library's cells (JSON; four pages also have a CSV) are an optical reading of a
faint printout, and are corrected here in two structural steps, each against the page
images (the 300-dpi render of the PDF, every cell of these 17 pages compared):

1. The printout's zeros are often read as the letters C, O, Q, D, U (and a comma for
   the decimal point, F for E): ``cell_text`` reads them as printed.
2. ``CORRECTIONS``: 396 cells the library misread otherwise (9 read as 5, 0 as 9, 8 as
   6, 1 as 0 ...), each as (library cell, cell on the page). Most were found by the
   redundancy checks below and then read on the image; the N-EFF, REFLECT% and
   ENERGYLOSS cells, which nothing else checks tightly, were read on the image row by
   row. On the faint last page of Table 4 (PDF p. 48) a few digits cannot be read; they
   are written "?" in ``printed`` (N-EFF, one reflectance and one loss cell only; no
   dataset value is affected). On PDF p. 45 (Table 3, 78 eV) the "1.-EPSL1" cell
   looks like "0.043" on the 300-dpi render; at 600 dpi the left halves of that row's
   glyphs are faded (its leading 0 reads as a 3) and the last digit is an 8 with its
   left half missing: 0.048, which 1 - n = 0.020 and k give (1 - epsilon_1 = 0.048).

The redundant columns are the check: for every row, epsilon_1 = n^2 - k^2, epsilon_2 =
2nk, R = |(N - 1)/(N + 1)|^2, -Im 1/epsilon = epsilon_2/|epsilon|^2 and mu = 2 E k /
(hbar c) must hold to the printed rounding (interval arithmetic over the last printed
digits) within ``TOLERANCE`` (3 %) beyond it; energies must increase and N-EFF must not
decrease. The rows that fail on the page itself are listed in ``INCONSISTENT`` with the
failing checks, and the script insists the list is exact. They are printed so (read
again on the image), not misread; the report computed each function on the analysis
grid and printed them at the table energies, and the largest departures are at the
first rows (the Drude extrapolation) and at steep structure:

- Table 2, Al: at 72.9 eV (p. 40) and 73.1-80 eV (p. 41) the ABSCOEFF and ENERGYLOSS
  cells are typed in another typeface (pasted in); mu there is 4-13 % below
  2 E k / (hbar c) from the printed k, which is consistent with epsilon_2. 3.00 eV:
  epsilon_2 printed 5.68 against 2nk = 5.48 (3 %). 14-16 eV (the plasmon): loss 4-10 %
  off.
- Table 3, Cu: 2 meV, mu 13 % below 2 E k / (hbar c) (epsilon_1 and n are blank in that
  row on the page).
- Table 4, Ag: 25 meV, epsilon_2 4 % above 2nk; 3.80 eV, R 7 % above.
- Table 5, Au: 5 meV and 0.1 eV, epsilon, loss and mu 4-20 % off n and k; 2.2 eV, loss
  5 %.
- Table 9, Al2O3: 6 eV, epsilon_2, loss and mu 9-28 % off n and k (k = 2.20E-03 at the
  absorption edge); 77.5 and 80 eV (the Al L edge), loss; 250-500 eV, mu 3-8 % below.

Table 6 has no row beyond 3 %. Below 3 %, the worst rows per table are in the record's
notes.
"""

import itertools
import re
from collections.abc import Callable, Sequence
from typing import Final

from _json_table import json_table
from _library import Point, write

DIRECTORY: Final = "Hagemann_1974_OpticalConstantsFar"

# Table -> (material, dataset suffix, pages in order: library table ids).
TABLES: Final = {
    "2": ("Al", "Al", ("tbl_5eadb4abaadf", "tbl_e9bd35b2f535", "tbl_7ad069048df5")),
    "3": ("Cu", "Cu", ("tbl_52ae5a332285", "tbl_163d5f56d649", "tbl_a4df553b9b72")),
    "4": ("Ag", "Ag", ("tbl_5954ad9a5d51", "tbl_a9377026a90a", "tbl_44c1d033263c")),
    "5": ("Au", "Au", ("tbl_e95227627cb4", "tbl_62a5fb0c7bc4", "tbl_a46c4afa9069")),
    "6": (
        "Au",
        "Au_table6",
        ("tbl_376d0c04b097", "tbl_28efb140cb95", "tbl_5ab449433de5"),
    ),
    "9": ("Al2O3_amorphous", "Al2O3", ("tbl_7cb0a953178f", "tbl_73348bbbe4f6")),
}

# The printout's columns.
ENERGY, EPS1, EPS2, N, K, MU, REFLECT, LOSS, NEFF = range(9)
# Dataset kind -> column: n and epsilon_1 where the page prints them, delta (the
# "1.-N" column) where it prints 1 - n instead.
DATASETS: Final = {"n": N, "k": K, "eps1": EPS1, "eps2": EPS2, "mu": MU, "delta": N}
PLAIN_ONLY: Final = frozenset({"n", "eps1"})
COMPLEMENT_ONLY: Final = frozenset({"delta"})

TOLERANCE: Final = 0.03
"""Relative departure allowed beyond the printed rounding in the redundancy checks."""

HBARC_EV_CM: Final = 1.973269804e-5
"""hbar c in eV cm (CODATA 2018), for mu = 2 E k / (hbar c) in cm-1."""

# (library table id, data row, column) -> (library cell, cell on the page).
CORRECTIONS: Final[dict[tuple[str, int, int], tuple[str, str]]] = {
    # PDF page 40, Table 2 page 1 (tbl_5eadb4abaadf)
    ("tbl_5eadb4abaadf", 2, 4): ("3.82E 02", "3.83E 02"),
    ("tbl_5eadb4abaadf", 5, 2): ("1.05E 03", "1.09E 03"),
    ("tbl_5eadb4abaadf", 8, 3): ("1.546", "1.946"),
    ("tbl_5eadb4abaadf", 12, 6): ("96.930", "96.980"),
    ("tbl_5eadb4abaadf", 15, 2): ("2.57E 01", "2.97E 01"),
    ("tbl_5eadb4abaadf", 19, 5): ("1.43E 06", "1.45E 06"),
    ("tbl_5eadb4abaadf", 19, 8): ("2.178", "2.108"),
    ("tbl_5eadb4abaadf", 23, 8): ("2.393", "2.383"),
    ("tbl_5eadb4abaadf", 25, 2): ("3.75E-01", "3.79E-01"),
    ("tbl_5eadb4abaadf", 26, 2): ("1.58E-01", "1.98E-01"),
    ("tbl_5eadb4abaadf", 27, 0): ("1.19E 01", "1.10E 01"),
    ("tbl_5eadb4abaadf", 34, 0): ("1.65E 01", "1.60E 01"),
    ("tbl_5eadb4abaadf", 37, 8): ("2.639", "2.638"),
    ("tbl_5eadb4abaadf", 38, 6): ("0.006", "0.086"),
    ("tbl_5eadb4abaadf", 42, 0): ("6.73E 01", "6.70E 01"),
    ("tbl_5eadb4abaadf", 43, 0): ("6.93E 01", "6.90E 01"),
    ("tbl_5eadb4abaadf", 45, 0): ("7.13E 01", "7.10E 01"),
    ("tbl_5eadb4abaadf", 47, 8): ("2.934", "2.834"),
    ("tbl_5eadb4abaadf", 48, 8): ("2.935", "2.835"),
    ("tbl_5eadb4abaadf", 49, 8): ("2.938", "2.838"),
    # PDF page 41, Table 2 page 2 (tbl_e9bd35b2f535)
    ("tbl_e9bd35b2f535", 1, 8): ("2.948", "2.848"),
    ("tbl_e9bd35b2f535", 2, 2): ("4.52E-02", "4.92E-02"),
    ("tbl_e9bd35b2f535", 17, 2): ("6.75E-02", "6.79E-02"),
    ("tbl_e9bd35b2f535", 26, 8): ("3.876", "3.870"),
    ("tbl_e9bd35b2f535", 33, 6): ("0.031", "0.021"),
    ("tbl_e9bd35b2f535", 36, 8): ("4.391", "4.381"),
    ("tbl_e9bd35b2f535", 37, 2): ("5.57E-02", "5.97E-02"),
    ("tbl_e9bd35b2f535", 47, 2): ("5.55E-02", "5.95E-02"),
    ("tbl_e9bd35b2f535", 47, 8): ("4.967", "4.957"),
    # PDF page 42, Table 2 page 3 (tbl_7ad069048df5)
    ("tbl_7ad069048df5", 0, 2): ("5.75E-02", "5.79E-02"),
    ("tbl_7ad069048df5", 8, 2): ("3.99E-02", "3.89E-02"),
    ("tbl_7ad069048df5", 12, 8): ("6.906", "6.806"),
    ("tbl_7ad069048df5", 14, 8): ("7.409", "7.499"),
    ("tbl_7ad069048df5", 16, 0): ("2.02E 02", "2.00E 02"),
    ("tbl_7ad069048df5", 20, 0): ("2.89E 02", "2.80E 02"),
    ("tbl_7ad069048df5", 21, 7): ("5.70E-03", "5.20E-03"),
    ("tbl_7ad069048df5", 28, 0): ("1.55E 03", "1.50E 03"),
    ("tbl_7ad069048df5", 28, 8): ("11.884", "11.886"),
    ("tbl_7ad069048df5", 32, 0): ("1.73E 03", "1.70E 03"),
    ("tbl_7ad069048df5", 33, 5): ("6.58E 03", "6.55E 03"),
    ("tbl_7ad069048df5", 33, 7): ("6.67E-05", "6.47E-05"),
    ("tbl_7ad069048df5", 38, 0): ("1.12E 04", "1.10E 04"),
    ("tbl_7ad069048df5", 40, 2): ("5.38E-09", "9.38E-09"),
    ("tbl_7ad069048df5", 42, 7): ("7.59E-10", "7.50E-10"),
    ("tbl_7ad069048df5", 47, 0): ("1.23E 05", "1.20E 05"),
    # PDF page 43, Table 3 page 1 (tbl_52ae5a332285)
    ("tbl_52ae5a332285", 1, 2): ("1.55E 05", "1.59E 05"),
    ("tbl_52ae5a332285", 2, 3): ("1.85E 02", "1.86E 02"),
    ("tbl_52ae5a332285", 4, 8): ("0.370", "0.310"),
    ("tbl_52ae5a332285", 8, 2): ("1.58E 00", "1.98E 00"),
    ("tbl_52ae5a332285", 13, 6): ("91.026", "91.006"),
    ("tbl_52ae5a332285", 13, 8): ("0.541", "0.561"),
    ("tbl_52ae5a332285", 14, 8): ("0.553", "0.563"),
    ("tbl_52ae5a332285", 18, 6): ("57.793", "57.703"),
    ("tbl_52ae5a332285", 24, 1): ("-1.475", "-1.476"),
    ("tbl_52ae5a332285", 26, 6): ("33.759", "33.569"),
    ("tbl_52ae5a332285", 34, 2): ("3.65E 00", "3.69E 00"),
    ("tbl_52ae5a332285", 37, 5): ("8.92E 05", "8.52E 05"),
    ("tbl_52ae5a332285", 39, 1): ("-0.067", "-0.007"),
    ("tbl_52ae5a332285", 43, 2): ("1.76E 00", "1.70E 00"),
    # PDF page 44, Table 3 page 2 (tbl_163d5f56d649)
    ("tbl_163d5f56d649", 3, 2): ("9.58E-01", "9.98E-01"),
    ("tbl_163d5f56d649", 4, 2): ("8.56E-01", "8.96E-01"),
    ("tbl_163d5f56d649", 7, 1): ("0.495", "0.695"),
    ("tbl_163d5f56d649", 8, 2): ("6.59E-01", "6.90E-01"),
    ("tbl_163d5f56d649", 15, 7): ("6.58E-01", "6.98E-01"),
    ("tbl_163d5f56d649", 16, 4): ("2.59E-01", "2.39E-01"),
    ("tbl_163d5f56d649", 17, 2): ("3.58E-01", "3.98E-01"),
    ("tbl_163d5f56d649", 27, 2): ("2.56E-01", "2.96E-01"),
    ("tbl_163d5f56d649", 28, 8): ("4.193", "4.183"),
    ("tbl_163d5f56d649", 40, 6): ("0.338", "0.398"),
    ("tbl_163d5f56d649", 46, 2): ("2.CBE-01", "2.08E-01"),
    ("tbl_163d5f56d649", 47, 8): ("5.490", "5.480"),
    ("tbl_163d5f56d649", 48, 2): ("2.CCE-01", "2.01E-01"),
    # PDF page 45, Table 3 page 3 (tbl_a4df553b9b72)
    ("tbl_a4df553b9b72", 0, 8): ("5.695", "5.685"),
    ("tbl_a4df553b9b72", 1, 0): ("6.62E 01", "6.60E 01"),
    ("tbl_a4df553b9b72", 7, 2): ("1.66E-01", "1.69E-01"),
    ("tbl_a4df553b9b72", 9, 8): ("6.293", "6.283"),
    ("tbl_a4df553b9b72", 11, 2): ("1.76E-01", "1.79E-01"),
    ("tbl_a4df553b9b72", 13, 1): ("0.043", "0.048"),
    ("tbl_a4df553b9b72", 15, 0): ("9.00E 01", "8.00E 01"),
    ("tbl_a4df553b9b72", 16, 0): ("9.10E 01", "8.10E 01"),
    ("tbl_a4df553b9b72", 17, 0): ("9.20E 01", "8.20E 01"),
    ("tbl_a4df553b9b72", 18, 0): ("9.30E 01", "8.30E 01"),
    ("tbl_a4df553b9b72", 19, 0): ("9.40E 01", "8.40E 01"),
    ("tbl_a4df553b9b72", 26, 5): ("6.55E 05", "6.65E 05"),
    ("tbl_a4df553b9b72", 28, 8): ("9.695", "9.605"),
    ("tbl_a4df553b9b72", 29, 2): ("6.CIE-02", "6.01E-02"),
    ("tbl_a4df553b9b72", 30, 8): ("13.375", "13.325"),
    ("tbl_a4df553b9b72", 31, 2): ("1.CIE-02", "1.01E-02"),
    ("tbl_a4df553b9b72", 32, 8): ("16.531", "16.931"),
    ("tbl_a4df553b9b72", 36, 2): ("5.24E-04", "5.26E-04"),
    ("tbl_a4df553b9b72", 37, 2): ("3.45E-04", "3.49E-04"),
    ("tbl_a4df553b9b72", 37, 6): ("0.00", "0.0"),
    ("tbl_a4df553b9b72", 41, 8): ("23.089", "23.080"),
    ("tbl_a4df553b9b72", 42, 2): ("1.5CE-04", "1.90E-04"),
    ("tbl_a4df553b9b72", 47, 3): ("-0.0C1", "-0.000"),
    ("tbl_a4df553b9b72", 48, 1): ("-0.0C1", "-0.000"),
    ("tbl_a4df553b9b72", 48, 2): ("6.15E-08", "6.19E-08"),
    ("tbl_a4df553b9b72", 48, 3): ("-0.0C1", "-0.000"),
    ("tbl_a4df553b9b72", 49, 1): ("0.0C1", "0.000"),
    ("tbl_a4df553b9b72", 49, 3): ("0.0C1", "0.000"),
    # PDF page 46, Table 4 page 1 (tbl_5954ad9a5d51)
    ("tbl_5954ad9a5d51", 2, 0): ("2.53E-02", "2.50E-02"),
    ("tbl_5954ad9a5d51", 2, 8): ("0.445", "0.465"),
    ("tbl_5954ad9a5d51", 4, 2): ("1.75E 03", "1.79E 03"),
    ("tbl_5954ad9a5d51", 8, 3): ("0.665", "0.669"),
    ("tbl_5954ad9a5d51", 8, 6): ("99.298", "99.208"),
    ("tbl_5954ad9a5d51", 15, 2): ("5.54E-01", "5.99E-01"),
    ("tbl_5954ad9a5d51", 16, 0): ("3.62E 00", "3.60E 00"),
    ("tbl_5954ad9a5d51", 16, 1): ("-1.219", "-1.218"),
    ("tbl_5954ad9a5d51", 17, 1): ("-6.497", "-0.497"),
    ("tbl_5954ad9a5d51", 20, 0): ("3.50E CC", "3.90E 00"),
    ("tbl_5954ad9a5d51", 21, 2): ("1.54E 00", "1.94E 00"),
    ("tbl_5954ad9a5d51", 23, 0): ("4.23E CC", "4.20E 00"),
    ("tbl_5954ad9a5d51", 24, 0): ("4.33E CC", "4.30E 00"),
    ("tbl_5954ad9a5d51", 28, 6): ("25.696", "25.695"),
    ("tbl_5954ad9a5d51", 31, 2): ("2.45E 00", "2.49E 00"),
    ("tbl_5954ad9a5d51", 32, 2): ("2.68E 00", "2.08E 00"),
    ("tbl_5954ad9a5d51", 34, 8): ("1.690", "1.600"),
    ("tbl_5954ad9a5d51", 38, 0): ("1.33E 01", "1.30E 01"),
    ("tbl_5954ad9a5d51", 39, 0): ("1.43E 01", "1.40E 01"),
    ("tbl_5954ad9a5d51", 39, 2): ("2.65E 00", "2.69E 00"),
    ("tbl_5954ad9a5d51", 42, 2): ("2.55E 00", "2.59E 00"),
    ("tbl_5954ad9a5d51", 43, 1): ("1.322", "1.022"),
    ("tbl_5954ad9a5d51", 43, 4): ("9.64E-01", "8.64E-01"),
    ("tbl_5954ad9a5d51", 45, 0): ("1.93E 01", "1.90E 01"),
    ("tbl_5954ad9a5d51", 45, 2): ("1.95E 00", "1.90E 00"),
    ("tbl_5954ad9a5d51", 47, 0): ("2.13E 01", "2.10E 01"),
    # PDF page 47, Table 4 page 2 (tbl_a9377026a90a)
    ("tbl_a9377026a90a", 0, 7): ("3.89E-C1", "3.88E-01"),
    ("tbl_a9377026a90a", 1, 4): ("9.42E-01", "9.43E-01"),
    ("tbl_a9377026a90a", 1, 5): ("2.22E 06", "2.20E 06"),
    ("tbl_a9377026a90a", 4, 8): ("5.865", "5.965"),
    ("tbl_a9377026a90a", 5, 2): ("1.52E CC", "1.57E 00"),
    ("tbl_a9377026a90a", 6, 6): ("14,357", "14.367"),
    ("tbl_a9377026a90a", 11, 2): ("1.05E CC", "1.06E 00"),
    ("tbl_a9377026a90a", 12, 4): ("5.74E-01", "5.72E-01"),
    ("tbl_a9377026a90a", 13, 3): ("0.915", "0.919"),
    ("tbl_a9377026a90a", 14, 8): ("7.359", "7.369"),
    ("tbl_a9377026a90a", 17, 7): ("8.COE-C1", "8.03E-01"),
    ("tbl_a9377026a90a", 18, 8): ("9.332", "8.332"),
    ("tbl_a9377026a90a", 19, 8): ("9.553", "8.553"),
    ("tbl_a9377026a90a", 21, 2): ("7.CBE-01", "7.08E-01"),
    ("tbl_a9377026a90a", 23, 2): ("5.77E-01", "6.77E-01"),
    ("tbl_a9377026a90a", 29, 2): ("5.45E-01", "5.49E-01"),
    ("tbl_a9377026a90a", 32, 2): ("4.95E-01", "4.99E-01"),
    ("tbl_a9377026a90a", 38, 7): ("5.85E-01", "5.86E-01"),
    ("tbl_a9377026a90a", 41, 3): ("0.889", "0.883"),
    ("tbl_a9377026a90a", 42, 2): ("3.62E-01", "3.63E-01"),
    ("tbl_a9377026a90a", 44, 2): ("3.CBE-01", "3.08E-01"),
    ("tbl_a9377026a90a", 45, 3): ("0.852", "0.853"),
    ("tbl_a9377026a90a", 45, 7): ("4.57E-01", "4.97E-01"),
    ("tbl_a9377026a90a", 46, 2): ("2.65E-01", "2.69E-01"),
    ("tbl_a9377026a90a", 47, 6): ("1,237", "1.287"),
    ("tbl_a9377026a90a", 49, 2): ("1.68E-01", "1.88E-01"),
    ("tbl_a9377026a90a", 49, 5): ("9,55E 05", "9.55E 05"),
    # PDF page 48, Table 4 page 3 (tbl_44c1d033263c)
    ("tbl_44c1d033263c", 0, 6): ("0.897", "0.867"),
    ("tbl_44c1d033263c", 1, 0): ("9.5JE C1", "9.50E 01"),
    ("tbl_44c1d033263c", 3, 0): ("1.COE 02", "1.05E 02"),
    ("tbl_44c1d033263c", 4, 1): ("0.166", "0.186"),
    ("tbl_44c1d033263c", 4, 3): ("0.058", "0.098"),
    ("tbl_44c1d033263c", 5, 6): ("0.139", "0.1?9"),
    ("tbl_44c1d033263c", 7, 7): ("6.15E-03", "8.15E-03"),
    ("tbl_44c1d033263c", 9, 2): ("8.45E-03", "8.49E-03"),
    ("tbl_44c1d033263c", 10, 0): ("1.7JE 02", "1.70E 02"),
    ("tbl_44c1d033263c", 10, 2): ("8.55E-03", "8.59E-03"),
    ("tbl_44c1d033263c", 12, 1): ("0.035", "0.036"),
    ("tbl_44c1d033263c", 12, 2): ("7.46E-03", "7.40E-03"),
    ("tbl_44c1d033263c", 12, 7): ("7.55E-C3", "7.95E-03"),
    ("tbl_44c1d033263c", 12, 8): ("20.357", "20.397"),
    ("tbl_44c1d033263c", 14, 6): ("0.013", "0.003"),
    ("tbl_44c1d033263c", 14, 8): ("20.379", "20.879"),
    ("tbl_44c1d033263c", 16, 0): ("2.EQE 02", "2.80E 02"),
    ("tbl_44c1d033263c", 16, 7): ("5.95E-C3", "5.50E-03"),
    ("tbl_44c1d033263c", 17, 1): ("0.C04", "0.008"),
    ("tbl_44c1d033263c", 17, 8): ("21.473", "21.573"),
    ("tbl_44c1d033263c", 18, 2): ("3.34E-03", "3.36E-03"),
    ("tbl_44c1d033263c", 18, 8): ("22.256", "22.056"),
    ("tbl_44c1d033263c", 19, 0): ("3.7JE 02", "3.70E 02"),
    ("tbl_44c1d033263c", 19, 4): ("3.23E-03", "8.23E-03"),
    ("tbl_44c1d033263c", 20, 8): ("24.255", "24.??5"),
    ("tbl_44c1d033263c", 21, 1): ("0.C03", "0.006"),
    ("tbl_44c1d033263c", 23, 2): ("2.62E-03", "2.83E-03"),
    ("tbl_44c1d033263c", 23, 5): ("1.05E 05", "1.06E 05"),
    ("tbl_44c1d033263c", 23, 7): ("2.65E-C3", "2.85E-03"),
    ("tbl_44c1d033263c", 24, 1): ("0.C03", "0.005"),
    ("tbl_44c1d033263c", 25, 2): ("1.64E-03", "1.04E-03"),
    ("tbl_44c1d033263c", 25, 7): ("1.C5F-C3", "1.05E-03"),
    ("tbl_44c1d033263c", 25, 8): ("34.256", "36.0?4"),
    ("tbl_44c1d033263c", 26, 8): ("36.239", "38.2??"),
    ("tbl_44c1d033263c", 28, 2): ("4.75E-05", "4.79E-05"),
    ("tbl_44c1d033263c", 28, 5): ("6.29E 03", "6.05E 03"),
    ("tbl_44c1d033263c", 28, 6): ("0.000", "0.0"),
    ("tbl_44c1d033263c", 29, 2): ("2.85E-05", "2.86E-05"),
    ("tbl_44c1d033263c", 29, 5): ("4.31E 03", "4.34E 03"),
    ("tbl_44c1d033263c", 31, 2): ("6.55E-05", "6.50E-05"),
    ("tbl_44c1d033263c", 32, 8): ("41.613", "41.513"),
    ("tbl_44c1d033263c", 33, 7): ("4.42E-05", "4.40E-05"),
    ("tbl_44c1d033263c", 34, 3): ("-0.C001", "-0.001"),
    ("tbl_44c1d033263c", 34, 7): ("3.63E-05", "3.03E-05"),
    ("tbl_44c1d033263c", 35, 7): ("1.65E-05", "1.66E-05"),
    ("tbl_44c1d033263c", 35, 8): ("43.463", "43.863"),
    ("tbl_44c1d033263c", 36, 2): ("5.65E-06", "5.85E-06"),
    ("tbl_44c1d033263c", 36, 6): ("0.000", "0.0"),
    ("tbl_44c1d033263c", 36, 7): ("5.89E-C4", "5.85E-06"),
    ("tbl_44c1d033263c", 36, 8): ("44.467", "44.342"),
    ("tbl_44c1d033263c", 37, 3): ("0.C001", "0.001"),
    ("tbl_44c1d033263c", 37, 7): ("2.59E-C6", "2.58E-06"),
    ("tbl_44c1d033263c", 38, 6): ("0.000", "0.0"),
    ("tbl_44c1d033263c", 38, 8): ("45.441", "46.0?1"),
    ("tbl_44c1d033263c", 39, 0): ("2.59E 04", "2.50E 04"),
    ("tbl_44c1d033263c", 39, 8): ("48.181", "46.1?1"),
    ("tbl_44c1d033263c", 40, 0): ("2.COE 04", "2.60E 04"),
    ("tbl_44c1d033263c", 40, 8): ("49.261", "46.261"),
    ("tbl_44c1d033263c", 41, 5): ("3.59E 02", "3.99E 02"),
    ("tbl_44c1d033263c", 41, 6): ("0.000", "0.0"),
    ("tbl_44c1d033263c", 41, 8): ("49.540", "46.530"),
    ("tbl_44c1d033263c", 42, 1): ("-0.CC3", "-0.005"),
    ("tbl_44c1d033263c", 42, 8): ("49.967", "46.907"),
    ("tbl_44c1d033263c", 43, 1): ("-0.CC3", "-0.005"),
    ("tbl_44c1d033263c", 43, 2): ("3.66E-C8", "3.86E-08"),
    ("tbl_44c1d033263c", 43, 7): ("3.82E-C6", "3.82E-08"),
    ("tbl_44c1d033263c", 44, 0): ("7.5JE 04", "7.50E 04"),
    ("tbl_44c1d033263c", 44, 1): ("-0.CC3", "-0.000"),
    ("tbl_44c1d033263c", 44, 6): ("0.000", "0.0"),
    ("tbl_44c1d033263c", 44, 7): ("6.67E-09", "?.67E-09"),
    ("tbl_44c1d033263c", 45, 2): ("2.57E-C9", "2.97E-09"),
    ("tbl_44c1d033263c", 45, 7): ("2.59E-C9", "2.98E-09"),
    ("tbl_44c1d033263c", 46, 2): ("2.26E-1C", "2.29E-10"),
    ("tbl_44c1d033263c", 47, 1): ("-0.CC3", "-0.006"),
    ("tbl_44c1d033263c", 47, 2): ("7.54E-12", "7.94E-12"),
    # PDF page 49, Table 5 page 1 (tbl_e95227627cb4)
    ("tbl_e95227627cb4", 8, 8): ("0.925", "0.825"),
    ("tbl_e95227627cb4", 10, 3): ("0.876", "0.878"),
    ("tbl_e95227627cb4", 10, 8): ("0.948", "0.848"),
    ("tbl_e95227627cb4", 11, 8): ("0.366", "0.866"),
    ("tbl_e95227627cb4", 17, 1): ("0.823", "0.826"),
    ("tbl_e95227627cb4", 18, 3): ("2.C35", "2.039"),
    ("tbl_e95227627cb4", 18, 4): ("1.96E 00", "1.86E 00"),
    ("tbl_e95227627cb4", 19, 3): ("1.955", "1.995"),
    ("tbl_e95227627cb4", 22, 7): ("1.65E-01", "1.66E-01"),
    ("tbl_e95227627cb4", 26, 2): ("3.57E CC", "3.97E 00"),
    ("tbl_e95227627cb4", 28, 8): ("1.981", "1.881"),
    ("tbl_e95227627cb4", 32, 2): ("3.35E CC", "3.39E 00"),
    ("tbl_e95227627cb4", 34, 0): ("1.CO E C1", "1.05E 01"),
    ("tbl_e95227627cb4", 35, 0): ("1.19E 01", "1.10E 01"),
    ("tbl_e95227627cb4", 35, 2): ("2.54E 00", "2.94E 00"),
    ("tbl_e95227627cb4", 36, 3): ("1.565", "1.569"),
    ("tbl_e95227627cb4", 41, 6): ("14.938", "14.988"),
    ("tbl_e95227627cb4", 41, 8): ("4.362", "4.302"),
    ("tbl_e95227627cb4", 45, 2): ("3.36E CC", "3.38E 00"),
    ("tbl_e95227627cb4", 46, 0): ("2.CO E 01", "2.10E 01"),
    ("tbl_e95227627cb4", 46, 3): ("1.495", "1.499"),
    ("tbl_e95227627cb4", 47, 3): ("1.426", "1.406"),
    ("tbl_e95227627cb4", 48, 0): ("2.CO E 01", "2.20E 01"),
    # PDF page 50, Table 5 page 2 (tbl_62a5fb0c7bc4)
    ("tbl_62a5fb0c7bc4", 0, 2): ("2.65E 00", "2.69E 00"),
    ("tbl_62a5fb0c7bc4", 8, 3): ("0.867", "0.887"),
    ("tbl_62a5fb0c7bc4", 11, 2): ("1.25E 00", "1.29E 00"),
    ("tbl_62a5fb0c7bc4", 14, 2): ("9.61E-01", "9.01E-01"),
    ("tbl_62a5fb0c7bc4", 22, 2): ("6.C5E-01", "6.09E-01"),
    ("tbl_62a5fb0c7bc4", 28, 7): ("5.63E-01", "9.63E-01"),
    ("tbl_62a5fb0c7bc4", 29, 6): ("3.651", "3.851"),
    ("tbl_62a5fb0c7bc4", 29, 7): ("6.13E-01", "9.13E-01"),
    ("tbl_62a5fb0c7bc4", 32, 2): ("3.55E-01", "3.95E-01"),
    ("tbl_62a5fb0c7bc4", 37, 1): ("0.643", "0.648"),
    ("tbl_62a5fb0c7bc4", 38, 2): ("3.CEE-01", "3.08E-01"),
    ("tbl_62a5fb0c7bc4", 39, 2): ("2.62E-01", "2.82E-01"),
    ("tbl_62a5fb0c7bc4", 44, 2): ("1.95E-01", "1.99E-01"),
    ("tbl_62a5fb0c7bc4", 44, 6): ("1.091", "1.051"),
    ("tbl_62a5fb0c7bc4", 45, 2): ("1.66E-01", "1.86E-01"),
    ("tbl_62a5fb0c7bc4", 48, 5): ("7.58E 05", "7.98E 05"),
    # PDF page 51, Table 5 page 3 (tbl_a46c4afa9069)
    ("tbl_a46c4afa9069", 3, 2): ("1.55E-01", "1.50E-01"),
    ("tbl_a46c4afa9069", 9, 2): ("9.55E-02", "9.99E-02"),
    ("tbl_a46c4afa9069", 9, 5): ("5.45E 05", "5.49E 05"),
    ("tbl_a46c4afa9069", 15, 2): ("2.66E-02", "2.86E-02"),
    ("tbl_a46c4afa9069", 17, 8): ("21.962", "21.952"),
    ("tbl_a46c4afa9069", 18, 1): ("0.083", "0.088"),
    ("tbl_a46c4afa9069", 23, 2): ("2.1CE-02", "2.18E-02"),
    ("tbl_a46c4afa9069", 37, 2): ("2.5CE-04", "2.90E-04"),
    ("tbl_a46c4afa9069", 37, 7): ("2.50E-04", "2.90E-04"),
    ("tbl_a46c4afa9069", 38, 2): ("3.26E-04", "3.28E-04"),
    ("tbl_a46c4afa9069", 40, 2): ("5.67E-C6", "5.07E-06"),
    ("tbl_a46c4afa9069", 45, 1): ("0.0C1", "0.000"),
    ("tbl_a46c4afa9069", 45, 3): ("0.0C1", "0.000"),
    ("tbl_a46c4afa9069", 46, 3): ("-0.0C1", "-0.000"),
    ("tbl_a46c4afa9069", 47, 1): ("-0.0C1", "-0.000"),
    ("tbl_a46c4afa9069", 47, 3): ("-0.0C1", "-0.000"),
    ("tbl_a46c4afa9069", 48, 6): ("0.0C1", "0.000"),
    ("tbl_a46c4afa9069", 49, 1): ("-0.0C1", "-0.000"),
    ("tbl_a46c4afa9069", 49, 3): ("-0.0C1", "-0.000"),
    # PDF page 52, Table 6 page 1 (tbl_376d0c04b097)
    ("tbl_376d0c04b097", 1, 5): ("6.15E 05", "6.75E 05"),
    ("tbl_376d0c04b097", 2, 2): ("1.5CE 00", "1.90E 00"),
    ("tbl_376d0c04b097", 3, 2): ("2.64E 00", "2.66E 00"),
    ("tbl_376d0c04b097", 7, 4): ("1.50E 00", "1.90E 00"),
    ("tbl_376d0c04b097", 11, 2): ("6.E5E 00", "6.85E 00"),
    ("tbl_376d0c04b097", 15, 2): ("5.25E 00", "5.29E 00"),
    ("tbl_376d0c04b097", 15, 7): ("1.64E-C1", "1.84E-01"),
    ("tbl_376d0c04b097", 17, 8): ("C.564", "0.584"),
    ("tbl_376d0c04b097", 24, 3): ("1.352", "1.392"),
    ("tbl_376d0c04b097", 27, 0): ("1.CO E 01", "1.05E 01"),
    ("tbl_376d0c04b097", 27, 6): ("19.636", "18.636"),
    ("tbl_376d0c04b097", 28, 0): ("1.CO E 01", "1.10E 01"),
    ("tbl_376d0c04b097", 33, 8): ("2.599", "2.589"),
    ("tbl_376d0c04b097", 34, 1): ("0.563", "0.963"),
    ("tbl_376d0c04b097", 36, 7): ("4.19E-C1", "4.10E-01"),
    ("tbl_376d0c04b097", 38, 6): ("11.450", "11.460"),
    ("tbl_376d0c04b097", 41, 2): ("2.C5E 00", "2.09E 00"),
    ("tbl_376d0c04b097", 42, 5): ("2.04E 06", "2.06E 06"),
    ("tbl_376d0c04b097", 43, 2): ("1.75E CC", "1.79E 00"),
    ("tbl_376d0c04b097", 47, 6): ("11.618", "11.619"),
    ("tbl_376d0c04b097", 49, 2): ("1.C5E 00", "1.09E 00"),
    # PDF page 53, Table 6 page 2 (tbl_28efb140cb95)
    ("tbl_28efb140cb95", 0, 0): ("2.5CE 01", "2.90E 01"),
    ("tbl_28efb140cb95", 0, 5): ("1.36E 06", "1.76E 06"),
    ("tbl_28efb140cb95", 0, 7): ("7.61E-C1", "7.81E-01"),
    ("tbl_28efb140cb95", 1, 5): ("1.83E 06", "1.80E 06"),
    ("tbl_28efb140cb95", 2, 5): ("1.93E 06", "1.83E 06"),
    ("tbl_28efb140cb95", 2, 8): ("5.685", "6.689"),
    ("tbl_28efb140cb95", 4, 5): ("1.67E 06", "1.87E 06"),
    ("tbl_28efb140cb95", 4, 7): ("9.70E-C1", "8.70E-01"),
    ("tbl_28efb140cb95", 5, 7): ("9.79E-C1", "9.70E-01"),
    ("tbl_28efb140cb95", 6, 7): ("1.C0E CC", "1.05E 00"),
    ("tbl_28efb140cb95", 12, 7): ("9.89E-C1", "9.85E-01"),
    ("tbl_28efb140cb95", 13, 8): ("8.658", "8.659"),
    ("tbl_28efb140cb95", 15, 7): ("8.54E-C1", "8.94E-01"),
    ("tbl_28efb140cb95", 20, 2): ("3.65E-01", "3.85E-01"),
    ("tbl_28efb140cb95", 24, 2): ("3.15E-01", "3.19E-01"),
    ("tbl_28efb140cb95", 25, 1): ("0.667", "0.687"),
    ("tbl_28efb140cb95", 26, 2): ("2.53E-01", "2.93E-01"),
    ("tbl_28efb140cb95", 29, 2): ("2.65E-01", "2.69E-01"),
    ("tbl_28efb140cb95", 39, 7): ("1.57E-C1", "1.97E-01"),
    ("tbl_28efb140cb95", 43, 2): ("1.14E-01", "1.16E-01"),
    ("tbl_28efb140cb95", 45, 2): ("1.C5E-01", "1.09E-01"),
    ("tbl_28efb140cb95", 47, 2): ("9.51E-02", "9.91E-02"),
    ("tbl_28efb140cb95", 49, 0): ("8.50E 01", "8.90E 01"),
    # PDF page 54, Table 6 page 3 (tbl_5ab449433de5)
    ("tbl_5ab449433de5", 0, 3): ("0.1C5", "0.109"),
    ("tbl_5ab449433de5", 3, 0): ("1.CGE 02", "1.00E 02"),
    ("tbl_5ab449433de5", 10, 7): ("6.$9E-C3", "6.99E-03"),
    ("tbl_5ab449433de5", 16, 1): ("0.021", "0.031"),
    ("tbl_5ab449433de5", 17, 8): ("15.578", "15.978"),
    ("tbl_5ab449433de5", 20, 2): ("1.C5E-02", "1.09E-02"),
    ("tbl_5ab449433de5", 23, 3): ("C.C0S", "0.009"),
    # PDF page 60, Table 9 page 1 (tbl_7cb0a953178f)
    ("tbl_7cb0a953178f", 2, 0): ("9.00E 00", "8.00E 00"),
    ("tbl_7cb0a953178f", 2, 8): ("0.072", "0.022"),
    ("tbl_7cb0a953178f", 3, 0): ("9.50E 00", "8.50E 00"),
    ("tbl_7cb0a953178f", 3, 2): ("2.6CE-01", "2.06E-01"),
    ("tbl_7cb0a953178f", 4, 1): ("5.269", "5.268"),
    ("tbl_7cb0a953178f", 6, 2): ("3.68E 00", "3.08E 00"),
    ("tbl_7cb0a953178f", 7, 6): ("23.333", "23.833"),
    ("tbl_7cb0a953178f", 10, 2): ("4.26E 00", "4.36E 00"),
    ("tbl_7cb0a953178f", 12, 3): ("1.645", "1.649"),
    ("tbl_7cb0a953178f", 12, 8): ("3.949", "3.849"),
    ("tbl_7cb0a953178f", 13, 8): ("4.393", "4.399"),
    ("tbl_7cb0a953178f", 15, 1): ("0.553", "0.558"),
    ("tbl_7cb0a953178f", 16, 2): ("2.54E 00", "2.94E 00"),
    ("tbl_7cb0a953178f", 18, 2): ("2.77E 00", "2.07E 00"),
    ("tbl_7cb0a953178f", 18, 6): ("19.337", "19.387"),
    ("tbl_7cb0a953178f", 20, 8): ("9.057", "9.067"),
    ("tbl_7cb0a953178f", 21, 3): ("0.886", "0.880"),
    ("tbl_7cb0a953178f", 24, 6): ("14.553", "14.653"),
    ("tbl_7cb0a953178f", 25, 2): ("9.57E-01", "9.97E-01"),
    ("tbl_7cb0a953178f", 26, 8): ("11.669", "11.609"),
    ("tbl_7cb0a953178f", 28, 2): ("5.66E-01", "5.06E-01"),
    ("tbl_7cb0a953178f", 29, 8): ("13.593", "13.583"),
    ("tbl_7cb0a953178f", 31, 1): ("0.667", "0.697"),
    ("tbl_7cb0a953178f", 31, 7): ("4.67E-01", "4.57E-01"),
    ("tbl_7cb0a953178f", 33, 2): ("1.64E-01", "1.54E-01"),
    ("tbl_7cb0a953178f", 35, 2): ("9.56E-02", "9.96E-02"),
    ("tbl_7cb0a953178f", 36, 8): ("20.260", "20.240"),
    ("tbl_7cb0a953178f", 37, 2): ("6.5CE-02", "6.90E-02"),
    ("tbl_7cb0a953178f", 38, 7): ("7.58E-02", "7.98E-02"),
    ("tbl_7cb0a953178f", 38, 8): ("21.079", "21.009"),
    ("tbl_7cb0a953178f", 41, 2): ("5.55E-02", "5.95E-02"),
    ("tbl_7cb0a953178f", 44, 8): ("21.497", "21.407"),
    ("tbl_7cb0a953178f", 45, 2): ("5.SCE-02", "5.90E-02"),
    ("tbl_7cb0a953178f", 48, 2): ("6.58E-02", "6.98E-02"),
    ("tbl_7cb0a953178f", 49, 8): ("21.793", "21.703"),
    # PDF page 61, Table 9 page 2 (tbl_73348bbbe4f6)
    ("tbl_73348bbbe4f6", 0, 0): ("7.50E 01", "7.90E 01"),
    ("tbl_73348bbbe4f6", 0, 2): ("7.66E-02", "7.69E-02"),
    ("tbl_73348bbbe4f6", 0, 7): ("3.80E-02", "8.80E-02"),
    ("tbl_73348bbbe4f6", 1, 0): ("7.55E 01", "7.95E 01"),
    ("tbl_73348bbbe4f6", 1, 6): ("0.070", "0.079"),
    ("tbl_73348bbbe4f6", 1, 7): ("5.51E-02", "9.51E-02"),
    ("tbl_73348bbbe4f6", 2, 6): ("0.070", "0.090"),
    ("tbl_73348bbbe4f6", 2, 7): ("4.57E-02", "9.57E-02"),
    ("tbl_73348bbbe4f6", 2, 8): ("21.518", "21.919"),
    ("tbl_73348bbbe4f6", 3, 1): ("0.915", "0.916"),
    ("tbl_73348bbbe4f6", 3, 7): ("5.77E-02", "9.77E-02"),
    ("tbl_73348bbbe4f6", 5, 0): ("3.20E 01", "8.20E 01"),
    ("tbl_73348bbbe4f6", 7, 2): ("6.55E-02", "6.59E-02"),
    ("tbl_73348bbbe4f6", 10, 8): ("22.857", "22.807"),
    ("tbl_73348bbbe4f6", 11, 8): ("22.974", "22.924"),
    ("tbl_73348bbbe4f6", 12, 6): ("0.059", "0.056"),
    ("tbl_73348bbbe4f6", 17, 2): ("5.62E-02", "6.62E-02"),
    ("tbl_73348bbbe4f6", 20, 3): ("0.956", "0.966"),
    ("tbl_73348bbbe4f6", 20, 7): ("6.03E-02", "6.08E-02"),
    ("tbl_73348bbbe4f6", 24, 7): ("3.98E-02", "3.96E-02"),
    ("tbl_73348bbbe4f6", 29, 7): ("2.57E-02", "2.97E-02"),
    ("tbl_73348bbbe4f6", 37, 2): ("3.65E-03", "3.60E-03"),
    ("tbl_73348bbbe4f6", 47, 7): ("6.89E-05", "6.80E-05"),
    ("tbl_73348bbbe4f6", 48, 2): ("5.53E-05", "5.93E-05"),
}

# (library table id, data row) -> the checks that fail beyond TOLERANCE on the page
# itself (the cells read again on the image): see the module docstring.
INCONSISTENT: Final[dict[tuple[str, int], frozenset[str]]] = {
    ("tbl_5eadb4abaadf", 22): frozenset({"eps2"}),
    ("tbl_5eadb4abaadf", 30): frozenset({"loss"}),
    ("tbl_5eadb4abaadf", 34): frozenset({"loss"}),
    ("tbl_5eadb4abaadf", 49): frozenset({"loss", "mu"}),
    ("tbl_e9bd35b2f535", 3): frozenset({"loss"}),
    ("tbl_e9bd35b2f535", 5): frozenset({"loss", "mu"}),
    ("tbl_e9bd35b2f535", 6): frozenset({"loss", "mu"}),
    ("tbl_e9bd35b2f535", 7): frozenset({"loss", "mu"}),
    ("tbl_e9bd35b2f535", 8): frozenset({"loss", "mu"}),
    ("tbl_e9bd35b2f535", 9): frozenset({"loss", "mu"}),
    ("tbl_7ad069048df5", 30): frozenset({"mu"}),
    ("tbl_52ae5a332285", 0): frozenset({"mu"}),
    ("tbl_5954ad9a5d51", 2): frozenset({"eps2"}),
    ("tbl_5954ad9a5d51", 19): frozenset({"R"}),
    ("tbl_44c1d033263c", 26): frozenset({"mu"}),
    ("tbl_e95227627cb4", 1): frozenset({"eps1", "eps2", "loss", "mu"}),
    ("tbl_e95227627cb4", 4): frozenset({"loss", "mu"}),
    ("tbl_e95227627cb4", 9): frozenset({"loss"}),
    ("tbl_7cb0a953178f", 0): frozenset({"eps2", "loss", "mu"}),
    ("tbl_7cb0a953178f", 47): frozenset({"loss"}),
    ("tbl_73348bbbe4f6", 2): frozenset({"loss"}),
    ("tbl_73348bbbe4f6", 36): frozenset({"mu"}),
    ("tbl_73348bbbe4f6", 37): frozenset({"mu"}),
    ("tbl_73348bbbe4f6", 38): frozenset({"mu"}),
    ("tbl_73348bbbe4f6", 40): frozenset({"mu"}),
}

_ZEROS: Final = str.maketrans(dict.fromkeys("COQDUco", "0"))
_NUMBER: Final = re.compile(r"^(-?)(\d+)\.(\d*)(?:E([ -])(\d\d))?$")

type Interval = tuple[float, float]


def cell_text(raw: str) -> str:
    """Return a library cell as the printout prints it (zeros read as letters)."""
    text = raw.strip().translate(_ZEROS).replace(",", ".")
    text = re.sub(r"(\d)F([ -]?\s*\d)", r"\1E\2", text)
    text = re.sub(r"\s+E", "E", text)
    text = re.sub(r"E\s+(\d)", r"E \1", text)
    return re.sub(r"E-\s+", "E-", text)


def parse(text: str) -> tuple[float, float] | None:
    """Return a printed number and half a unit of its last digit, or None if blank."""
    if not text or "?" in text:
        return None
    match = _NUMBER.match(text)
    if match is None:
        msg = f"not a printed number: {text!r}"
        raise ValueError(msg)
    sign, whole, decimals, exponent_sign, exponent = match.groups()
    power = 0
    if exponent is not None:
        power = int(exponent) * (-1 if exponent_sign == "-" else 1)
    value = float(f"{sign}{whole}.{decimals or '0'}") * 10.0**power
    return value, 0.5 * 10.0 ** (power - len(decimals))


def stored(text: str) -> str:
    """Return a printed number as stored: ``E 05`` written ``E+05``."""
    return text.replace("E ", "E+")


def _span(function: Callable[..., float], *intervals: tuple[float, float]) -> Interval:
    """Return the range of ``function`` over the rounding intervals of its arguments."""
    values = [
        function(*point)
        for point in itertools.product(*((v - h, v, v + h) for v, h in intervals))
    ]
    return min(values), max(values)


def _excess(printed: tuple[float, float], span: Interval) -> float:
    """Return how far a printed number lies outside ``span``, relative to it."""
    value, half = printed
    low, high = span
    if value + half < low:
        return (low - value - half) / abs(low)
    if value - half > high:
        return (value - half - high) / abs(high)
    return 0.0


def failures(cells: Sequence[str], complement: bool) -> frozenset[str]:
    """Return the redundancy checks a row fails beyond ``TOLERANCE``.

    ``complement`` is true on the pages that print 1 - epsilon_1 and 1 - n.
    """
    energy, eps1, eps2, n, k, mu, reflect, loss = (parse(c) for c in cells[:NEFF])
    if energy is None or k is None or mu is None:
        msg = f"row without energy, k or mu: {cells}"
        raise ValueError(msg)
    excess = {"mu": _excess(mu, _span(lambda e, b: 2 * e * b / HBARC_EV_CM, energy, k))}
    if n is not None and eps1 is not None and eps2 is not None:
        if complement:
            eps1, n = (1 - eps1[0], eps1[1]), (1 - n[0], n[1])
        excess["eps1"] = _excess(eps1, _span(lambda a, b: a * a - b * b, n, k))
        excess["eps2"] = _excess(eps2, _span(lambda a, b: 2 * a * b, n, k))
        if reflect is not None and reflect[0] > 0.05:
            low, high = _span(
                lambda a, b: 100 * ((a - 1) ** 2 + b * b) / ((a + 1) ** 2 + b * b), n, k
            )
            if n[0] - n[1] <= 1 <= n[0] + n[1]:
                low = 0.0
            excess["R"] = _excess(reflect, (low, high))
        if loss is not None:
            excess["loss"] = _excess(
                loss, _span(lambda a, b: b / (a * a + b * b), eps1, eps2)
            )
    return frozenset(name for name, value in excess.items() if value > TOLERANCE)


def page_rows(table_id: str) -> tuple[bool, list[list[str]]]:
    """Return whether a page prints 1 - epsilon_1 and 1 - n, and its rows as printed."""
    grid = json_table(DIRECTORY, table_id)
    header = next(i for i, row in enumerate(grid) if row[0].startswith("ENE"))
    complement = grid[header][EPS1].startswith("1.-")
    rows = []
    for index, raw in enumerate(grid[header + 1 :]):
        if any(cell.strip() for cell in raw[NEFF + 1 :]):
            msg = f"{table_id}: row {index} has cells beyond N-EFF: {raw}"
            raise ValueError(msg)
        cells = [cell_text(cell) for cell in raw[: NEFF + 1]]
        for column in range(NEFF + 1):
            if (table_id, index, column) in CORRECTIONS:
                library, page = CORRECTIONS[table_id, index, column]
                if raw[column] != library:
                    msg = f"{table_id}: row {index} col {column} is {raw[column]!r}"
                    raise ValueError(msg)
                cells[column] = page
        found = failures(cells, complement)
        if found != INCONSISTENT.get((table_id, index), frozenset()):
            msg = f"{table_id}: row {index} fails {sorted(found)}: {cells}"
            raise ValueError(msg)
        rows.append(cells)
    return complement, rows


def extract() -> dict[str, list[Point]]:
    """Return n (or delta), k, epsilon_1, epsilon_2 and mu of each table."""
    files: dict[str, list[Point]] = {}
    for material, suffix, pages in TABLES.values():
        table: dict[str, list[Point]] = {kind: [] for kind in DATASETS}
        previous_energy, previous_neff = 0.0, 0.0
        for table_id in pages:
            complement, rows = page_rows(table_id)
            for cells in rows:
                energy = parse(cells[ENERGY])
                if energy is None or energy[0] <= previous_energy:
                    msg = f"{table_id}: energies out of order at {cells}"
                    raise ValueError(msg)
                previous_energy = energy[0]
                neff = parse(cells[NEFF])
                if neff is not None:
                    if neff[0] < previous_neff:
                        msg = f"{table_id}: N-EFF decreases at {cells}"
                        raise ValueError(msg)
                    previous_neff = neff[0]
                row = stored(cells[ENERGY])
                printed = " | ".join(cells)
                for kind, column in DATASETS.items():
                    if kind in (PLAIN_ONLY if complement else COMPLEMENT_ONLY):
                        continue
                    if not cells[column]:
                        continue
                    table[kind].append(
                        Point(row, row, stored(cells[column]), "", printed)
                    )
        for kind, points in table.items():
            if points:  # Table 9 prints no "1.-N" page, so it has no delta
                files[f"{material}/hagemann_1974__{kind}_{suffix}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()

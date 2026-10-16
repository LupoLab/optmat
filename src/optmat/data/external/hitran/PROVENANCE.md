# HITRAN: infrared band strengths of CH4, N2O, SF6, CO2, O3 and CO

Derived from HITRANonline (https://hitran.org), the current edition being HITRAN2024:
I. E. Gordon, L. S. Rothman, R. J. Hargreaves et al., "The HITRAN2024 molecular
spectroscopic database", J. Quant. Spectrosc. Radiat. Transfer 353, 109807 (2026). HITRAN's
citation policy (https://hitran.org/citepolicy/) asks for that citation and encourages
citing the original sources of the data points: the HITRAN reference codes of the
intensities are kept in `intensity_refs` (the second two-character field of each line's
reference block), resolvable in HITRANonline's source list.

Queried on 2026-09-25 through the HITRAN API version 2 with a HITRAN API key (not
kept in the repository), by `tools/hitran/band_strengths.py`:

| File | Source | Query | Returned |
|---|---|---|---|
| `CH4_bands.csv` | line list, all isotopologues (global ids 32-35) | 0-11000 cm^-1 | 456922 lines (10:25:12 UTC) |
| `N2O_bands.csv` | line list, all isotopologues (global ids 21-25) | 0-11000 cm^-1 | 514508 lines (10:24:30 UTC) |
| `CO2_bands.csv` | line list, all isotopologues | 0-11000 cm^-1 | 540624 lines (queried later on 2026-09-25, when CO2 was taken up) |
| `O3_bands.csv` | line list, all isotopologues | 0-11000 cm^-1 | 487446 lines, 286 bands (queried on 2026-09-25, when O3 was taken up) |
| `CO_bands.csv` | line list, all isotopologues | 0-11000 cm^-1 | 5210 lines, 83 bands (queried on 2026-09-25, when CO was taken up) |
| `SF6_bands.csv` | cross-sections, Ref. [1164] (298.1 K, 760 Torr N2, 0.112 cm^-1, 560-6500 cm^-1) | `SF6_298.1K-760.0Torr_560.0-6500.0_0.11_N2_220_43.xsc` | integrated in three windows |

SF6 has no line list in HITRAN; its cross-sections are measured, and the one file of Ref.
[1164] covers every band. Its 800-1075 cm^-1 window (nu3) integrates to 1.9619e-16
cm/molecule against 1.9811e-16 from the independent Ref. [1638] file
(`SF6_296.1K_745.7Torr_800.0-1075.0_0.040_air_140.xsc`): -0.97 %.

Only the derived band totals are committed: per band (grouped by isotopologue and the
global upper and lower quanta) the number of lines, the intensity-weighted wavenumber and
range, `S` = sum of the line intensities (cm/molecule at 296 K, per molecule of the
natural isotopic mixture, HITRAN's definition, which includes the lower-state population
and stimulated emission) and sum of `S / nu^2`. The line lists and cross-sections
themselves are fetched on demand into `~/.cache/optmat/hitran` and not
redistributed (hitran.org states no redistribution terms; this keeps to the citation
policy and the API's use).

Checked against Shelton's phase-matching density ratios (TC-104): for SF6 the two
wavelengths are consistent with HITRAN's strengths scaled by 1.007, for CH4 by 0.88
(+-0.08).

## Derived windows, band rules and bound (OEP 3)

`tools/hitran/windows.py` derives, from the same cached line lists of CO, CO2, CH4 and
N2O and the parent fits' electronic terms, four tables per gas:

| File | Holds |
|---|---|
| `<gas>_band_rules.csv` | Every band system (a band, or its part between two windows; all pure rotational lines as one system; the weak bands inside each window folded onto two poles outside it), its lines' totals M_0 and M_-2, and the nodes and strengths of its Gaussian rule |
| `<gas>_windows.csv` | The qualified intervals from 12.5 um to the parent fit's shortest wavelength, and the gaps between them, each with the condition that failed there; the absorption at the edges, the largest single-line Lorentz peak and far-wing absorption, and the bound over the tolerance |
| `<gas>_bound.csv` | The representation bound, the largest \|D_m[model] - D_m[reference]\| for m = 0, 2, 3 in atomic units, per geometric cell of the windows |
| `<gas>_absorption.csv` | The pure gas's absorption at 101325 Pa and 296 K, smoothed by the Gaussian, its maximum per 1 cm^-1 where it is at least 1e-4 m^-1 (for the windows of a mixture) |

Each file begins with `#` lines naming the tool's version and commit, the SHA-256 of the
cached line list and of its query record, the parent fit's configuration digest and
every setting. Only these derived values are committed; the line lists are not
redistributed. The method, the settings and the checks are those of OEP 3, sections 4
and 7, and the register records each finding.

## The fits' gaps (OEP 3, Q13)

`tools/hitran/gap_bound.py` measures, for the stored fit of CO, CO2, CH4 and N2O, the
same quantities in the fit's gaps, from 12.5 um to the top of the absorption grid
(11 050 cm^-1), where a mixture holding the gas as a trace may still be qualified:

| File | Holds |
|---|---|
| `<gas>_gap_bound.csv` | Per geometric cell (ratio 1.001) of the gaps: the largest \|D_m[fit] - D_m[reference]\| for m = 0, 2, 3 (atomic units) and the largest smoothed absorption of the pure gas at 101325 Pa and 296 K (m^-1), on the 0.25 cm^-1 grid and the half-step offset grid; `pole = 1` where a pole of the fit lies in the cell, whose bound is then unbounded |

Each begins with the same `#` lines as the tables above, and the stored fit's SHA-256.
The air model reads them (`optmat.gas.air`; TC-809).

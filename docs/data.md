# The curated data

This reference describes source records, datasets, corrections and the conversion
from printed values to physical quantities. Start with {doc}`/introduction/curation`
for an overview, or browse {doc}`/the-data/sources` to find a paper. The
{doc}`/the-data/quantities` defines the fields used by each quantity family.

## Sources and datasets

A **source** is one paper, compilation, theory table or database. Its record,
`src/optmat/data/sources/<source>.toml`, holds:

| Key | Meaning |
|---|---|
| `citation` | The citation as it should be printed. |
| `doi` | The DOI, or omitted. |
| `library` | The key of the paper in the curators' private paper library, from which the extraction scripts read its tables (`$PAPER_LIBRARY/<library>/tables/`); omitted if the work is not in the library. Nothing but the extraction needs it. |
| `identity` | How firmly the library copy's bibliographic identity was established (`VALIDATED`, `VALIDATED_WITH_WARNINGS`, `UNVERIFIED`). |
| `kind` | `measurement`, `theory`, `compilation`, `fit` or `database`. |
| `technique` | How the values were obtained: instrument, method, source of light, conditions. |
| `notes` | Everything else a reader needs: how the values were reduced, what the paper says about accuracy (quoted, with the page), what is odd. |
| `[[datasets]]` | One entry per dataset (below). |
| `[[errata]]` | One entry per corrected printed value (below). |
| `[[formulas]]` | Recorded equations, coefficients, ranges, conventions, checks and any formula corrections. |

The source id is the record's file stem, `<first author>_<year>` in lower-case ASCII
(`bideau_mehu_1981`), with a letter suffix if one author has two papers in a year.

A **dataset** is a set of values of one quantity for one material from one source:

| Key | Meaning |
|---|---|
| `id` | Unique within the source. |
| `material` | The material's id in the registry (`He`, `N2`, `SiO2`, `CaF2`, ...; see *Materials* below). |
| `file` | The extracted data file, `<material>/<source>__<id>.csv`. It lives in `data/values/<family>/`, the family of the dataset's quantity (below), so a file is always where its quantity puts it. Omitted when the values are given inline as `points`. |
| `points` | Inline values `[[row, x, value], ...]` (or with a fourth element, the per-point uncertainty, and then the extra columns in their declared order), for values that exist only in the text. The page goes in `notes`. |
| `table` | The library table id (`tbl_...`); for values read from a figure, the library figure id (`fig_...`); `"data/<file>"` for a data file attached byte for byte to the paper in the library (a supporting-information file, an author's data file, a database's table), with its origin recorded there; or `"text"`. |
| `digitised` | For values read from a figure, where the source prints no number (TC-699): `"vector"`, read from the figure's vector drawing in the PDF (marker and line coordinates, the axes calibrated from the drawing's ticks), or `"raster"`, read from its image. Such a dataset states its uncertainty (the digitising uncertainty with the source's), and is `use = "check"`; raster values are never `primary`. |
| `quantity` | What the values are (see {doc}`/the-data/quantities`). |
| `unit` | The unit after `scale`: `"1"` for dimensionless (the default), or a unit named by the quantity. |
| `scale` | Printed values times `scale` are in `unit`: `1e-6` for a column of `(n - 1) x 10^6`. |
| `abscissa` | What the `x` column is (the abscissae in {doc}`/the-data/quantities`): `wavelength` (the default), `order`, `temperature`, `pressure`, `density`, `intensity`, `duration` (a pulse length), `field` (an electric field), `label` or `none`. |
| `abscissa_unit` | The unit of the `x` column. For a wavelength: `m`, `nm`, `um`, `angstrom`, a wavenumber `cm-1` or `um-1`, or a photon energy `eV`, `hartree` or `rydberg` (R∞hc; a source in reduced-mass units says so in `parameters`). The cell `inf` (or an energy of 0) is the static limit, and `not stated` a wavelength the source does not give. Omitted for a `label` or `none`. |
| `wavelength_medium` | `vacuum` (the default) or `standard_air` (15 °C, 101 325 Pa; converted with Edlén 1966), for wavelengths and wavenumbers (a wavenumber in air is the inverse of the air wavelength). |
| `columns` | Extra per-row columns and their units, `{ x_upper = "eV" }`: `x_upper` is the upper edge of an interval, in the abscissa's unit; any other column is a per-row condition (a temperature, a probe wavelength) with its own unit. A quantity may require some. |
| `state` | The state the values refer to or are reduced to: `"1atm_0C"`, `"1bar_0C"`, `"1atm_15C"`, `"1atm_20C"`, or a table `{ pressure_pa = ..., temperature_k = ... }` (also `pressure_atm`, `temperature_c`). Required for a gas's quantities that depend on density; a solid has no state (its temperature, where stated, is `parameters.temperature_k`). |
| `density` | How the source related its values to the state: `real` (a real-gas reduction: virial or compressibility tables, van der Waals), `ideal` (the ideal-gas law), or `unknown` (the default; the only value for a solid). |
| `[datasets.uncertainty]` | `relative`, `absolute` (printed unit), `column = true` for a per-point `sigma` column, `level` (`"1 sigma"`, `"probable error"`, `"90 % confidence"`, `"95 % confidence"`, `"99 % confidence"`, `"2 sigma"`, `"3 sigma"`, `"rms"`, `"bound"`, or `"not stated"`), `systematic_relative` (a scale error common to all points, one sigma), `relative_rows = { "<row>" = ... }` for rows with their own relative uncertainty. |
| `use` | `primary` (original data eligible for selection by a fit), `reproduced` (a copy of another source's values), `check` (compared, not fitted: a paper's own fit column, say), or `excluded` (kept for the record; `notes` says why). |
| `reproduced_from` | The source ids of reproduced values, or their citations when the original is not in the library. Required for `use = "reproduced"`. |
| `[datasets.parameters]` | Quantity-specific metadata, listed with the vocabulary; `beyond_resonance = true` marks values on the far side of a resonance, where $n - 1$ may be negative. |
| `excluded_rows` | `{ "<row>" = "<reason>" }`: rows of the file left out of the dataset (the paper's own exclusions, points it marks extrapolated). |
| `notes` | Anything specific to this dataset. |

A use designation records the intended role of a dataset; the model configuration
determines which datasets actually enter a fit. A source may contain only formulas
or notes, with no datasets.

A reproduced dataset is never fitted beside its original. When the original is not in
the library (Chashchina's vacuum-ultraviolet line shifts, for example, survive only in
Bideau-Mehu's tables), the reproduction is the best copy available and may be fitted.

## Materials

Every dataset names its material by an id from the registry, `data/materials.toml`
({mod}`optmat.materials`). A material is registered, with its identity, before data for
it are recorded. An entry gives:

- `name`, and the `formula` where it has one;
- `kind`: `gas`, `ion` (a gas-phase cation, whose data never enter the neutral gas's
  models), `glass`, `crystal` or `metal` (TC-543);
- `optics`: `isotropic` (the default), `uniaxial` or `biaxial`, and for an anisotropic
  crystal the `axes` a dataset's `parameters.axis` may name (`o` and `e`; `n_p`, `n_m`
  and `n_g`, or `x`, `y` and `z`, for a biaxial crystal);
- `aliases`, other names {func}`optmat.materials.as_material` accepts ("fused silica").
  A record always uses the id.

The id is the directory of the material's files. A glass's grade (Corning 7980,
Suprasil 1, Lithosil Q2) is a dataset parameter, not a separate material: grades
differ in their impurities and thermal history, which the data record and the models
state. Crystalline quartz (`SiO2_quartz`) is a different material from silica glass
(`SiO2`). A gas's state and density convention have no meaning for a solid, and the
loader refuses them there (TC-318).

## Data files

The curators' paper library is prepared with
[paperextract](https://github.com/jtravs/paperextract), which keeps structured tables,
equations and figures linked to their source pages. Source-specific extraction scripts
and checks against the original tables and published numerical results establish the
values used here; {doc}`about` describes the extraction and verification process.

A data file, `data/values/<family>/<material>/<source>__<id>.csv`, has the columns `row`,
`x`, `value`, `sigma` and `printed`, then the dataset's extra columns in their declared
order. The numbers are **as printed**, with only a column factor such as
`\times 10^{-5}` removed (the dataset's `scale` accounts for it); `printed` keeps the raw
cells of the table row, for audit. Files are written only by the source's extraction
script, `tools/extract/<source>.py`, from the library's `tables/<id>.csv`: never by
hand. The script names a file as the record does (`<material>/<source>__<id>.csv`) and the
helper puts it in its family's directory. A test re-runs every script against the
library and compares (`poe test-library`); another checks that every file belongs to
exactly one dataset and lies in its quantity's family.

## Errata

A value the paper prints wrongly is corrected by an erratum in the source record, never
by editing the data file:

```toml
[[errata]]
dataset = "fit_Ne"
row = "253.7"
column = "value"          # x, value, sigma or an extra column
printed = "69.60"         # checked against the extracted value on loading
corrected = "68.60"
reason = "..."
evidence = "..."
```

`printed` must equal what the file holds, or loading fails: an erratum cannot silently
outlive a change in the extraction. A misprinted column factor is corrected by an
erratum of the dataset's scale, `row = ""` and `column = "scale"`, with the scale as
recorded in `printed` (Börzsönyi et al.'s 10^-19 for 10^-20).

A misprinted coefficient of a published formula is corrected in the formula's entry:

```toml
[[formulas.corrections]]
term = "C"                # a named coefficient of `coefficients`
printed = "6.043e-31"     # checked against the recorded coefficient
corrected = "6.043e-32"
reason = "..."
evidence = "..."
```

A formula with a `form` names a coefficient list and an `index` in it instead
(`term = "B"`, `index = 0`). Structural problems of an extracted table (a caption
attached to the wrong table, rows misaligned under a brace, markers lost) are handled in
the extraction script, with a comment saying what the table image shows.

## Quantities

The complete {doc}`quantity vocabulary </the-data/quantities>` lists families, units,
abscissae, required metadata and third-order conventions.

(third-order-quantities)=
See also [third-order quantities](the-data/quantities.md#third-order-quantities).

## From printed values to polarizability

`optmat.gas.reduction.polarizability(dataset)` converts a refractivity at the dataset's state with the
number density the source meant (real-gas unless the source used the ideal gas),
through Lorentz-Lorenz; a polarizability in atomic units with $4\pi\varepsilon_0 a_0^3$;
a molar refractivity with $3\varepsilon_0/N_A$. `Dataset.standard_uncertainty()` turns
the stated uncertainty into one standard deviation (a probable error is 0.6745 sigma, a
90 % interval 1.645 sigma, a stated bound a uniform distribution's half-width).

## Infrared bands (HITRAN)

The infrared-active molecules have vibrational bands whose strengths add to the
polarizability: a large vibrational static polarizability for SF6 (13.9 a.u.) and
N2O (1.19 a.u.), a small one for CH4 (0.24 a.u.), and at optical frequencies a small
negative term that matters in differences of polarizabilities, such as phase-matching
conditions. The band strengths come from HITRAN
(`data/external/hitran`, written by `tools/hitran/band_strengths.py` with a HITRAN
API key; the line lists themselves are not kept): line lists summed by band
for CH4 and N2O, measured cross-sections integrated for SF6. {mod}`optmat.gas.infrared`
provides both a single-oscillator band representation and the band rules described
below. The older band representation keeps each band’s total strength and static
contribution and excludes pure rotation (TC-156). The current CO, CO2, CH4 and N2O
fits use the band rules; CO and N2O carry rotation as a separate system, whose
static contribution is the orientation polarizability ({doc}`conventions`).

Shelton's phase-matching ratios at 1319 and 1064 nm test them: the infrared share of
SF6's dispersion differs by a factor of 2.4 between the two wavelengths, and the two are
consistent with HITRAN's strengths to 1 % (CH4, whose share is 6 %, to 12 +- 8 %). The same
test shows that He et al.'s formulas overstate the electronic dispersion of SF6 by 23 %
and of CH4 by 3 % (TC-104).

A few effective poles in place of a molecule's bands, for print, are its compact model
(`compact/<gas>`; {doc}`fitting`, *Compact infrared terms*).

**Band rules and windows** ({doc}`OEP 3 <oeps/archive/oep-0003-infrared-windows>`).
`tools/hitran/windows.py` derives, from the same line lists of CO, CO2, CH4 and N2O,
the tables used by the molecular fits (`<gas>_band_rules.csv`,
`<gas>_windows.csv`, `<gas>_bound.csv`, `<gas>_absorption.csv`; see `PROVENANCE.md`
there). Each band, split at the windows, becomes a Gaussian quadrature rule of up to
eight poles that keeps its lines' moments (Golub and Welsch); the weak bands inside a
window are folded onto two poles outside it, or kept, which splits the window; the
pure rotational lines are one more system. The reference is the gas's refractivity
over every line, smoothed by a Gaussian of 1 cm^-1 (a 5.3 ps gate). A window is where
the pure gas's smoothed absorption at 1 atm and 296 K is below 0.1 m^-1 and, beyond
the fit's optical range, the model is within the tolerance of that reference: a tenth
of the fit's statistical sigma in the polarizability, and 1 % of the larger of the
total and the electronic GVD. Inside the optical range only the absorption splits it,
and the bound is reported (TC-749). The loaders are
{func}`~optmat.gas.infrared.band_systems`, {func}`~optmat.gas.infrared.infrared_windows`
and {func}`~optmat.gas.infrared.bound_cells`. The molecular fits are made with them
(TC-797). {func}`~optmat.gas.infrared.gap_cells` loads each fit's error and absorption in
its gaps (`<gas>_gap_bound.csv`), with which a mixture's intervals widen across a trace
constituent's bands (TC-809).

The band rules, windows and gap tables are what
[OEP 3](oeps/archive/oep-0003-infrared-windows.md) proposed and implemented: its theory,
decisions and phase reports are the record.

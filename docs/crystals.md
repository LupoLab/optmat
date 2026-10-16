# The window crystals

CaF2, MgF2, BaF2 and sapphire (`Al2O3`) are the window and substrate crystals of
ultrafast optics, and N-BK7 its commonest glass. Each recommended model is the Sellmeier
formula of its source, made absolute (TC-545, TC-546, TC-651). This page explains the
reference-medium corrections, source agreement and temperature limits. Read
{doc}`solid` for the common interface and {doc}`/quickstart` for a first calculation.

```python
from optmat.solid.models.registry import default_model, default_reason, models

caf2 = default_model("CaF2")  # daimon_2002/sellmeier_20C, absolute
mgf2_o = default_model("MgF2", "o")  # a uniaxial crystal needs its axis
n = caf2.refractive_index(1030e-9)
beta2 = mgf2_o.gvd(800e-9)  # s^2/m
default_reason("CaF2")  # why this model
models("CaF2")  # every model of CaF2, by name
```

## The recommended models

<!-- crystals:start -->
| Material | Axis | Model | Range (nm) | Sample | Index relative to | n at 500 nm | n at 1 µm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CaF2 |  | `daimon_2002/sellmeier_20C` | 138-2326 | 20 °C | nitrogen, 1013.25 hPa | 1.43690 | 1.42931 |
| MgF2 | o | `dodge_1984/sellmeier_o` | 203-7042 | 19 °C | room air (inferred) | 1.38016 | 1.37396 |
| MgF2 | e | `dodge_1984/sellmeier_e` | 203-7042 | 19 °C | room air (inferred) | 1.39207 | 1.38557 |
| BaF2 |  | `malitson_1964/sellmeier` | 265-10349 | 25 °C | room air (inferred) | 1.47824 | 1.46895 |
| Al2O3 | o | `weber_2003/sellmeier_o` | 220-5001 | 20 °C | room air (inferred) | 1.77475 | 1.75616 |
| Al2O3 | e | `weber_2003/sellmeier_e` | 220-5001 | 20 °C | room air (inferred) | 1.76655 | 1.74828 |
| N-BK7 |  | `schott_2023/sellmeier` | 313-2326 | 20 °C, catalogue | air, 1013.3 hPa | 1.52184 | 1.50791 |
<!-- crystals:end -->

**Relative to a gas, made absolute.** A prism refractometer measures the index relative
to the gas around the prism.
- Daimon & Masumura measured CaF2 in nitrogen at 1013.25 hPa and state it.
- SCHOTT's catalogue indices are relative to air at 20 °C and 1013.3 hPa, as its
  Technical Information TIE-29 states.
- The NBS papers (Malitson's CaF2, BaF2 and sapphire; Dodge's MgF2; Malitson &
  Dodge's sapphire of 1972) measured in room air and do not say so. The reading is supported by CaF2: made absolute, Malitson's
  formula agrees with Daimon's within 3.3e-5 at one temperature, and read as absolute
  it would be 3.9e-4 low.

Each model multiplies by the index of its gas:
- nitrogen: optmat's `fit/N2` at the ideal-gas density;
- air: Edlén's dry air.

Treating coefficients measured relative to a gas as absolute indices omits this
reference-medium correction (TC-545, TC-546).

**Sapphire and N-BK7 come from reprints.**
- **Sapphire.** Malitson & Dodge's 1972 formulas of both rays appeared only as a talk
  abstract and in handbooks. optmat reads them in Weber's *Handbook of Optical
  Materials* (2003).
  - They are recommended as a pair, so that the birefringence (0.0091 at 0.3 µm to
    0.0056 at 5 µm) comes from one determination, and they reach 220 nm.
  - Malitson's 1962 formula, fitted to 46 measured ordinary indices at 24 °C, is the
    check (`malitson_1962/sellmeier_o`). It lies 5-8e-5 above the ordinary ray at
    0.4-2 µm, about what the 4 K difference in temperature gives.
- **N-BK7.** SCHOTT's data sheet gives the catalogue formula of a melt of nominal nd
  and vd. A given melt departs from it by up to its tolerance grade (±5e-4 in the
  widest). The data sheet's internal transmittance and dn/dT are curated
  (`schott_2023`), but the model is lossless: its index has no imaginary part, and a
  propagation code that needs the glass's absorption takes it from the curated
  transmittance. `BK7` (the older catalogue glass) and `N-BK7` are distinct materials.

**How well each formula holds its paper.** With its reference gas removed, each formula
reproduces its paper's table (`tests/solid/test_crystals.py`):
- Daimon's reproduces the 69 measured indices within 1.9e-6 (rms 6e-7).
- Dodge's and Malitson's reproduce their computed columns to the printed digit, apart
  from a few printed values their authors do not explain (TC-546).

The papers state their accuracy as follows:
- Daimon: 6e-7 (rms of the fit), with a stated uncertainty per wavelength;
- Dodge: "better than 2e-5";
- Malitson: 2e-5 to 3e-5.

**Temperature.** Each model is at its paper's temperature, between 19 and 25 °C. For
CaF2 near room temperature, dn/dT is about −1.1e-5 per K: 4 K moves n by 4e-5, more
than the formulas' uncertainty. `daimon_2002/sellmeier_25C` gives CaF2 at 25 °C, and
Daimon's measured dn/dT is `daimon_2002/dndt`.

**Not yet recommended:** CaF2 beyond 2.33 µm. Use `malitson_1963/sellmeier` (to
9.7 µm); Daimon's infrared term is weakly determined.

# The metals

Silver, gold, copper and aluminium: the mirror coatings and hollow-waveguide claddings
of ultrafast optics. A metal's optical constants are measured piecewise, each range by a
different technique. So each recommended model is a **table model**: curated tables
joined by range, one source per range, with the reason for each (TC-544, TC-619,
TC-646, TC-647; {doc}`solid`, *Tables of optical constants*).

```python
from optmat.solid.models.registry import default_model, default_reason

silver = default_model("Ag")  # table/Ag
n = silver.refractive_index(800e-9)  # n + ik, k > 0
eps = silver.permittivity(800e-9)
default_reason("Al")  # what the aluminium table rests on
```

Read {doc}`solid` for the complex-index convention and interpolation interface;
the {doc}`quickstart` gives a first calculation. Source identifiers below can be
looked up in the {doc}`the-data/sources` catalogue.

## The recommended tables

<!-- metals:start -->
| Metal | Range (nm) | 100 nm | 400 nm | 800 nm | 1030 nm | 1550 nm | 10 um |
|---|---|---|---|---|---|---|---|
| Ag | 4.1-12398 | 1.274 + 0.594i | 0.053 + 2.140i | 0.054 + 5.606i | 0.073 + 7.387i | 0.139 + 11.311i | 4.155 + 71.291i |
| Au | 3.5-12398 | 1.460 + 0.745i | 1.569 + 1.838i | 0.079 + 4.836i | 0.102 + 6.609i | 0.190 + 10.353i | 5.863 + 68.134i |
| Cu | 3.1-55556 | 1.112 + 0.731i | 1.239 + 2.078i | 0.061 + 5.135i | 0.084 + 7.003i | 0.186 + 10.949i | 7.296 + 73.222i |
| Al | 0.8-200000 | 0.034 + 0.718i | 0.490 + 4.862i | 2.795 + 8.449i | 1.266 + 9.967i | 1.444 + 15.957i | 25.280 + 89.828i |
<!-- metals:end -->

Each join lies where the two neighbouring sources agree to a few per cent. The one
exception is copper at 12.4 µm, where there is a step of 11 % in k (below).

### Silver, gold and copper

| Range | Silver | Gold | Copper |
|---|---|---|---|
| 70.5 nm-12.4 µm | Babar & Weaver | Babar & Weaver | Babar & Weaver (from 41 nm) |
| 25-70 nm | Werner et al., REELS | Werner et al., REELS | |
| X-ray end to 25 nm | Hagemann et al., Table 4 | Hagemann et al., Table 6 | Hagemann et al., Table 3 (to 41 nm) |
| Beyond 12.4 µm | | | Ordal et al. (1985), to 55.6 µm |

- **Babar & Weaver (2015)** analysed the absorptivity of polished single crystals, 0.1-6
  eV, by Kramers-Kronig. Their paper prints 0.1-6 eV. Weaver's website publishes the
  same analysis to 30 eV (the paper's Table 1 is those values rounded; TC-634).
  - Above 6 eV it rests on other authors' reflectance: Leveque, Olson & Lynch's for
    silver, Olson & Lynch's for gold, and Johnson & Christy's then Hagemann et al.'s for
    copper.
  - Gold's reflectance below 4 eV is Weaver's unpublished absorptivity. The 1981
    compilation lists it at 4.2 K, which the table itself does not state (TC-635).
- **Werner, Glantschnig & Ambrosch-Draxl (2009)** derived the permittivity from
  reflection electron energy-loss spectra (REELS). They give about 10 % accuracy above
  5 eV.
- **Hagemann, Gudat & Kunz's DESY report SR-74/7 (1974)** is the numerical tables
  behind their JOSA letter. Its core is their own transmission of evaporated films,
  13-150 eV, analysed by Kramers-Kronig.
  - Above 65-91 eV it prints 1 - n instead of n (`refractive_index_decrement`). There its
    three decimals limit the tables: they end where the printed 1 - n is still 0.004 or
    more (300 eV for silver, 350 eV for gold, 400 eV for copper).
  - For gold the report gives two versions. Table 6, joined to Canfield et al.'s measured
    reflectance, agrees with REELS. Table 5, fitted to the transmission alone, is 30-40 %
    higher in k.

**The vacuum ultraviolet is not well known.** From 40 to 200 nm the curated sources for
silver and gold differ by 20-50 % in k. Silver at 100 nm illustrates the spread:

| Silver at 100 nm | n + ik | Method |
|---|---|---|
| Babar & Weaver (the table) | 1.274 + 0.594i | Kramers-Kronig of reflectance |
| Werner et al. | 1.288 + 0.919i | REELS |
| Hagemann et al. | 1.627 + 0.608i | Kramers-Kronig of film transmission |

At 128-206 nm Stahrenberg et al.'s ellipsometry of Ag(110) is 14-27 % above the table in
k, and at 6 eV the silver sources give:

| Silver at 6.0 eV (207 nm) | n + ik |
|---|---|
| Babar & Weaver | 1.079 + 1.247i |
| Johnson & Christy | 1.18 + 1.31i |
| Stahrenberg et al. | 1.053 + 1.422i |
| Werner et al. | 0.949 + 1.548i |

Gold is similar. The joins were placed where independent methods agree:
- From 41 to 70 nm Babar & Weaver's values and REELS agree to 3-16 % in n and k (29 %
  in silver's k near 60 nm), and to 4 % at 70 nm, the join.
- The DESY films' k is 25-100 % above both at 40-50 nm: two methods against one. Film
  transmission can include an oxide, as Shiles et al. suspected for aluminium.
- Babar & Weaver's own k is likely least reliable at the top of their range, where their
  analysis extrapolates the reflectance beyond 30 eV. So the table leaves them at 70 nm.

**Copper** is simpler. Babar & Weaver's copper above 6 eV is built on Hagemann et al.'s,
and the two meet at 30 eV (41 nm) to 2 % in n and 5 % in k.

**Beyond 12.4 µm** only copper has a measurement: Ordal et al.'s (1985) infrared
reflectance, 180-800 cm⁻¹. It is 9-10 % below Babar & Weaver in k where both tabulate,
hence the step at the join.
- Silver and gold have only Drude curves there: Bennett & Bennett's, which Ordal et al.
  (1983) reprint, and the Drude fits of Ordal et al. (1983, 1985).
- The tables hold their 12.4 µm values beyond it, with a warning. A caller that needs a
  refusal instead takes `model.with_out_of_range("raise")`: outside its knots the table
  then raises `StateError` (OEP 1, section 6.10). The status is kept, since the numbers
  it returns are the same.

### Aluminium

| Range | Source |
|---|---|
| 0.79-31 nm | Rakic (1995), Table 2 |
| 31.8 nm-31 µm | Shiles, Sasaki, Inokuti & Smith, as tabulated by Weaver et al. (1981) |
| 33-200 µm | Ordal et al. (1988) |

**Shiles et al.'s analysis is the core.** It is a self-consistent Kramers-Kronig
analysis of the measured reflectance, ellipsometry, transmission and energy-loss data of
bulk aluminium, 0.04-1551 eV. It is tested against the sum rules, the dc conductivity
and the stopping power.
- The paper itself (Phys. Rev. B 22, 1612 (1980)) says the table "will be given
  elsewhere". It was printed in Weaver et al.'s Physics Data 18-2 (pp. 79-86).
- To make the sum rules hold, Shiles et al. reduced the measured absorption above the L
  edge by 14 %, suspecting the films' oxide.
- Ordal et al.'s (1988) table rests on it above 320 cm⁻¹ and lies 1-3 % above it.

**Rakic (1995) continues it into the EUV.** Above the plasma edge Rakic's k is Shiles
et al.'s, from Palik's handbook. Its n is Rakic's own Kramers-Kronig result, within
0.3 % of Shiles's. Its five figures resolve what the three-decimal tabulation cannot:
- above 40 eV the tabulation's k is 0.00x, and it is 0.000 from 500 eV;
- Rakic samples the L edge (72.7 eV) every 0.05-0.1 eV.

## Films against bulk

A mirror's coating is an evaporated film, and films differ from bulk metal. The curated
film data are kept as datasets, not as tables:

| Source | Metals | Range | Against the tables |
|---|---|---|---|
| Johnson & Christy (1972), with stated errors | Cu, Ag, Au | 188-1937 nm | k within -7/+4 % (Cu), -5/+23 % (Ag), -13/-1 % (Au); n 2-4 times higher in the infrared |
| McPeak et al. (2015), template-stripped | Ag, Au, Cu | 300-1700 nm | Ag and Cu ε₁ within 3 %; Au \|ε₁\| 15-18 % larger; ε₂ higher (Au by 20-70 %) |
| McPeak et al. (2015), template-stripped | Al | 200-1700 nm | \|ε\| 21-31 % smaller than Shiles et al.'s |
| Rakic (1995), intrinsic constants from films | Al | 0.12 nm-200 µm | within 1-2 % in k at 5-10 eV; 2-4 % in k and 9 % in n below 2 eV |

Babar & Weaver attribute the difference to the films' morphology: polycrystalline,
rough and "far from bulk-like". For copper they find the films' infrared ε₂ 2-4 times
the crystal's.

<span id="against-solprop"></span>

## Source and interpolation differences

The choice of source and interpolant can change the optical constants appreciably:

- **Silver, 70–207 nm:** Werner et al.'s REELS data to 124 nm and Stahrenberg et al.'s
  ellipsometry to 207 nm differ from Babar & Weaver's analysis, used by the recommended
  table. The latter is up to 18 % lower in $k$.
- **Aluminium, 150 nm–1.7 µm:** McPeak et al.'s film is 12–17 % below the recommended
  table in $k$.
- **Aluminium, 1.8–31 µm:** Ordal et al.'s data are 1–2 % above the table.
- **Interpolation:** Babar & Weaver tabulate in 0.05 eV steps. Near 10 µm the
  wavelengths are 4 µm apart; a cubic spline and the model's monotone interpolant
  differ by 17 % in $n$ and 0.4 % in $k$. The latter dominates reflectance. The data
  do not establish which interpolant is closer to the material's response; this
  spread indicates sensitivity to interpolation.

## Limits

- **Uncertainty.** Few sources state it per point:
  - Johnson & Christy give Δn and Δk;
  - Weaver et al.'s compilation claims ±10 % "in most cases";
  - Werner et al. claim about 10 % above 5 eV.

  In the vacuum ultraviolet the spread between sources (20-50 %) is the better measure.
- **Stahrenberg et al.'s numbers** reached the database privately. Whether they are the
  roughness-corrected dielectric function is not established (TC-637).
- **Henke et al.'s atomic scattering factors** are kept for the four metals
  (`data/external/henke/`), but no table uses them.
  - For aluminium with Rakic's atom density, they follow Rakic's k to about 2 % at
    1 keV.
  - Near the plasma edge and the L edge they are 20-50 % off, as an independent-atom
    model must be.
  - The DESY report and Weaver et al. print no densities for the noble metals.

  A table can take a Henke segment (`henke = "<element>"` with an `atom_density`
  dataset) where that is wanted.
- **Group-velocity dispersion** of a metal from a table is the interpolant's (see
  {doc}`solid`). For a reflection the phase on reflection matters more, and it follows
  from n + ik.

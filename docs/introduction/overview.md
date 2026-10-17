# Functionality overview

**Curated, validated optical properties of materials, with the provenance of every
number.**

optmat provides the optical response of materials used in ultrafast and nonlinear
optics, in gases, glasses, crystals and metals:
- the refractive index and its dispersion, from the extreme ultraviolet to the
  infrared;
- the complex index across absorption thresholds and bands;
- the Kerr and third-harmonic coefficients;
- the Raman response of the molecular gases, dry air and fused silica;
- the strong-field ionisation rates of the noble gases and of molecules, and the
  parameters those rates are built from;
- two-photon continuum cross sections and bound-excitation records for noble gases.

Each value can be traced to the paper, table and page it came from, together with:
- the technique and the conditions of the measurement;
- the state the value refers to;
- any correction applied to what was printed.

The package provides published formulas, curated tables, and fitted models. Their
descriptions state which sources they rest on and what is known about their accuracy.

The source code, issue tracker and contribution history are at
[LupoLab/optmat on GitHub](https://github.com/LupoLab/optmat).

## Why optmat

Reliable simulations need material data whose conditions, conventions and limitations
are understood. Our motivation is to fix transcription errors, incorrect readings of
measurement conditions and inappropriate applications of models, and to evaluate data
that have not been systematically assessed. optmat provides reference implementations
of published fits and develops new fits from the best available data and analysis.
Their provenance and the calculations behind them remain available for inspection and
reproduction.

[refractiveindex.info](https://refractiveindex.info) is an invaluable resource for the
optics community, and we do not seek to replace it. Our goals are different and
complementary: systematic evaluation of source data, explicit physical conventions,
reproducible fits and validated implementations for simulations.

## From sources to predictions

**Every number keeps its provenance.** Values are taken exactly as printed, by script
from the paper's tables or by hand from its text with the page. The source record says
how each was measured or calculated, at what state, with what uncertainty, and whether
it is someone else's value reproduced. A misprint is corrected by an *erratum* that
records the printed value, the correction, the reason and the evidence; the printed value
itself is never edited ({doc}`/data`).

**Every check is on record.** The review register ({doc}`/register`) records what was
transcribed and how it was checked, every correction, every decision taken in fitting,
and every question left open. The code, the data and this manual cite it by entry.

**The models are physical.** A gas's dispersion is modelled as the dipole
oscillator-strength distribution of its atoms or molecules:
- the resonance lines at their spectroscopic energies;
- effective oscillators for the rest of the distribution;
- the vibrational bands of the infrared-active molecules from HITRAN;
- ab initio values and sum rules where they exist.

The spectroscopic lines constrain the model alongside the fitted data. Read
{doc}`/fitting` for the construction and {doc}`/accuracy` for the supported accuracy;
the physical form alone does not establish accuracy outside the qualified range.

For fused silica, Franta et al.'s universal dispersion model represents the complex
dielectric function through its excitations:
- interband and excitonic transitions;
- the Urbach tail;
- phonon bands.

Its construction preserves causality, passivity, and the sum rule. Other solid
models implement published formulas or interpolate curated tables; {doc}`/solid`
explains these representations and their derivatives.

**Density is handled carefully.** A measurement is reduced with the density its source
meant, real-gas or ideal. The conversion to a refractive index goes through
Lorentz–Lorenz with the refractivity virial coefficient and the magnetic term
({doc}`/conventions`). At the level of parts per million these corrections bring
precision metrology and ab initio theory into agreement.

**Uncertainties are stated and tested:**
- every fit stores its covariance;
- leave-one-dataset-out refits show the sensitivity of fitted models to their data;
- measurements not used in a fit are kept as validators;
- the accuracy pages of this manual are checked against the stored fits by the test
  suite.

**The nonlinear response is process-aware.** The second hyperpolarizability is modelled
as a function of all four frequencies of a process, so one model gives the Kerr effect,
third-harmonic generation and four-wave-mixing channels. The fitted data and
resonance limits determine its support ({doc}`/nonlinear`).

**Responses are checked against experiment and exact theory.** The noble gases'
ionisation rates are validated against numerical solutions of the Schrödinger equation
and against measured yields ({doc}`/ionisation_rates`). The Raman lines are checked
against measured gains ({doc}`/raman`).

## Available data and models

- **The curated records:** 762 source records: 8481 datasets and 197,804 values, 259
  errata and 780 published formulas.
- **15 gases:** He, Ne, Ar, Kr, Xe, H2, D2, N2, O2, CO, CO2, CH4, SF6, N2O and O3.
  - A fitted refractive-index model for every gas:
    - helium and argon at the accuracy of ab initio theory;
    - ozone as a complex index built from its measured absorption;
    - how well each model is known, and how it compares with the formulas in use,
      on {doc}`/accuracy`.
  - The complex index of the noble gases in the extreme ultraviolet ({doc}`/euv`).
  - A model of the second hyperpolarizability of every gas, with its accuracy on
    {doc}`/nonlinear_accuracy`.
  - The Raman response:
    - the Raman lines of N2, O2, H2, D2, CO2, CH4 and SF6, and N2O's vibrational
      bands ({doc}`/raman`);
    - the response functions of the gases and of dry air ({doc}`/raman_response`).
  - Ionisation:
    - the parameters of every rate in use ({doc}`/ionisation`);
    - validated rates of the noble gases ({doc}`/ionisation_rates`);
    - a rate for the molecules, with the forms of the literature beside it
      ({doc}`/molecular_ionisation_rates`);
    - two-photon ionisation cross sections and curated bound-level information for
      the noble gases ({doc}`/two_photon`).
- **Solids:** 95 of the records hold data of glasses, crystals and metals: 2027
  datasets for 163 materials:
  - the linear index and the band gap;
  - n2 and the Raman response;
  - two-photon and multiphoton absorption;
  - ionisation and laser damage.
- **Fused silica** ({doc}`/silica`):
  - its complex index from 25 nm to 125 µm, Franta et al.'s universal dispersion model
    refitted to the absolute refractometry, which it follows to about 2e-5 in the
    transparent window;
  - its n2 from 350 nm to 3.5 µm, with the Raman fraction set by the measured Raman
    gain and the dependence on pulse duration;
  - its Raman response function and consistent nonlinear bundles ({doc}`/silica`).
- **Crystalline silicon** ({doc}`/silicon`): its complex index from 31 nm to 143 µm at
  20 °C, Franta et al.'s universal dispersion model refitted to the absolute
  refractometry and the measured absorption edge.
- **The window crystals and N-BK7** ({doc}`/crystals`): CaF2, MgF2, BaF2, sapphire and
  SCHOTT's N-BK7, each its source's Sellmeier formula made absolute.
- **The metals** ({doc}`/metals`): silver, gold and copper from 3-4 nm to 12.4 µm
  (copper to 56 µm), aluminium from 0.8 nm to 200 µm, as curated tables joined by
  range.

Where the package goes next is on the {doc}`/roadmap`:
- the models of the other solids, and their two-photon absorption, ionisation and
  damage;
- the molecules in the extreme ultraviolet;
- line mixing in the dense molecular gases, and the remaining reference rates.

## Finding a starting point

Follow {doc}`installation` and the {doc}`/quickstart` to evaluate the properties
above. The {doc}`/user-guide/index` develops their conventions, physical models,
and validation in reading order. {doc}`curation` introduces the source records,
datasets, and review register; {doc}`/the-data/index` gives the full reference.

The {doc}`/api/index` is organised into core, gas, and solid interfaces. Package
architecture and contribution workflows belong to the {doc}`/developers/index`.

## AI-assisted development and verification

AI agents assist with literature curation, code, tests, and documentation under the
direction of John C. Travers. Source-specific extraction scripts and checks against
published values support the curated records. The review register records
corrections, decisions, and unresolved questions.

The tests check physical identities, numerical references, and executable examples.
Separate verification runs reproduce extractions and stored fits. {doc}`/about`
explains the roles of AI assistance, source verification, and repeatable
calculations.

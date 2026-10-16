# OEP 7 — Strong-field dissociation of small molecules: the products a pulse leaves

| | |
|---|---|
| **OEP** | 7 |
| **Title** | Strong-field dissociation of small molecules: the products a pulse leaves |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 (2026-10-09) |
| **Status** | Draft |
| **Review status** | Not reviewed |
| **Reviewers** | — |
| **Type** | Feature |
| **Created** | 2026-10-09 |
| **Resolution** | — |
| **Requires** | [OEP 1](archive/oep-0001-propagation-interface.md): model identity, digests and the declared public surface; [OEP 6](archive/oep-0006-two-photon-ionisation.md): the rule that array evaluation returns NaN with an availability mask at unavailable points |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | — |
| **Manual pages** | — |

## Part I — Proposal

### 1. Abstract

An intense near-infrared pulse leaves more than ions and electrons in a molecular gas. It
leaves molecular ions in excited states, some of which fall apart; fragment ions;
dications that explode; and, in fragile molecules such as ozone, neutral fragments from
absorption alone. optmat gives the strong-field ionisation rate of a molecule into its
ion's ground state only. This OEP proposes a model of the **products of one molecule at
the end of a pulse** for O2, N2, N2O, NO2, O3 and N2O5, at 800 nm from a few cycles to
30 fs and up to about 10¹⁵ W cm⁻². It has four parts:
- ionisation into each state of the ion;
- each state's fragmentation;
- the dications;
- neutral channels.

It also gives one labelled **recombination-limited bound** on the atoms the products can
become. The plasma kinetics after the pulse stays outside optmat; optmat supplies its
curated rate data. The model succeeds if every product yield is traceable to curated
records with a stated status, if it reproduces the published model tables it is built
from, and if it says plainly where it fails. It is known to fail the measured N2⁺ state
populations without post-ionisation coupling (section 4), and it has no data for NO2,
O3 and N2O5 ionisation.

### 2. Motivation and scope

**The need.** A model of the atomic oxygen and nitrogen a femtosecond pulse produces
needs, for each molecule, what the pulse leaves behind:
- how many ions;
- in which states;
- how many of them dissociate;
- what neutral fragments appear without ionisation.

This applies to air, to nitrogen-oxide traces, and to ozone and N2O5 in atmospheric or
plasma-chemistry settings. The maintainer's priority is the total atomic yield after
the pulse, with all six molecules equally wanted.

The literature curated for this OEP (TC-954 to TC-981) shows the yield has two parts:
- **Within the pulse.** Each molecule's own products, including unimolecular breakup
  within about a picosecond.
- **After the pulse.** Collisions and recombination over picoseconds to microseconds,
  which depend on pressure, electron temperature and composition.
  - For N2 the second part dominates. Kleine et al. see only a weak prompt N⁺ signal,
    and N + N⁺ appearing collisionally with a 19.8 ps time constant at 4 bar backing
    pressure (TC-958).
  - Peters et al.'s kinetics model of femtosecond-ionised N2 makes its N atoms by
    dissociative recombination of N2⁺ and its cluster ions (TC-977).

The first part is a property of the molecule and the pulse, which optmat can own. The
second is a kinetics problem for the user's model or a propagation code (decided by the
maintainer, 2026-10-09; TC-954).

**In scope:**
- **The six molecules:** O2, N2, N2O, NO2, O3 and N2O5, each with a status, and
  unavailable channels marked as such.
- **Ionisation into individual ionic states:** the inner orbitals (HOMO−1, HOMO−2) as
  well as the HOMO, with curated vertical energies and structure coefficients.
- **State to products:** each ionic state's fragmentation, from curated photoionisation
  branching ratios and dissociation limits.
- **Post-ionisation laser coupling** among the ion's states, in a later phase, for the
  measured N2⁺ populations.
- **Dications:** their states and Coulomb-explosion channels, with kinetic energy
  release.
- **Neutral channels:** one-photon absorption into dissociative states (ozone's
  Chappuis and Wulf bands), and multiphoton channels where a cross section is curated.
- **The products a pulse leaves per molecule**, with the energy the pulse deposits.
- **The recombination-limited atom bound** (section 4), labelled `provisional`.
- **A data accessor** for the curated dissociative-recombination coefficients.
- **The work's own deliverables:** tests, docstrings, a manual page, a gallery example
  and the register entries.

**Out of scope:**
- **Plasma kinetics after the pulse:** recombination dynamics, clustering, charge
  transfer, attachment, electron-impact dissociation, and transport. optmat states the
  rate data; it does not integrate them.
- **Focal-volume averaging** and propagation (Q5).
- **Molecular alignment and its revivals.** The rates are orientation-averaged, as
  `molecular_rate_model`'s are.
- **Circular polarisation**, which has no validated molecular rate in optmat.
- **Higher charge states** than the dication.
- **A change to `molecular_rate_model`'s recommended total rates** (section 10).

### 3. Context

- **Manual:**
  - *Molecular ionisation rates* ({doc}`../molecular_ionisation_rates`): the PPT/MO-ADK
    rate of a molecule into its ion's ground state, its Coulomb factor, the Tong–Lin
    correction, the electron count and its validation against static all-electron
    rates and measured yields.
  - *Ionisation rates* ({doc}`../ionisation_rates`): the pulse integration
    (`pulse_yield`, `cycle_pulse`).
  - *The curated data* ({doc}`../data`): the vocabulary.
  - *Propagation* ({doc}`../propagation`): the surface a propagation code builds on,
    OEP 1.
- **Modules:**
  - `optmat.gas.molecular_ionisation_rates` (`molecular_rate_model`,
    `MolecularRateModel`, `moadk_rate_au`);
  - `optmat.gas._tunnelling` (the structure factors and their orientation average);
  - `optmat.gas.ionisation_validation` (`pulse_yield`, `cycle_pulse`);
  - `optmat.gas.raman_response` (`raman_response_status`, the precedent for asking a
    gas's status without building its model);
  - `optmat.identity` (`Identified`, `Status`).
- **Data:**
  - **The quantities added for this work** (TC-954):
    - `dissociation_energy`;
    - `ionic_state_population`;
    - `fragmentation_branching_ratio`;
    - `kinetic_energy_release`;
    - `dissociation_time`;
    - `dissociative_recombination_rate_coefficient` (its law in `parameters.law`);
    - `photodissociation_quantum_yield`;
    - `transition_dipole_moment`.
  - **Existing quantities:**
    - `ionization_energy` with a final state;
    - `tunnelling_structure_coefficient`;
    - `appearance_energy`;
    - `absorption_cross_section`.
  - **Materials registered by reference only** (TC-954): NO2, N2O5, NO, NO2⁺, NO⁺,
    N2²⁺ and O2²⁺.
- **Register:**
  - TC-446, TC-675: Zhao et al.'s structure coefficients (2010, 2017) behind the HOMO
    rates. N2O has only Zhao 2017's.
  - TC-508, TC-510: Talebpour et al. (1999) and Guo et al. (1998). Their measured
    molecular yields include dissociative channels; Talebpour gives N⁺/N2⁺ and
    O⁺/O2⁺ ≈ 10⁻² in their experiment.
  - TC-954 to TC-981: the curation for this OEP.
- **Roadmap:** none. This is new.

### 4. Theory

**4.1 Ionisation into individual ionic states.** Each state $s$ of the ion is reached by
removing an electron from one orbital. Its static rate is the MO-ADK/PPT rate of that
orbital, with:
- the state's vertical ionisation energy $I_s$;
- its asymptotic coefficients $C_{lm}^{(s)}$;
- its electron count $N_s$.

These are the same formulas `molecular_rate_model` uses for the HOMO
({doc}`../molecular_ionisation_rates`, *The rate*):

$$
w_s(F) = N_s\, w_{\mathrm{PPT}}\!\left(F;\, I_s,\, \{C^{(s)}_{lm}\},\, m_s\right),
$$

orientation-averaged as $\langle B^2(m')\rangle = \sum_l (C^{(s)}_{l})^2 Q^2/(2l+1)$.

In the **sudden limit** the ion's states are populated independently and not coupled
afterwards. For neutral population $P_0$ and ionic-state populations $P_s$:

$$
\frac{dP_0}{dt} = -\sum_s w_s(F(t))\,P_0, \qquad \frac{dP_s}{dt} = w_s(F(t))\,P_0 .
$$

The HOMO channel with the existing default parameters reproduces `molecular_rate_model`
exactly, and the inner channels add to it (Q1, resolved).

The inner orbitals' rates fall steeply with $I_s$. For N2 the sudden limit puts N2⁺(B)
roughly an order of magnitude below N2⁺(X); Kleine et al. state this tunnelling
expectation. Kleine et al. measured instead X : A : B = 0.47 : 0.07 : 0.46 of the ions
at 800 nm and 4.5×10¹⁴ W cm⁻² (TC-958). **The sudden limit is therefore known to fail
for N2.** That is why section 4.2 exists.

**4.2 Post-ionisation coupling (Phase 2).** Yuen and Lin's density-matrix model (TC-955)
keeps the ion's states coherent after ionisation. Their reduced density matrix $\rho$
evolves under the field-dressed Hamiltonian
$H_{ss'}(t) = I_s\,\delta_{ss'} - \mathbf{d}_{ss'}\cdot\mathbf{F}(t)$, with a source
from the neutral and a loss to the dication:

$$
\dot\rho = -\tfrac{i}{\hbar}[H,\rho] + \mathrm{diag}\!\left(w_s P_0\right) - \Gamma^{(2)}\rho .
$$

Here $\Gamma^{(2)}$ holds the state-resolved rates of the second ionisation. The
transition dipoles $\mathbf{d}_{ss'}$ are curated:
- **N2⁺:** X–A and X–B in Yuen and Lin; against bond length in Jia et al. (TC-957).
- **O2⁺:** a–b in Yuen et al. (TC-956).

Yuen et al. (2023) and Jia et al. (2024) extend it to O2 and to dication-state coupling.
Their published populations are figures, not tables (TC-955 to TC-957), so reproducing
them needs vector digitising under rule 8.

**4.3 From ionic states to products.** At the end of the pulse each state $s$ goes to
products $k$ with a branching $B_{s\to k}$, $\sum_k B_{s\to k} = 1$:

$$
Y_k = \sum_s P_s(t_\mathrm{end})\, B_{s\to k}.
$$

A state wholly below its lowest dissociation limit has $B_{s\to\mathrm{parent}} = 1$.
N2⁺ X, A and B and O2⁺ X are expected to be such states, but their limits are not yet
curated (section 5). A state above it takes its branching from state-resolved
photoionisation data where curated. N2O⁺ is the one held case: Berkowitz and Eland give
the states' fragments, ratios and appearance energies (TC-965).

These are photoionisation branchings, and they are taken to hold after strong-field
ionisation into the same state. The assumption is checked against:
- Kübel et al.'s N–N over N–O bond-breaking ratio of 3.5 (TC-966);
- Zhu et al.'s NO2⁺ : NO⁺ ≈ 3 at 810 nm (TC-967; NO2 has no ionisation model yet).

The O2⁺ states above the O⁺ + O limit (b ⁴Σg⁻ and above) need their predissociation
branching from a source not yet held (Q3). De et al. show that bound O2⁺ a ⁴Πu
wavepackets dissociate when a later field couples them to f ⁴Πg (TC-962). That is a
coupling effect, outside the sudden limit.

**4.4 Dications (Phase 3).** The dication's states are populated from the ion's states
by the second ionisation, at the effective binding energies Yuen et al. define. These are
curated as `ionization_model_parameter`, not ionisation energies (TC-956). The dication
states and their energies above the neutral are curated for N2 (Yuen and Lin; Jia et al.
against bond length) and for O2 (Yuen et al.). So are their dissociation limits and
kinetic energy releases. Metastable dication states count as N2²⁺ or O2²⁺; the rest give
two ions.

**4.5 Neutral channels.** A neutral molecule $M$ absorbing one photon with cross section
$\sigma(\omega)$ and dissociating with quantum yield $\phi(\omega)$ does so with
probability

$$
P_1 = 1 - \exp\!\left[-\phi\,\sigma(\omega)\,\Phi\right], \qquad
\Phi = \int \frac{I(t)}{\hbar\omega}\,dt ,
$$

where $\Phi$ is the photon fluence. For a Gaussian pulse of peak intensity $I_0$ and
intensity FWHM $\tau$, $\Phi = I_0\tau\sqrt{\pi/(4\ln 2)}/\hbar\omega$.

For O3 at 800 nm, optmat's absorption model gives
$\sigma = 1.5\times10^{-22}$ cm², which agrees with El Helou et al. to about 1 %
(TC-972). Ozone's Chappuis and Wulf bands dissociate to O2 + O. The quantum yield of
order one is the subject of Grebenshchikov et al.'s review; its printed yields are
reviews of other work and are not curated (TC-973, Q4). At $I_0 = 10^{14}$ W cm⁻² and
$\tau = 30$ fs, $\Phi \approx 1.3\times10^{19}$ cm⁻² and $P_1 \approx 2\times10^{-3}$.

Two other one-photon cases need care:
- **NO2** absorbs one 800 nm photon, but its NO + O dissociation energy exceeds the
  photon's energy, so one photon cannot dissociate it. The limit is not yet curated
  (section 5) and is quoted here from general knowledge.
- **N2O5:** its curated cross sections end at 398 nm (TC-974); no absorption at 800 nm
  is curated.

Multiphoton dissociation into a channel of order $n$ would have the rate
$W_n = \sigma_n (I/\hbar\omega)^n$. No $\sigma_n$ for these molecules is curated; Song
et al. give only O2's photon order (about 9 to 11) and timescale (TC-960). Neutral
multiphoton channels are therefore **unavailable**, never zero (section 6).

**4.6 The recombination-limited bound.** Each product ion $j$ recombining dissociatively
with an electron gives $a_{j,X}$ atoms of element $X$, by its curated dominant channel
(TC-976):
- O2⁺ → O + O;
- N2⁺ → N + N;
- NO⁺ → N + O;
- NO2⁺ → NO + O.

The bound on atoms of $X$ per molecule is

$$
A_X^{\max} = A_X^{\mathrm{prompt}} + \sum_j Y_j\, a_{j,X} ,
$$

assuming every ion recombines dissociatively and nothing else competes.

The bound fails, as an approximation, where other processes win. For example, in air
at one atmosphere O2⁺ clusters to O4⁺, whose recombination gives O2 + O2. The
manual will state the competing processes, with their coefficients where curated
(Kossyi et al., TC-976, includes the clusters in its notes). The recombination time
$1/(k n_e)$ uses the curated $k(T_e)$, e.g. $2\times10^{-7}(300/T_e)$ cm³ s⁻¹ for O2⁺.

**4.7 Energy deposited.** The pulse deposits $\sum_s P_s I_s$ per molecule into ionic
states (the HOMO channel's adiabatic energy, as `energy_per_event`). It deposits
$\hbar\omega$ per absorbed photon in the neutral channels. The kinetic energy of
electrons is not counted, as now.

**Units and conventions:** SI at the interface; populations and yields per molecule;
energies in J (eV in the tables). Intensity is the cycle-averaged
$I = \tfrac12 c\varepsilon_0 F_0^2$ of a linearly polarised field, as in
{doc}`../conventions`.

### 5. Primary sources and data

**In hand (curated):**

| Need | Sources | Register |
|---|---|---|
| Ionic-state vertical energies, inner-orbital $C_{lm}$ | Yuen and Lin 2022 (N2); Yuen et al. 2023 (O2, three orbitals); Jia et al. 2024 (N2, against $R$); Zhao et al. 2010, 2017 (HOMO) | TC-955–957, TC-446, TC-675 |
| Transition dipoles of ions and dications | Yuen and Lin; Yuen et al.; Jia et al. | TC-955–957 |
| Dication states, limits, KER | Yuen and Lin; Yuen et al.; Jia et al.; Wu et al. (triple ionisation KER) | TC-955–957, TC-959 |
| N2O⁺ state → fragments | Berkowitz and Eland 1977 | TC-965 |
| Validators: measured state populations | Kleine et al. 2022 (N2⁺ X/A/B) | TC-958 |
| Validators: strong-field channel ratios | Kübel et al. 2014 (N2O); Zhu et al. 2006 (NO2); Singhal et al. 1998 (NO2, figures); Talebpour et al. 1999 and Guo et al. 1998 (O⁺/O2⁺, N⁺/N2⁺) | TC-966–968, TC-508, TC-510 |
| Dissociation times | Kleine et al. (collisional N2⁺); Song et al. (O2 superexcited); Zhao et al. 2023 (N2O⁺) | TC-958, TC-960, TC-971 |
| O3 neutral channel | El Helou et al. 2005; optmat's O3 absorption; Grebenshchikov et al. 2007 (limits) | TC-972, TC-973 |
| N2O5 photochemistry | Harwood et al. 1998; Ravishankara et al. 1986 (disagree at 248 nm) | TC-974, TC-975 |
| Recombination coefficients | Kossyi et al. 1992 (and its copies in Peters et al. and Hooshyar and Dumitrache) | TC-976, TC-977, TC-979 |

**Missing, and requested:** these are not curated as source records, and each is to be
identified and checked before use.
- **O2⁺ predissociation branching.** The state-resolved fragmentation of O2⁺ above the
  O⁺ + O limit, from threshold photoelectron–photoion coincidence data (Q3).
- **Structure coefficients for NO2 and O3** (weak-field asymptotic theory or MO-ADK),
  and for N2O's inner orbitals. Zhao 2017 covers N2O's HOMO only.
- **N2O5 ionisation:** whether N2O5⁺ survives, and its fragment appearance energies.
- **NO2⁺ and O3⁺ state-resolved fragmentation** from photoionisation.
- **Absolute fragment yields** against intensity at 30 fs or shorter, for O2 and NO2.
- **Validators of the total atom yield:** ozone and NOx produced by femtosecond
  filaments; absolute atomic O and N densities after femtosecond excitation.
- **The Supplemental Material of Kleine et al.**, for the pulse duration and
  polarisation of their measurement (TC-958).
- **Dissociation limits of N2⁺, O2⁺ and NO2** (thermochemical), to show which states
  are bound.

**Figures to digitise** (rule 8, `use = "check"`): the register entries of TC-955 to
TC-971 list them. The first needed are:
- Yuen and Lin's yields against intensity (Fig. 5);
- Yuen et al.'s Fig. 9;
- Jia et al.'s Fig. 2;
- Wu et al.'s Fig. 4 channel ratios (raster).

**Uncertainty.** Each channel states:
- its source uncertainty where printed;
- the spread between curated coefficient sets. For O2's HOMO, Yuen et al.'s
  0.69/0.06 against Zhao 2010's 0.52/0.03 and Zhao 2017's 0.64/0.04; Jia et al.'s N2
  values also differ (TC-956, TC-957; Q2).

Each model carries an `uncertainty_factor` with its reason, as the molecular rates do.
A product yield's uncertainty is not a sum in quadrature of these. It is stated per
channel, and validators decide the status.

### 6. Proposed interface

A new module, `optmat.gas.dissociation`.

| Name | Meaning | Units |
|---|---|---|
| `fragmentation_model(gas, wavelength, model="default")` | the model of one molecule at one carrier wavelength | m |
| `.states` | the ion's states: label, vertical energy, structure coefficients, electron count, products table, sources | J |
| `.state_rates(field_strength)` | $w_s$ for each state, instantaneous, linear polarisation | s⁻¹ |
| `.pulse_products(peak_intensity, fwhm=None, *, cycles=None)` | a `Products` record: yield per molecule of each species (parent ion by state, fragment ions, neutral fragments, dications, electrons), the energy deposited, an `available` mask per species, and the model's status | 1; J |
| `.neutral_channels` | each neutral channel: order, cross section or its absence, quantum yield, sources | |
| `.recombination_limited_atoms(products)` | $A_X^{\max}$ per element, labelled `provisional`, with its assumptions | 1 |
| `.identity`, `.digest`, `.status`, `.uncertainty_factor` | OEP 1 | |
| `fragmentation_status(gas)` | `(Status, reason)` or None, without building the model (as `raman_response_status`) | |
| `dissociative_recombination_coefficient(ion, electron_temperature)` | the curated $k(T_e)$ and its source; no kinetics | m³ s⁻¹ |

**Warnings and refusals:**
- **Unavailable channels:** a channel without curated data is unavailable. Evaluating it
  returns NaN with `available = False` (OEP 6's rule), never 0. Examples are NO2's
  ionisation in Phase 1, every neutral multiphoton channel, and O2⁺ b-state branching
  until Q3.
- **Unknown molecule:** a molecule outside the six raises `OptmatError`.
- **Range:** where the HOMO model warns (`ModelRangeWarning`: below 350 nm, above its
  `correction_range`), so does this one. The curated validators are all near 800 nm,
  and the manual will say so.
- **The sudden limit's status:** any model without coupling carries status `provisional`,
  with the reason that it fails Kleine's N2⁺ populations.
- **The bound:** the recombination-limited bound always carries status `provisional`.
  It is never offered as a prediction of atoms at a time.

A sketch (fenced as `text` until it runs; values to be generated):

```text
>>> from optmat.gas.dissociation import fragmentation_model
>>> o2 = fragmentation_model("O2", 800e-9)
>>> [s.label for s in o2.states]
['X 2Pi_g', 'a 4Pi_u', 'b 4Sigma_g-']
>>> p = o2.pulse_products(1e18, fwhm=30e-15)   # 1e14 W/cm2, in W/m2
>>> p.yields["O2+"], p.yields["O+"], p.available["O"]
...
>>> o2.recombination_limited_atoms(p)["O"]     # provisional bound
...
```

**Composition:**
- **With `molecular_rate_model`:** the HOMO channel is the existing model and is not
  changed. The inner channels add to it. A user who wants the total ionisation keeps
  using `molecular_rate_model`; the fragmentation model's total is reported beside it.
  Measured yields used to validate the HOMO model already include dissociative
  channels (TC-508, TC-510). The difference is of order 10⁻², within its uncertainty
  factor.
- **With optmat's absorption models:** the O3 neutral channel uses the gas's curated
  absorption cross section. It does not hold a copy.
- **With a propagation code:** the energy deposited and the electrons are per molecule
  per pulse. A time-resolved loss current would come from `state_rates` and the
  neutral channels' rates. The split between channels is the code's choice.

### 7. Implementation strategies

**Where it belongs:** the gas subpackage (`optmat.gas.dissociation`), beside the
molecular rates, reusing:
- `optmat.gas.ionisation_rates` and `molecular_ionisation_rates` for each state's rate,
  as a `MolecularRateModel` built `from_parameters` with the state's energy and
  coefficients;
- `ionisation_validation`'s pulses for the integration.

**Configuration:** a TOML file per model, `data/gas/dissociation.toml`, in the pattern of
`molecular_ionisation_rates.toml`. Each state names its datasets and its reason
(TC-n). Each products table names the dataset of each branching or the reason for a
fixed value, e.g. bound below the limit, citing the `dissociation_energy` or
`appearance_energy` dataset that shows it.

**Option A: sudden limit first (recommended).** Phase 1 is the rate-equation model of
4.1 with the tables of 4.3. It is cheap and directly traceable. Its known failure for N2
is stated in its status and tested against Kleine as a validator that fails by design.

**Option B: density matrix from the start.** It is correct for N2⁺ but needs:
- the transition dipoles against geometry;
- the field-dressed integration;
- digitised figures to check it.

It would delay a usable model for all six molecules for the sake of the two best served
by data. It is recommended as Phase 2, not as the start.

**Option C: empirical fragment fractions.** For example, a constant 25 % of NO2⁺
dissociating, from Zhu et al. It is rejected as a default: those fractions are measured
under one pulse with doubtful intensity calibration. They are kept as validators.

**The coefficient sets:** where Yuen et al. or Jia et al. and Zhao et al. differ, the
model takes one set per orbital and records the other as a comparison model (Q2,
resolved: the HOMO on the default's set, the inner orbitals on the defining paper's).

**Cost:** the sudden limit integrates a handful of rates over a pulse: milliseconds per
intensity. The density matrix is a few-state ODE per pulse: well under a second.

**Vocabulary:** done (TC-954). The species of products are strings from a controlled
list in the module (O, N, NO, NO2, NO3, O2, N2 and their ions). Atoms are not registered
materials (Q9).

### 8. Verification and validation

- **Reproductions:**
  - The HOMO channel equals `molecular_rate_model`'s default to floating-point
    tolerance, for every molecule that has one.
  - Each state's static rate reproduces `moadk_rate_au` with that state's curated
    $I_s$ and $C_{lm}$.
  - Phase 2 reproduces Yuen et al.'s and Jia et al.'s published populations (digitised,
    vector). The tolerance is the digitising uncertainty, combined with the difference
    between their MO-ADK form and optmat's (Tong–Lin and Coulomb factor). That
    difference is stated and switched off for the comparison where possible.
- **Invariants:**
  - $P_0 + \sum_s P_s + \sum_d P_d = 1$.
  - $\sum_k B_{s\to k} = 1$ for every state.
  - Element balance: atoms in equal atoms out, for every products record.
  - The bound is at least the prompt atoms.
- **Analytic limits:**
  - The neutral one-photon channel against $1-\exp(-\phi\sigma\Phi)$ for a Gaussian
    pulse.
  - The weak-field limit of each state's rate.
- **Validators:**
  - **Kleine's N2⁺ X : A : B.** Phase 1 reports the miss; Phase 2 is judged on it, within
    Kleine's stated errors and the unknown pulse duration (Q7).
  - **Talebpour et al.'s fragment-to-parent ratios** of about 10⁻² for O2 and N2.
  - **Kübel's N2O channel ratios**, and Zhu's NO2 ratio once NO2 ionises.
- **Ozone:** the linear probability against El Helou's cross section at its tabulated
  maxima, interpolated with optmat's own model.
- **solprop:** no parity target exists. A paired solprop proposal, if any, consumes
  `Products` and `state_rates`.

### 9. Documentation

- **A new manual page, *Strong-field dissociation*.** It covers:
  - the two parts of the atom yield and the boundary;
  - the sudden limit and why it fails for N2;
  - the products tables, with their sources;
  - the neutral channels;
  - the bound and its failure modes, with the recombination coefficients;
  - a table of each molecule's status and its unavailable channels;
  - open questions.
- **Docstrings** with units and ranges.
- **A gallery example:** O2 and N2 products against intensity at 30 fs, with the
  bound.
- **Links:** the molecular-rates page links to the new page, and the data page's
  quantities are already documented (TC-954).

### 10. Compatibility and impact

- **No recommended value changes.** `molecular_rate_model` and every default are
  untouched.
- **No stored fit is re-run.**
- **New public names only:** `optmat.gas.dissociation` and its config file.
- **No solprop parity fixture changes.**
- **No register entry is superseded.**

### 11. Alternatives considered

- **A kinetics solver in optmat.** Declined by the maintainer (TC-954). Its result
  depends on pressure, mixture and electron temperature, which are the user's
  conditions, not the material's. Its rate data would be mostly secondary compilations.
- **Fitting fragment fractions to measured ratios against intensity.** The measured
  ratios carry estimated or unstated intensity calibrations (TC-967, TC-508), and some
  would be digitised from raster figures, which optmat never fits.
- **Taking published strong-field model curves (e.g. Yuen et al.'s) as data.** Those
  curves are model outputs. They are used as reproductions of the method, not as
  references, as for the ionisation rates (TC-439's tiers).
- **Leaving the triatomics out.** Rejected by the maintainer's priority (all six
  equally). They are included with explicit unavailable channels instead, so that a
  user sees what is missing.

### 12. Phases

1. **The sudden-limit framework.** Deliverables:
   - `optmat.gas.dissociation` with the interface of section 6;
   - `data/gas/dissociation.toml`;
   - the state channels for O2 (X, a, b) and N2 (X, A, B), and N2O (X only, until
     its inner orbitals are curated);
   - the products tables from section 4.3, with unavailable entries for O2⁺'s b state
     until Q3;
   - the O3 neutral one-photon channel;
   - NO2, O3 and N2O5 models whose ionisation is unavailable, with their status
     `untested`;
   - the recombination accessor and bound;
   - tests of section 8 except Phase 2's;
   - the manual page, the gallery example, register entries and a changelog line.
2. **Post-ionisation coupling.** Deliverables:
   - the density-matrix evolution of 4.2 for N2⁺ and O2⁺, with the curated dipoles;
   - the digitised populations of Yuen and Lin, Yuen et al. and Jia et al. as checks.

   It is judged against Kleine's populations. It changes N2's status only if it passes.
3. **Dications.** Deliverables:
   - the second ionisation from the ion's states;
   - the dication states, channels and KER from the curated tables;
   - Wu et al.'s and Kübel et al.'s ratios as validators.
4. **Triatomic ionisation and neutral multiphoton channels.** Each of these follows its
   curation:
   - NO2, O3 and N2O5 ionisation, once their structure coefficients and ionic
     fragmentation are curated;
   - the neutral multiphoton channels, once a cross section is.

### 13. Open questions

- **Q1 (the total; resolved).** The fragmentation model's total ionisation is the HOMO
  model plus the inner channels; the inner channels are not renormalised to keep the
  total equal to `molecular_rate_model`'s. The difference is the inner orbitals' share,
  small in the sudden limit, and it is reported beside the HOMO model's total (the
  maintainer, 2026-10-09; TC-982).
- **Q2 (structure coefficients; resolved).** Yuen et al. cite Zhao et al. 2010 for O2's
  coefficients but print different values, and Jia et al.'s N2 values differ from both
  of Zhao's sets (TC-956, TC-957). Decided (the maintainer, 2026-10-09; TC-982):
  - the HOMO stays on the default model's set, so that its channel is
    `molecular_rate_model` unchanged;
  - the inner orbitals take the set of the paper that defines them: Yuen et al. (2023)
    for O2's HOMO−1 and HOMO−2; Jia et al. (2024) for N2's, at R = 1.10 Å, the
    tabulated distance nearest equilibrium;
  - the other sets are offered as comparison models.
- **Q3 (O2⁺ above the dissociation limit; before Phase 1 completes O2).** The branching
  of O2⁺ b ⁴Σg⁻ and higher states into O⁺ + O needs a curated photoionisation source.
  Until then those entries are unavailable.
- **Q4 (ozone's quantum yield; Phase 1).** The Chappuis/Wulf quantum yield is taken as
  one, from the predissociative character of the bands. Grebenshchikov et al. review
  but do not measure it. Accept 1 with a stated reason, or wait for a curated
  measurement?
- **Q5 (focal-volume averaging; resolved recommendation).** Leave it downstream. A
  helper may be offered in a later phase if a validator needs it.
- **Q6 (a kinetics validator; resolved).** None in optmat (TC-954).
- **Q7 (Kleine's pulse; Phase 2).** The duration and polarisation of the measurement are
  only in the Supplemental Material, which is not held (TC-958). Request it before
  Phase 2 is judged.
- **Q8 (pulse dependence; Phase 1).** Singhal et al. report that shorter pulses favour
  NO2⁺ over NO⁺ (TC-968). The sudden limit's branching does not depend on the pulse.
  Record this as a limitation, or wait for Phase 2's coupling?
- **Q9 (product species; Phase 1).** Products are named by a controlled list of species
  strings. Should atoms and radicals (O, N, NO3) become registered materials instead?
  Recommended: no, until data are curated for them.

### 14. References

Curated as source records unless marked (the record's id in brackets).

1. C. H. Yuen and C. D. Lin, "Density-matrix approach for sequential dissociative double
   ionization of molecules," *Phys. Rev. A* **106**, 023120 (2022).
   doi:10.1103/PhysRevA.106.023120 [yuen_2022]
2. C. H. Yuen, P. Modak, Y. Song, S.-F. Zhao and C. D. Lin, "Modeling the sequential
   dissociative double ionization of O2 by ultrashort intense infrared laser pulses,"
   *Phys. Rev. A* **107**, 013112 (2023). doi:10.1103/PhysRevA.107.013112 [yuen_2023]
3. Y.-W. Jia, C. H. Yuen, W.-Q. Jing, Z.-Y. Zhou, C. D. Lin and S.-F. Zhao, "Improved
   model for the dissociative sequential double ionization of nitrogen molecules by an
   intense few-cycle infrared laser pulse," *Phys. Rev. A* **110**, 023112 (2024).
   doi:10.1103/PhysRevA.110.023112 [jia_2024]
4. C. Kleine et al., "Electronic state population dynamics upon ultrafast strong field
   ionization and fragmentation of molecular nitrogen," *Phys. Rev. Lett.* **129**,
   123002 (2022). doi:10.1103/PhysRevLett.129.123002 [kleine_2022]
5. J. Berkowitz and J. H. D. Eland, "Photoionization of N2O: Mechanisms of
   photoionization and ion dissociation," *J. Chem. Phys.* **67**, 2740 (1977).
   doi:10.1063/1.435188 [berkowitz_1977]
6. M. Kübel et al., "Strong-field control of the dissociative ionization of N2O with
   near-single-cycle pulses," *New J. Phys.* **16**, 065017 (2014).
   doi:10.1088/1367-2630/16/6/065017 [kubel_2014]
7. J. Zhu, W. Guo, Y. Wang and L. Wang, "Dissociation of NO2 in femtosecond intense
   fields," *Chin. Sci. Bull.* **51**, 1185 (2006). doi:10.1007/s11434-006-1185-2
   [zhu_2006]
8. R. P. Singhal et al., "Comment on 'On the ionisation and dissociation of NO2 by short
   intense laser pulses'," *Chem. Phys. Lett.* **292**, 643 (1998).
   doi:10.1016/S0009-2614(98)00718-0 [singhal_1998]
9. J. Wu, H. Zeng, J. Wang and C. Guo, "Dynamics of triple-ionization-induced
   dissociation in diatomic molecules in strong fields," *Phys. Rev. A* **73**,
   051402(R) (2006). doi:10.1103/PhysRevA.73.051402 [wu_2006]
10. D. Song et al., "Neutral dissociation of superexcited oxygen molecules in intense
    laser fields," *J. Phys. Chem. A* **114**, 3087 (2010). doi:10.1021/jp907662z
    [song_2010]
11. S. De et al., "Following dynamic nuclear wave packets in N2, O2, and CO with
    few-cycle infrared pulses," *Phys. Rev. A* **84**, 043410 (2011).
    doi:10.1103/PhysRevA.84.043410 [de_2011]
12. X. Zhao et al., "Probing electron localization during molecular dissociation by
    femtosecond strong-field ion momentum spectroscopy," *Commun. Phys.* **6**, 124
    (2023). doi:10.1038/s42005-023-01248-3 [zhao_2023]
13. A. Talebpour, J. Yang and S. L. Chin, "Semi-empirical model for the rate of tunnel
    ionization of N2 and O2 molecule in an intense Ti:sapphire laser pulse," *Opt.
    Commun.* **163**, 29 (1999). doi:10.1016/S0030-4018(99)00113-3 [talebpour_1999]
14. C. Guo, M. Li, J. P. Nibarger and G. N. Gibson, "Single and double ionization of
    diatomic molecules in strong laser fields," *Phys. Rev. A* **58**, R4271 (1998).
    doi:10.1103/PhysRevA.58.R4271 [guo_1998]
15. S.-F. Zhao, C. Jin, A.-T. Le, T. F. Jiang and C. D. Lin, "Determination of structure
    parameters in strong-field tunneling ionization theory of molecules," *Phys. Rev.
    A* **81**, 033423 (2010). doi:10.1103/PhysRevA.81.033423 [zhao_2010]
16. S.-F. Zhao, J.-K. Li, G.-L. Wang, P.-C. Li and X.-X. Zhou, "Accurate structure
    parameters for tunneling ionization rates of gas-phase linear molecules," *Commun.
    Theor. Phys.* **67**, 289 (2017). doi:10.1088/0253-6102/67/3/289 [zhao_2017]
17. Z. El Helou et al., "Absolute cross
    sections of ozone at atmospheric temperatures for the Wulf and the Chappuis bands,"
    *J. Chem. Phys.* **122**, 244311 (2005). doi:10.1063/1.1937369 [el_helou_2005]
18. S. Yu. Grebenshchikov, Z.-W. Qu, H. Zhu and R. Schinke, "New theoretical
    investigations of the photodissociation of ozone in the Hartley, Huggins, Chappuis,
    and Wulf bands," *Phys. Chem. Chem. Phys.* **9**, 2044 (2007).
    doi:10.1039/b701020f [grebenshchikov_2007]
19. M. H. Harwood, J. B. Burkholder and A. R. Ravishankara, "Photodissociation of BrONO2
    and N2O5: Quantum yields for NO3 production at 248, 308, and 352.5 nm," *J. Phys.
    Chem. A* **102**, 1309 (1998). doi:10.1021/jp9729829 [harwood_1998]
20. A. R. Ravishankara, P. H. Wine, C. A. Smith, P. E. Barbone and A. Torabi, "N2O5
    photolysis: Quantum yields for NO3 and O(3P)," *J. Geophys. Res.* **91**, 5355
    (1986). doi:10.1029/JD091iD05p05355 [ravishankara_1986]
21. I. A. Kossyi, A. Yu. Kostinsky, A. A. Matveyev and V. P. Silakov, "Kinetic scheme of
    the non-equilibrium discharge in nitrogen-oxygen mixtures," *Plasma Sources Sci.
    Technol.* **1**, 207 (1992). doi:10.1088/0963-0252/1/3/011 [kossyi_1992]
22. C. J. Peters, M. N. Shneider and R. B. Miles, "Kinetics model of femtosecond laser
    ionization in nitrogen and comparison to experiment," *J. Appl. Phys.* **125**,
    243301 (2019). doi:10.1063/1.5098306 [peters_2019]

## Part II — Implementation record

None yet.

## Part III — History

- 2026-10-09: Drafted, after the curation of TC-954 to TC-981.
- 2026-10-09: Q1 and Q2 resolved by the maintainer, as recommended (TC-982).

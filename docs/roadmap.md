# Roadmap

optmat began with the gases used in gas-filled fibres and cells, and with the properties
a simulation of nonlinear propagation in them needs first: the refractive index across
the spectrum and the Kerr response. The aim is the same for every property and every
material that follows. Take the data from the primary literature with their provenance,
check them on record, fit physical models with honest uncertainties, and keep the
published alternatives available for comparison.

A major item starts as an {doc}`enhancement proposal <oeps/index>` (OEP). That means a new
medium, a new family of properties or models, a changed default or convention, or a
redesigned interface. Accepted proposals are listed first, then drafts, then the informal
items.

## Accepted enhancement proposals

OEPs that have been accepted and are not yet final, each waiting to be picked up or
being implemented.

None at present: OEPs 1-6 are Final (see the {doc}`OEP index <oeps/index>`).

## Draft enhancement proposals

OEPs under discussion. Their scope may change, and implementation has not been accepted.
The {doc}`OEP index <oeps/index>` lists every proposal and its current status.

- [OEP 7 — Strong-field dissociation of small molecules: the products a pulse
  leaves](oeps/oep-0007-strong-field-dissociation.md) (Draft). It covers the products one molecule is left with at the end
  of an 800 nm pulse, from a few cycles to 30 fs:
  - ionisation into each state of the ion, then each state's fragmentation (the sudden
    limit first);
  - post-ionisation coupling;
  - dications;
  - neutral channels.

  The molecules are O2, N2, N2O, NO2, O3 and N2O5. It also gives a labelled
  recombination-limited bound on the atoms. The plasma kinetics after the pulse stays
  outside optmat. The data are curated (TC-954 to TC-981); what is missing is listed in
  its section 5.
- [OEP 8 — The bound two-photon resonances of the noble gases for
  propagation](oeps/oep-0008-bound-two-photon-response.md) (Draft; Reviewed).
  Q8–Q10 are resolved by the maintainer (TC-1055). Phase 1 supplies He, Xe and Ar
  coefficients and weak-field applications, plus the local He 1s2s nonlinear
  reference model with its light shifts. It includes the approved decay quantities,
  Drake's on-shell perturbative kernel and Bondy's helium decay rate, after extraction.
  Xe and Ar strong-field multilevel dynamics and general frequency-dependent or
  two-colour dynamics remain deferred. The empirical bandwidth estimate is not a
  universal error bound. The xenon dispersive validation failure and missing shifts,
  photoionisation and pressure data remain explicit; Ne and Kr remain refused until
  their strengths are resolved.

  The initial application is the disappearance of an Ar VUV RDW near the two-photon
  ionisation threshold under a self-compressing 6 fs, 800 nm pump, with He as a second
  target (TC-1056). Start the solprop investigation with OEP 6's continuum coefficients:
  compare pump-ionisation baselines and evaluate perturbative VUV pair-spectral loss
  on the generated fields. OEP 8's lowest bound lines are at different wavelengths.
  This diagnostic need not wait for OEP 8, and it does not provide the continuum's
  missing time-dependent loss current or decide which mechanism dominates.
  optmat supplies material properties and basic data; solprop owns the resulting
  polarisation, electron accounting, plasma coupling and propagation energy balance
  (TC-1057). The internal atomic fixtures test material coefficients and conventions.

The informal items below record work not fully covered by these proposals. When a
proposal is added or revised, an informal item is removed only if the proposal covers it
in full. Otherwise the item links the proposal, states what it covers and keeps the
remaining work explicit. Moving work into a proposal does not mean it is implemented.

## Gases: what comes next

- **VUV RDW suppression in argon and helium** (TC-1056). A solprop investigation,
  using OEP 6's existing continuum coefficients and the bound coefficients as OEP 8
  provides them. Use the evolving complex VUV field, spatial overlap and neutral
  density, not the input pump duration. Resolve the impact of near-threshold gaps and
  unresolved structure; do not replace them by zero. Compare pump ionisation,
  VUV two-photon ionisation, and both, including absorption and plasma feedback.
  A causal continuum loss/ionisation rule belongs in a solprop SEP (OEP 8 Q6;
  TC-1057), with a coupled optmat OEP only for additional material quantities, models
  or interfaces it requires. An effective continuum matched to OEP 6's cross sections
  is a candidate for the SEP, not an approved response; integrated yield agreement
  alone does not validate its transient polarisation or electron accounting.
  Distinguish pump-driven detuning from the ultimate observable cut-off, and the
  spectral centroid from the blue edge and detection limit. The pressure scan alone
  does not distinguish the hypotheses. Fix the dataset's actual input duration before
  specifying simulation fixtures.
- **Beyond OEP 8's approved first scope** (Q8–Q9; TC-1055). Derive and qualify a
  general frequency-dependent/two-colour response in solprop, with a coupled optmat
  proposal for any additional material inputs (TC-1057), and establish when
  excited-state Raman couplings, shared-continuum
  loss and correlated dephasing matter for strong multilevel excitation. Include
  pump-induced dressing where the RDW overlaps the intense pump. These are deferred
  scientific tasks, not consequences of approving Q8–Q10. Higher bound series near
  threshold still need strengths; the existing low levels do not bound their loss.

- **Raman and ionisation parameters for propagation (TC-439).** The aim is that
  solprop holds no material parameters of its own. optmat takes over the conversion from
  the literature, and solprop keeps only the response kernels and the molecules' rate
  formulas.
  - **Done:** the Raman lines of seven gases and N2O's bands ({doc}`raman`), the
    response functions of the gases, dry air and fused silica ({doc}`raman_response`,
    TC-685 to TC-696), the ionisation parameters ({doc}`ionisation`), and the rates of
    the noble gases and the molecules ({doc}`ionisation_rates`, TC-652;
    {doc}`molecular_ionisation_rates`, TC-673 to TC-684).
  - **Next:**
    - the ionisation rates: the disagreement of the two measured yields (TC-660), for
      which a one-electron calculation of the 42 fs pulse is wanted; hydrogen, 400 nm
      and 1030 nm; all-electron static rates of Kr and Xe, and of Ne and Ar below
      1.7e14 W/cm², where the tabulated ones are in doubt (TC-669); for the
      molecules, a calculation at a finite frequency, all-electron static rates of
      O2, CO2, N2O and CH4, and SF6's coefficients (TC-680, TC-681);
    - N2O's rotational lines, once its rotational constants are curated;
    - modern H2 and D2 level energies (Komasa et al.) in place of Stoicheff's;
    - the 11-39 % excess of the model's H2 and D2 gains over the measured ones
      (TC-542);
    - the Q-branch widths of CO2 and N2O, and N2O's depolarisation ratios (stopgaps
      now);
    - **Q-branch line mixing** (TC-527, TC-688). Above about 2 amagat for N2 (and dry
      air), 4 for O2 and 13 for CO2, the Q(J) lines overlap. Collisions then transfer
      coherence between them, and the branch narrows towards a single line. Today the
      response warns there and treats the lines as independent, which overstates the
      branch's dephasing. The plan:
      - a relaxation matrix W for each band, whose diagonal is today's per-J widths
        (so nothing changes at low density), with off-diagonal rates from a fitting
        law such as the modified exponential gap (Koszykowski et al. 1987, curated as
        `ExponentialGapLaw`) or the energy-corrected sudden law, with detailed
        balance imposed;
      - the response from the eigen-decomposition of Omega + iW. This is still a sum of
        damped oscillators, now with complex amplitudes, so the response's form and
        the propagation code's use of it are unchanged. Its eigenvalues depend on the
        density, which the low-rank tabulation already absorbs;
      - validation against the collapsed branch: Mikhailov's widths of N2's Q branch at
        15-114 atm (curated), Rosasco et al.'s first-order mixing coefficients at 0.2-2
        atm (curated, but they fail the detailed-balance sum rule, TC-527), and Millot
        et al.'s O2 spectra;
      - the same for the rotational S lines of CO2, which overlap near 13 amagat;
    - the structure coefficients of CH4, SF6 and N2O;
    - the polarizabilities of the remaining cations;
    - reference rates for the gases without them: the papers requested below, and TDSE
      calculations where the literature leaves a gap (TC-444).
- **The molecules in the extreme ultraviolet.** Brion & Thomson's compilation of
  absolute photoabsorption (H2, N2, O2, CO, CO2, N2O, CH4 and SF6, 5-80 eV, 5 %) and
  Olney et al.'s dipole sums are curated. A complex index of the molecules above their
  ionisation thresholds, built as the noble gases' is, comes next.
- **The infrared beyond OEP 3.** What OEP 3 leaves out:
  - **SF6's infrared treatment.** Its input is a measured, broadened cross-section, so
    it needs a Kramers-Kronig reference with its resolution, baseline and truncation
    errors. Its cross-sections also show weak absorption at 1.5-5 µm, inside its
    declared range (OEP 3, §13 Q14).
  - **A complex index for the weak bands** (OEP 3, Option F). Damped oscillators or
    smoothed band shapes would qualify pure CO2 at its laser lines, and the weak-band
    regions of N2O and CH4.
  - **Air across its trace gases' bands** (OEP 3, Q13; TC-809). The full CIPM-2007 air
    widens only to 1499-1516 and 1548-1562 nm beyond its constituents' intersection:
    the trace CH4, N2O and CO2's band-rule poles fill their gaps. A representation of
    the traces with their bands folded, and N2O's ultraviolet absorption below
    320 nm, would widen it further.
- **Completing the curated records.** Values that were recorded in a record's notes
  before the format could hold them become datasets (TC-316): magnetic
  susceptibilities, higher refractivity virials, molar Kerr constants and Kerr virials,
  vibrational corrections to the polarizability, band-integrated cross sections.
- **High pressure.** The density dependence of the third-order response (the pair
  coefficients of χ⁽³⁾), and the higher refractivity virials, for fibres filled to
  tens of bar. The second refractivity virial of the gases that have none (O2, D2, N2O,
  CO, O3); OEP 4 records their omission in mixtures.

## solprop

optmat supplies the material properties used by solprop. The integration interface
and remaining responsibilities on the propagation side are:
- the recommended model of each gas;
- Lorentz–Lorenz density scaling;
- the process-aware Kerr and third-harmonic coefficients, through the model's exact
  time-domain form (see {doc}`nonlinear`). For solprop's SEP 6, which evaluates every
  four-wave quartet leg by leg, the manual now gives where the fits are supported: the
  legs of the data, each gas's resonance limits and the evidence on bounds for the legs
  beyond (TC-864 to TC-868); how to use them is solprop's decision. A public, tested
  record of it on the model, `HyperpolarizabilityModel.support`, so that solprop does
  not copy the numbers, is implemented by
  [OEP 5](oeps/archive/oep-0005-hyperpolarizability-support.md) (TC-870); the polyatomic
  gases' missing resonances stay with the papers we are looking for (below).
- the Raman responses of the gases, dry air and fused silica ({doc}`raman_response`),
  and the ionisation rates and parameters, through solprop's SEP 3. The snapshots
  `tests/data/solprop_raman.json`, `tests/data/solprop_raman_response.json` and
  `tests/data/solprop_ionisation.json` record historical solprop snapshots for the parity
  targets (TC-439, TC-690). The page's section *Using it in a propagation code* gives
  the cost of per-line widths along a density gradient and the low-rank tabulation that
  removes it.
- solprop's solids, so that it keeps no solid data either:
  - its window crystals and glass (CaF2, MgF2, BaF2, sapphire, N-BK7);
  - its silver and aluminium;
  - fused silica, with its Raman response and its n2 paired (TC-692, TC-696).

  The primary papers are now curated (TC-543 to TC-625): the crystals have their
  primary formulas made absolute ({doc}`crystals`), and silver, gold, copper and
  aluminium are table models ({doc}`metals`). Further sources for the metals (TC-626 to
  TC-648) took the tables to the soft x-ray range and aluminium to 0.8 nm-200 µm, and
  SCHOTT's data sheet and Weber's handbook brought N-BK7 and sapphire's extraordinary
  ray (TC-649 to TC-651). Each solid covered by those snapshots has a recommended model.

- the preparation on optmat's side is
  [OEP 1](oeps/archive/oep-0001-propagation-interface.md). That covers:
  - the declared surface, identities and flags;
  - overrides;
  - the polarizability derivative;
  - the EOS inverse and transport;
  - compositions and the air model;
  - the Raman groups;
  - the rate energy;
  - the silica bundle;
  - the fixtures, the first release and the integration page.

  The adapters themselves are solprop's.
- the finite-density mixture law of the gas index, for solprop's SEP 11, is
  [OEP 4](oeps/archive/oep-0004-mixture-index.md), with the mixture's density.

optmat will never depend on solprop.

## Solids: glasses and crystals

Glasses and crystals share one subpackage, `optmat.solid` (TC-318). **Fused silica is
the first** (see {doc}`silica`):
- its complex index from 25 nm to 125 µm (Franta et al.'s model refitted to the
  absolute refractometry);
- its n2 with the Raman part, the Raman fraction set by the measured Raman gain
  (TC-696);
- its Raman response: Hollenbeck & Cantrell's (the default), Lin & Agrawal's with its
  isotropic and anisotropic parts, and Blow & Wood's ({doc}`raman_response`).

The curated records already hold the data for what comes next (TC-325 to TC-392,
TC-543 to TC-651), for silica and the other 162 solids they cover:

- **Two-photon absorption:** beta at 193-800 nm of silica grades, the fluorides,
  sapphire, quartz and BBO, and a model of its dispersion with the band gap.
- **Ionisation and damage:**
  - multiphoton and field (Keldysh) ionisation rates and their inputs;
  - avalanche coefficients;
  - the laser-damage thresholds against pulse length.
- **Temperature:** a temperature-dependent recommended model of silica. Leviton & Frey's
  absolute indices and dn/dT (Corning 7980, 30-300 K) are curated, with their formula as
  a model (`leviton_2008/sellmeier`), and so are Matsuoka et al.'s (108-356 K) and
  Waxler & Cleek's (73-293 K) indices; Daimon & Masumura's dn/dT of CaF2.
- **More materials:** CaF2, MgF2, BaF2, sapphire and N-BK7 now have their formulas
  ({doc}`crystals`), silver, gold, copper and aluminium their tables ({doc}`metals`),
  and crystalline silicon its refitted universal dispersion model at 20 degrees C
  ({doc}`silicon`). LiF, YAG and KGW next:
  - their indices from primary measurements;
  - Franta's models where they exist (MgF2, Al2O3 films);
  - the n2 already curated (Adair, DeSalvo, Ensley, Patwardhan, Kabacinski and others).
- **Silicon, further:**
  - a temperature-dependent complex index: Franta et al.'s 2017 model needs its
    unprinted parameters from the authors (TC-1088); Frey et al.'s formula covers the
    window from 20 to 300 K;
  - the band-to-band absorption beyond 1.15 µm, which the model's tail cannot follow,
    and doped silicon's free carriers;
  - Herzinger et al.'s Si table (an AIP PAPS deposit) and Dash & Newman (1955), not
    yet in the library (TC-1083, TC-1085).
- **Catalogues:** vendor catalogues as comparison sources, and the spread between melts
  and vendors as part of the stated uncertainty.
- **Open questions (TC-393, TC-399):**
  - the grade dependence near the 2.73 µm OH band, where the recommended model departs
    from Malitson's glasses by up to 5e-5;
  - the electrostrictive part of nanosecond n2 data;
  - silica's Raman gain below 25 cm⁻¹, where Hollenbeck & Cantrell's shape is 1.6-2.4
    times Stolen et al.'s low-frequency law (TC-695).

## Crystals

For second-order nonlinear optics:
- the principal refractive indices of uniaxial and biaxial crystals, and their
  temperature dependence;
- the second-order susceptibility tensor ($d_{ij}$), and its dispersion through
  Miller's rule;
- the third-order response;
- the conventions for axes and tensor components that such data need, made explicit in
  the vocabulary.

## Infrastructure

- **Releases:** on PyPI, and documentation hosted online. The first tagged version, and
  the test that a clean install ships its data, were part of
  [OEP 1](oeps/archive/oep-0001-propagation-interface.md). Publishing and hosting remain
  open.
- **Read-only registries** (OEP 1, Q31): the registry functions return plain dicts of
  read-only models. Returning read-only mappings would change their return types.
- **Numerical regression across platforms** (OEP 1, Q32): gas Raman golden checks
  now compare canonical identity records with an arithmetic tolerance; public digests
  remain exact. Nonlinear refit checks compare predictions, the full objective and
  uncertainty-scaled parameter displacements. See {doc}`developers/testing`.
  Verification on other numerical platforms, including the exact digest checks for
  other model families (libm, CoolProp), remains part of setting up CI.
- **Continuous integration:** CI running the full gate on every pull request.

## Papers and data we are looking for

The papers below would improve the weakest parts of the models. Each is cited by a paper
already curated. If you have access to one, or know of better data, please open an
issue. Many of the papers requested earlier are now curated (TC-543 to TC-625); what
they left open is below.

**Gases:**

- **Vacuum-ultraviolet refractometry of the noble gases** between their first
  resonance and 140 nm:
  - G. I. Chashchina, V. I. Gladushchak and E. Ya. Shreider, Opt. Spectrosc. 24, 542
    (1968), argon;
  - G. I. Chashchina and E. Ya. Shreider, Opt. Spectrosc. 27, 79 (1969), krypton and
    xenon;
  - R. Abjean, A. Mehu and A. Johannin-Gilles, C. R. Acad. Sci. B 271 (1970).
- **The electronic absorption onsets and two-photon states of the polyatomic gases**
  (TC-865): CO2's weak absorption below its first band near 147.5 nm, SF6's onset (only
  absorption evidence is curated: none detected at Ityaksov et al.'s sampled wavelengths
  down to 197.7 nm, TC-869), and the lowest two-photon states of CO2, CH4, SF6, N2O and
  O3. Without them these gases' support records give no E2 and no derived band edge
  (TC-870).
- **An accuracy margin below each resonance** (OEP 5, Q8): how far below $E_1$ or
  $E_2$ the fitted polynomial stays within its stated accuracy. It would be published as
  information on the support record, never as a rule, if optmat had evidence for one;
  it has none today.
- **Two-photon ionisation and absorption of the noble gases** (the curation of
  TC-872 to TC-942; {doc}`two_photon`). The papers requested in the first round are
  curated; what they left open is:
  - **krypton's two-photon ionisation cross section in a vector-drawn or tabulated
    form** ([OEP 6](oeps/archive/oep-0006-two-photon-ionisation.md), Q2): every curve so far is a scan (L'Huillier, Tang and
    Lambropoulos 1989; Fink and Johnson 1990; McGuire 1981). Leads they cite include:
    - A. L'Huillier, L. Jönsson and G. Wendin, Phys. Rev. A 33, 3938 (1986);
    - B. Ritchie, Phys. Rev. A 16, 2080 (1977);
    - K. C. Kulander, Phys. Rev. A 36, 2726 (1987);
  - **xenon between its two thresholds** (Q3): only a scanned MQDT curve
    (Gangopadhyay et al. 1986);
  - **any measured two-photon ionisation cross section of argon or neon;**
  - **measured bound two-photon strengths of He, Ne and Ar, and the higher members of
    every series** (OEP 8);
  - **material inputs for the continuum's time-resolved two-photon loss** (TC-903;
    OEP 8, Q6; TC-1057): identify them with a solprop SEP for the propagation response;
    use a coupled optmat OEP for new material quantities, models or interfaces.
    Tier-3 TDSE runs support the pair-frequency rule for the helium and argon continua
    (TC-1049), but no published small-split amplitude is curated. Explicit approximate
    responses are possible candidates; their transient behaviour is not thereby
    qualified;
  - **Dakka et al.'s photoionisation rate factor** (OEP 6, Q7; TC-892).
- **Converged ab initio dispersion, and the magnetic susceptibilities, of krypton and
  xenon.** These would give Kr and Xe the theory backbone that He, Ne and Ar have.
- **Ultraviolet refractometry of SF6** below 458 nm, where its fit has no refractive
  data but Shelton's harmonics: Wilmouth & Sayres's scattering at 264-363 nm lies about
  2 % above the model, and Ityaksov et al.'s at 198-270 nm about 5 % (TC-712, TC-150).
- **Neon's refractivity virial** (TC-744): Achtermann et al.'s 1993 optical B_R of neon
  and Huot and Bose's B_ε, both as cited by Garberoglio and Harvey (2020). The package
  uses Orcutt and Cole's static value, about three times theirs and the ab initio one.
  Both papers are curated (TC-771) and confirm the comparison. The choice of B_R is the
  author's (OEP 4, Q18).
- **The mixtures' cross second virials and densities near ambient**, for OEP 4's Phase
  2. The citations are as printed by Tkaczuk et al. (2020) and Gernert (2013):
  - J. Brewer and G. W. Vaughn, J. Chem. Phys. 50, 2960 (1969);
  - M. L. Martin et al., Aust. J. Chem. 35, 1525 (1982);
  - W. F. Vogl and K. R. Hall, Physica 59, 529 (1972);
  - N. K. Kalfoglou and J. G. Miller, J. Phys. Chem. 71, 1256 (1967);
  - Vatter et al., Ber. Bunsenges. Phys. Chem. 100, 73 (1996);
  - Blancett et al.'s (1970) He–Ar and Kestin et al.'s (1981) He–Ne and Ne–Ar mixture
    densities.

  All are curated (TC-780 to TC-784). Still wanted: a complete copy of Blancett et al.
  (the library's scan clips pp. 82-88), and density data of the Ne–Ar mixtures and of
  air's pairs (N2–O2, Ar–N2, Ar–O2) near ambient, whose absence leaves those defaults
  provisional (TC-779).
- **Infrared refractivity of the molecular gases**, OEP 3's validators (TC-756). The
  five candidates of OEP 3's literature search are curated (TC-802 to TC-806):
  Rollefson and Havens's CH4 (1.7-15 um), Matsumoto's CO2 at 10.6 um and Tsankova et
  al.'s static CO are validators of the fits; Legay (1955) measured HCl alone;
  Mitrofanov et al.'s group delay of air at 3.7-4.15 um is recorded without a dataset,
  for want of a quantity for a group delay whose zero is not stated and of air as a
  registered material. Still wanted: Statescu, Phil. Mag. 30, 737 (1915), CO2's
  infrared index by the prism method. No measurement of CO or N2O between 2 and 12 um
  was found.
- **Hohm's measured polarizabilities against temperature:** U. Hohm, Dissertation, TU
  Braunschweig (1989), the α(ω, T) of Ar, H2, O2 and CH4 that Hohm & Kerl fit but do
  not tabulate; and which of Hohm's two 1994 xenon measurements, 0.31 % apart, is
  right (TC-606).
- **The molecules' polarizability away from their reference temperature:** an opt-in
  correction, linear in T - T_ref where a coefficient is measured (N2, O2, CO2; TC-730),
  with its wavelength dependence checked, and a model for H2's rotational average far
  from room temperature. OEP 1 exposes the reference temperature and the coefficients
  only (its §7.13).
- **Molecular photoabsorption:** the JILA Data Center Report No. 32 tabulations.
  Olney et al.'s S(−2) of N2O is 4.7 % above the one recomputed from Chan et al.'s
  published spectrum (TC-606).
- **Third-order response:**
  - the calibration of the Maryland group's absolute n2: the supplemental material of
    J. K. Wahlstrand, Y.-H. Cheng and H. M. Milchberg, Phys. Rev. Lett. 109, 113904
    (2012);
  - any measured hyperpolarizability of N2O.
- **Ionisation:**
  - static polarizabilities of Kr+, Xe+ and the molecular cations;
  - tunnelling structure coefficients (MO-ADK or weak-field asymptotic) for CH4, SF6
    and N2O: L. B. Madsen et al., Phys. Rev. A 85, 053404 (2012) and 87, 013406 (2013);
  - D. W. Turner, A. D. Baker, C. Baker and C. R. Brundle, Molecular Photoelectron
    Spectroscopy (Wiley, 1970): the original of the N2O and CO2 vertical energies;
  - a modern threshold of SF6 (Frost et al.'s 1967 value is of unstated kind);
  - the source of solprop's PPT Stark shift and dipole correction, which are not in
    Ilkov et al. 1992 (perhaps L. Bergé et al., Rep. Prog. Phys. 70, 1633 (2007));
  - the PPT originals: A. M. Perelomov, V. S. Popov and M. V. Terent'ev, Sov. Phys.
    JETP 23, 924 (1966); A. M. Perelomov and V. S. Popov, Sov. Phys. JETP 25, 336
    (1967), which would settle the factor of 2 in Ilkov et al.'s Eq. (8);
  - R. H. Garvey, C. H. Jackman and A. E. S. Green, Phys. Rev. A 12, 1144 (1975): the
    model potentials behind Batishchev et al.'s and Trinh et al.'s rates.
- **Raman:**
  - **for Q-branch line mixing:** relaxation-matrix (exponential-gap or
    energy-corrected sudden) parameters of N2, O2 and CO2 fitted to Q-branch spectra up
    to tens of atmospheres, and tabulated Q-branch widths or time-domain dephasing
    measured above 10 atm, where the branch has collapsed;
  - the vibrational O and S branches (their anisotropic elements beyond H2's and D2's),
    and the line shapes below the Dicke cutoff (TC-686, TC-688);
  - the rotational constants of N2O (a primary source; the Weber volume has none);
  - the Q-branch widths of CO2 and N2O, and N2O's depolarisation ratios;
  - Komasa et al.'s 2011 Supporting Information (the per-level uncertainties), and a
    vocabulary quantity for rovibrational level energies to hold their tables;
  - whether Murray & Javan's (1972) D2 widths are half or full widths (TC-586).
- **Ozone:** an independent zero-point vibrational correction to its polarizability, and
  the supplementary material of Christodouleas et al. (2010). The higher dispersion
  orders of its measured-absorption term (the GVD and beyond), which are repeated
  gradients of the grid rather than derivatives of the evaluated spline (TC-732).

**Solids** (the fused-silica models and the data curated for the next ones; TC-399,
TC-393, TC-325 onwards):

- **solprop's solids** (TC-546, TC-619, TC-648):
  - **silver and gold beyond 12.4 µm:** measured far-infrared data, such as M. A.
    Ordal et al., Appl. Opt. 26, 744 (1987) (gold); only Drude curves are curated;
  - **the vacuum ultraviolet of silver and gold** (40-200 nm), where the curated sources
    spread by 20-50 % in k (TC-646): a modern measurement of both n and k;
  - **MgF2's temperature dependence:** A. Feldman et al., NBS Tech. Note 993 (1979);
- **Silica's linear index:**
  - Popova et al. (1972), the IR data behind Kitamura's model;
  - the numbers behind Tan, Lemon & French's 2003 VUV constants (figures only);
  - J. H. Wray and J. T. Neu, J. Opt. Soc. Am. 59, 774 (1969): silica to 1100 K.
- **Franta's model:** the supplement of the 2025 Ta2O5 paper.
- **Raman gain:** a second absolute measurement of fused silica's Raman gain, bulk or
  in a pure-silica-core fibre (TC-696 rests on Stolen & Ippen's one), and its polarised
  and depolarised spectra, to test Lin & Agrawal's isotropic and anisotropic parts.
- **n2:** D. Milam, M. J. Weber and A. J. Glass, Appl. Phys. Lett. 31, 822 (1977), the
  fluoride crystals' n2; the other originals behind Weber, Milam & Smith's 1978 tables.
- **Damage and ionisation:** B. C. Stuart et al., J. Opt. Soc. Am. B 13, 459 (1996);
  L. Sudrie et al., Phys. Rev. Lett. 89, 186601 (2002); Keldysh, Sov. Phys. JETP 20,
  1307 (1965).

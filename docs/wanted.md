# Papers and data we are looking for

The papers below would improve the weakest parts of the models. Each is cited by a paper
already curated. If you have access to one, or know of better data, please use the
[source form](https://github.com/LupoLab/optmat/issues/new?template=offer_a_source.yml).
Cite the paper there rather than attaching it. Many of the papers requested earlier are
now curated (TC-543 to TC-625); what they left open is below. An item leaves this page when
it is curated, with its register entry
([#31](https://github.com/LupoLab/optmat/issues/31)).

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

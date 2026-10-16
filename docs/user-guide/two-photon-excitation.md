# Bound two-photon excitation

Start with {doc}`/two_photon` for the cross-section convention and model interface. This chapter describes the recorded evidence and its limits; source details are in the {doc}`/the-data/sources` catalogue.

## Bound two-photon excitation and the excited levels' ionisation

Below $I_p/2$, two photons excite the ground state to even-parity levels: $J = 0$ and 2
for equal, linearly polarised photons; $J = 1$ only with unequal photons or another
polarisation. {mod}`optmat.gas.two_photon_bound` states what is curated about them, per
gas (OEP 6, Phase 2a: curation only, TC-945).

```python
from optmat.gas.two_photon_bound import two_photon_bound_levels

xenon = two_photon_bound_levels("Xe")
xenon.levels  # every even J = 0, 1, 2 level below I_p, from NIST
level = xenon.level("5p5.(2P*<3/2>).6p", "2[1/2]", 0)  # NIST's strings
level.wavelength  # each of two equal photons, m: 249.63 nm
xenon.strengths_of(level)  # calculated and measured, in one convention
xenon.photoionisation_of(level)  # the excited level's ionisation, m^2
xenon.assessment, xenon.assessment_reason  # "checked" or "unresolved", and why
```

**The levels** are every even level with $J = 0, 1, 2$ below the ionisation energy in
the NIST lists held as returned; helium's singlets only.

**One convention for the strengths** (TC-952). A calculation's two-photon strength
$\delta^{TP}$ is per level, summed over $M_J$, for parallel linear polarisation and
equal photons, in atomic units. A frequency-integrated cross section follows from it:

$$
\sigma_0^{(2)} = 2\pi^3\alpha^2 a_0^4\,\omega_\mathrm{au}^2\,\delta^{TP}
= 2.589\times10^{-44}\ \mathrm{m^4}\ \omega_\mathrm{au}^2\,\delta^{TP},
$$

for coherent light, with the line profile normalised in the two-photon angular
frequency. The rate is $W = \sigma_0^{(2)} g(2\omega)\Phi^2$. This is the formula of
every curated measurement (Drag et al., Shekhtman et al., Khambatta et al., Pindzola et
al.), and the Rabi frequencies of Drag et al., Dakka et al. and Kumar et al. follow from
it. A measured $\sigma_0^{(2)}$ is held corrected for the photon statistics its source
used.

- **Per level, not per component:** for a $J = 2$ level, the strength of one tensor
  component, which Hättig et al., Sundholm et al. and Hettema et al. print, is a fifth of
  the level's. Henriksson et al. and Moccia and Rizzo print the level's. Xenon decides it.
  Within one experiment, where calibration and photon statistics cancel, the measured
  $J = 0$ to $J = 2$ ratios follow the per-level reading: Drag et al.'s 6p′ ratio is 1.38,
  where per level gives 1.41 and one component 7.0.
- **Placement:** each assigned strength's $\sigma_0^{(2)}$ is evaluated at its level's
  NIST photon frequency. A calculation's own energy can be off NIST's by 0.1 eV (Xe, Kr)
  to 1.4 eV (Ne).

**What the records show:**
- **Helium, checked:** Hättig et al.'s coupled-cluster (CCSD) 1s2s 1S strength agrees
  within 1 to 4 % with:
  - Drake's exact two-photon decay rate;
  - Kumar et al.'s and Kaprálová-Žďánská's two-photon couplings;
  - Bassi and Baluja's partial closure (TC-952).

  The simpler calculations (CCS, CC2, Henriksson et al.'s TDDHF) are 0.85 to 0.94 of it.
  There is no measured cross section.
- **Xenon, checked:** Henriksson et al.'s relativistic strengths are 1.30, 0.83 and 0.98
  of Kröll and Bischel's measurements at the 6p levels, and 1.17 and 1.19 of Drag et
  al.'s at 6p′.
- **Krypton, unresolved:** the measurements and calculations disagree by up to 94 times
  (Shekhtman et al. against Khambatta et al. at 6p[3/2]₂). Henriksson et al. are 0.16 of
  Dakka et al. at 5p[3/2]₂.
- **Argon, unresolved:** its only measurement, McCann et al.'s, is uncertain by a factor
  of about 4.5 within its own laser widths, which spans the factor of 5 at issue.
- **Neon, unresolved:** calculations only.

Values that cannot be tied to one level are held unassigned, with the reason: LS-coupled
$^1S$ strengths, cross sections averaged over $J$.

**No evaluator.** A line profile, the excited population and the (2+1) ionisation that
follows need the theory OEP 6 asks for in Q9 (Phase 2b, TC-945). The record supplies
their inputs.

**Missing:**
- measured bound strengths of He, Ne and Ar;
- a resolution of argon's and krypton's measurements;
- the photoionisation of most excited levels;
- the higher members of every series, which have positions but no strengths.

These are listed among the papers the curation requests (*Papers and data we are
looking for*, on the roadmap).


See {doc}`/oeps/archive/oep-0006-two-photon-ionisation` for the completed proposal and its phase reports.

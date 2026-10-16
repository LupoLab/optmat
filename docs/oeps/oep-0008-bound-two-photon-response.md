# OEP 8 — The bound two-photon resonances of the noble gases for propagation

| | |
|---|---|
| **OEP** | 8 |
| **Title** | The bound two-photon resonances of the noble gases for propagation |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 (2026-10-09); GPT-6 (review revisions, 2026-10-12) |
| **Status** | Draft |
| **Review status** | Reviewed |
| **Reviewers** | John Travers (human); GPT-6 (AI) |
| **Type** | Feature |
| **Created** | 2026-10-09 |
| **Resolution** | — |
| **Requires** | [OEP 1](archive/oep-0001-propagation-interface.md): model identity, digests and the declared public surface; [OEP 5](archive/oep-0005-hyperpolarizability-support.md): the fitted $\gamma$ and its resonance limits; [OEP 6](archive/oep-0006-two-photon-ionisation.md): the continuum model, the bound-level records and the pair-spectral identity |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | — |
| **Manual pages** | — |

## Part I — Proposal

### 1. Abstract

Below half its ionisation energy a noble gas absorbs two ultraviolet photons only at its
even-parity bound levels: He 1s2s, Xe 6p and the like. The excited atom may then absorb
a third photon and ionise. OEP 6 curated those levels and their two-photon strengths
(Phase 2a), but deferred the response a propagation code needs (its question Q9). The
reason was that no derivation yet tied a time-dependent response to the curated
amplitudes, its units, its residues and its energy balance. This OEP supplies that
derivation and proposes what optmat should publish from it.

- **The model.** An effective-level Bloch model, obtained by eliminating the odd
  intermediate states. Its two-photon coupling is $-\tfrac12M_f(E^+)^2$, and its
  conservative polarisation is $-\partial\langle H_R\rangle/\partial E^-$, with
  photoionisation loss accounted for separately. The local model obeys the energy
  identity in section 4.4. General frequency-dependent dynamics are deferred (Q8).
- **What optmat publishes, per qualified level:**
  - the residue $\delta^{TP}$, with its uncertainty and polarisation factors;
  - an empirical estimate of the error of the pair-frequency form;
  - the excited level's photoionisation cross section;
  - where curated, the two-colour amplitude.

  The propagation code owns the resulting polarisation, electron accounting,
  propagation energy balance and integration. The derivation here checks the material
  coefficients and their conventions; it is not an optmat propagation interface.
- **Coefficient qualification proposed:**
  - He 1s2s strengths and photoionisation; strong-field dynamics additionally need
    the light shifts to be curated in Phase 1;
  - Xe 6p and 6p′, for weak-field absorption and populations. Their excited-level light
    shifts are not curated, and the dispersive part of their response fails Lehmberg
    et al.'s measured $n_2$ by a factor of 2.2 to 2.8;
  - Ar 4p, the $J=0$ and $J=2$ levels, for weak-field absorption and populations within
    the limits of section 6, on calculated strengths only (Q7). Their radiative widths
    are measured, and their photoionisation at the drive is calculated. Their light
    shifts and pressure broadening are not curated.
- **Refused:** Ne and Kr.

**Success** is judged by three reproductions and one honesty test:
- the weak-field limit reproduces OEP 6's pair-spectral yield;
- a published broadband He 1s2s calculation (28 % transfer) is reproduced;
- the specified model closes its energy budget to a converged numerical tolerance;
- every unqualified case is refused.

### 2. Motivation and scope

**The need.** Resonant dispersive-wave emission and soliton self-compression in
gas-filled hollow-core fibres reach the vacuum ultraviolet. Their spectra cross the
two-photon resonances of the filling gas:
- He 1s2s at 120.3 nm;
- Xe 6p[1/2]₀ at 249.6 nm;
- the 6p and 6p′ levels of Xe between 222 and 256 nm.

**What a resonance does to the pulse.**
- It depletes the pulse.
- It leaves excited atoms that a third photon ionises: the (2+1) channel, which a
  plasma model must count.
- Near the resonance, it changes the Kerr response.

**Why linear absorption is not enough.** A femtosecond pulse at these intensities
drives the level coherently. In Kumar et al.'s case [1] the peak two-photon Rabi
frequency exceeds both the level's decay rate and the pulse's bandwidth, so rate
equations fail (section 4.5).

**The paired solprop proposal** (solprop's two-photon ionisation proposal) asked optmat
for:
- the continuum cross section, which OEP 6 supplies;
- the bound resonances, deferred until optmat can say what they are.

**In scope:**
- the theory needed to establish the conventions and qualification of the curated
  two-photon amplitudes and material coefficients;
- the coefficients optmat publishes for it, their uncertainties and qualification;
- the reproducible checks;
- the refusals.

The first dynamical scope is the local one-beam model (one field, one polarisation
state across its supported spectrum), with its spectral and field-strength conditions
checked. A pulse may be broad compared with the final-state linewidth while narrow
enough relative to intermediate-state detunings. The two-colour helium kernel is
restricted to perturbative on-shell yields, not a general pump-plus-ultraviolet
evolution model (Q8, TC-1055).

**Out of scope:**
- **The propagation response.** solprop owns the equations and implementation for
  the resulting polarisation, electron accounting, plasma current and propagation
  energy balance, as well as integration. optmat supplies material properties and
  basic data, as it supplies ionisation rates without a plasma polarisation and Raman
  responses without the resulting field polarisation (TC-1057). The single-atom
  derivations and internal fixtures here check material conventions and coefficients;
  they do not transfer ownership of the propagation response to optmat.
- **The continuum's time-resolved loss.** No published small-split continuum amplitude
  is curated. Tier-3 TDSE runs now support the pair form for helium and argon (section
  4.2), but its propagation response belongs in a solprop SEP. A coupled optmat OEP
  is needed only if that SEP requires new material quantities, models or interfaces
  beyond OEP 6. This does not extend OEP 8's phases (Q6, TC-1057).
- **The plasma kinetics after the pulse.**
- **Higher members of the even series,** without curated strengths.
- **Molecules.**
- **The ionic fine structure of the (2+1) electron,** a matter for the plasma model.

#### 2.1 The initial experiment and the simplest useful test

**The maintainer's use case (2026-10-12; TC-1056).** A 6 fs, 800 nm pump
self-compresses in argon and emits a weaker VUV resonant dispersive wave (RDW). When
pressure and phase matching tune the RDW towards the two-photon ionisation threshold,
the RDW disappears. Pump-driven ionisation disrupting phase matching can explain the
observations. The scientific question is whether VUV two-photon absorption or
ionisation can also explain them, and which contribution dominates. Helium is a
second target. These are the maintainer's reported context and competing hypotheses,
not a curated experimental dataset or a conclusion of this OEP.

Distinguish departure from a phase-matching curve as the pump blue-shifts from the
shortest observable RDW across optimised operating conditions. Evidence for the
first mechanism does not by itself establish the cause of the latter limit.
Frequency-selective loss can also shift the measured RDW centroid; distinguish that
centroid, the blue edge and the detection limit when testing the hypotheses. Reconcile
the input duration for each experimental or simulated dataset before defining a
fixture; the 6 fs context above is not a universal input specification (TC-1057).

**Continuum first.** The manual's NIST-derived thresholds ({doc}`../two_photon`;
TC-943) give 157.34 nm for two equal photons in Ar and 100.85 nm in He. These are
OEP 6's continuum thresholds. OEP 8's first bound levels are elsewhere: Ar 4p at
about 185–192 nm and He 1s2s at 120.28 nm. Include them when the actual pair spectrum
overlaps them. The higher bound series approaching threshold still lack curated
strengths; a null result for the retained low levels does not exclude all bound
absorption. The continuum diagnostic need not wait for this OEP's implementation.

At the Ar threshold, an 800 nm photon plus one 157.34 nm photon carries about
9.43 eV, below the 15.76 eV ground-state ionisation energy. Thus that particular
mixed two-photon channel is closed for field-free atoms; the near-threshold direct
channel uses two VUV photons. This does not exclude pump dressing, higher-order
pump-assisted processes, excited-state ionisation, or other spectral components.
The intense pump can shift levels even when this mixed channel is closed.

**Recommended solprop investigation, not a new optmat propagator:**

1. Generate matched pressure scans with pump-driven ionisation disabled and enabled,
   holding the other input parameters fixed. Save the complex VUV field versus
   propagation distance, its transverse distribution or mode normalisation, the local
   neutral density and the pump-induced electron density. The input's 6 fs is not a
   proxy for the RDW duration, its bandwidth or its intensity at generation.
2. On these prescribed fields, evaluate OEP 6's perturbative pair-spectral probability
   and energy loss. With $E_V^+$ the declared VUV band and $S_{2,V}$ the transform of
   $(E_V^+)^2$, use the existing kernel $K$ of OEP 6 section 4:

   $$
   P_2(z)=\frac{1}{2\pi}\int_0^\infty K(\Omega)|S_{2,V}(z,\Omega)|^2d\Omega,\qquad
   \Delta\mathcal E_2(z)=\frac{1}{2\pi}\int_0^\infty\hbar\Omega K(\Omega)
   |S_{2,V}(z,\Omega)|^2d\Omega .
   $$

   These are per-atom diagnostics within that approximation. Integrating
   $N_g\Delta\mathcal E_2$ over transverse area and distance estimates the lost pulse
   energy on the prescribed field. Likewise $N_gP_2$ estimates the added electron
   density only when depletion is negligible and the neutral population can be treated
   as fixed during the VUV pulse. A time-varying pump-depleted population needs a
   specified dynamical model; it cannot simply be inserted into an integrated yield.
   Track unsupported spectral weight and vary the band split and available source
   envelope. In particular, the Ar model has a gap immediately above threshold and
   unresolved fine-structure detail between 157.34 and 155.59 nm. Neither a gap nor an
   unqualified contribution is zero.
3. Where the whole relevant band has a slowly varying supported cross section, the
   cheaper estimate $\beta=2N_g\sigma^{(2)}/(\hbar\omega_V)$ and
   $dI_V/dz=-\beta I_V^2$ checks the scale of loss. Its dimensionless scale is
   $\beta I_VL$, with the actual intensity and interaction length. A centroid
   coefficient is not reliable when the pair spectrum straddles the threshold.
4. If the diagnostic indicates a potentially relevant effect, compare four propagation
   cases: neither channel, pump ionisation only, VUV two-photon ionisation only, and
   both. The two-photon channel includes its field-energy loss and its electron
   contribution; an additional loss-only run can isolate direct absorption from plasma
   feedback. Avoid counting the VUV band again in a pump strong-field rate. Record any
   assumed continuum timing and spectral-loss rule: OEP 6's integrated constraints do
   not supply a unique time-dependent current (Q6). A phenomenological closure in
   solprop is an exploratory model, not a newly qualified optmat response.

The post-processing step tests the scale and spectral onset; it cannot by itself
predict disappearance or establish the dominant mechanism. Small direct loss also
need not imply a small phase-matching change from the added electrons. If alternative
admissible continuum timing models change the conclusion, or pump dressing is large,
that is the evidence needed to prioritise the deferred dynamics. No particular outcome
is asserted without the simulated fields and experimental geometry.

**The continuum follow-up belongs in a solprop SEP (TC-1057).** Its deliverable is
the propagation response: continuum polarisation, electron accounting, plasma coupling
and energy-balance tests. optmat's contribution is the cross sections and other
material properties, their provenance, uncertainty and validity limits. First identify
which additional material inputs the SEP actually needs; a coupled optmat OEP would
cover only those additions, with material validation fixtures for parity.

An effective continuum driven by $(E_V^+)^2$, with coupling strengths matched to
OEP 6's cross sections, is a candidate for that SEP to derive and validate. The absence
of published nondegenerate amplitudes does not forbid such an approximate model.
However, matching integrated yields does not uniquely fix the transient response,
validate electron production during the pulse, or resolve missing near-threshold data.
The candidate is not an approved or qualified model. Its detailed design and tests
belong in the SEP; new material-data requirements may motivate the coupled OEP.

### 3. Context

**What optmat supplies today** ({doc}`../two_photon`):
- `optmat.gas.two_photon`: the degenerate continuum cross section $\sigma^{(2)}(\omega)$
  of He, Ne, Ar and Xe, with both thresholds, labelled regions, an envelope and an
  identity (OEP 6 Phase 1; TC-943, TC-944).
- `optmat.gas.two_photon_bound`: the even $J = 0, 1, 2$ levels below threshold from
  NIST, with every curated strength in one convention and the excited levels'
  photoionisation cross sections (OEP 6 Phase 2a; TC-952, TC-953).
  - **The convention** (TC-952): $\delta^{TP}$ per level, summed over $M_J$, parallel
    linear polarisation, equal photons, in atomic units, with
    $\sigma_0^{(2)}=2\pi^3\alpha^2a_0^4\omega_\mathrm{au}^2\delta^{TP}$.
  - **The assessments:** He and Xe checked; Ne, Ar and Kr unresolved.
- `optmat.gas.hyperpolarizability`: the fitted Taylor-series $\gamma$ with its support,
  below the band edge $E_2/2$ (OEP 5).
- **The near-infrared NIST lines of Ne, Ar, Kr and Xe** (TC-983): the transition
  probabilities among excited levels. They feed the sums over intermediate states and
  the radiative widths (section 5).

**What OEP 6 deferred:**
- its Phase 2b, the complex evaluator, until Q9 is answered (TC-945);
- its optional Phase 3, consistency of a real part with $\gamma$ (Q8);
- the widths and their conventions (Q7).

Its own fallback is to report the evaluator as deferred to a later proposal. This is
that proposal. It answers Q9 for the bound levels. It takes up Q7 and Q8 where they
concern the bound levels. OEP 6 is Final with Phases 1 and 2a implemented (TC-1053).

**The roadmap item answered:** OEP 6's Phase 2 evaluator, and the "bound resonances"
deferred by the paired solprop proposal.

### 4. Theory

SI throughout, with the manual's conventions ({doc}`../conventions`):
- $E=E^++E^-$, with $E^+(t)=\int_0^\infty\frac{d\omega}{2\pi}\tilde E(\omega)e^{-i\omega t}$;
- a monochromatic $E_0\cos\omega t$ has $E^+=(E_0/2)e^{-i\omega t}$;
- $\gamma$ is the Taylor-series hyperpolarizability, so $\mu_\omega=\gamma E_0^3/8$;
- $\Phi=I/\hbar\omega$, and $W=\sigma_0^{(2)}g(2\omega)\Phi^2$ for a bound line of unit-area
  profile $g$ in the two-photon angular frequency.

$N$ is the initial atomic number density. These equations use the dilute-gas
approximation in which the field acting on an atom equals the propagated field and
$I=2\varepsilon_0c|E^+|^2$. Finite-density local-field corrections are not established
by this derivation. Coherence rates describe homogeneous decay: Doppler broadening
requires an average over velocity-dependent detunings, not insertion of a Gaussian
width as $\gamma^\phi$.

The numerical results below were reported during drafting from source records and
NIST data. They are targets for reproducible Phase 1 fixtures, not tests already in
the public tree. Section 8 distinguishes the independent review checks (TC-1054).

#### 4.1 The exact second-order amplitude

**The amplitude.** With $H=H_0-\mathbf d\cdot\mathbf E$, take second order and keep the
$E^+E^+$ terms (the two-photon rotating-wave approximation). Symmetrising over the two
photons gives the final-state amplitude:

$$
c_f(\infty)=\frac{i}{2\hbar}\int_0^{\Omega_f}\frac{d\omega}{2\pi}\,
\tilde E_\alpha(\Omega_f-\omega)\tilde E_\beta(\omega)\,M^f_{\alpha\beta}(\Omega_f-\omega,\omega),
\tag{1}
$$

**The kernel:**

$$
M^f_{\alpha\beta}(\omega_a,\omega_b)=\sum_n\Big[\frac{d^\alpha_{fn}d^\beta_{ng}}{E_{ng}-\hbar\omega_b}
+\frac{d^\beta_{fn}d^\alpha_{ng}}{E_{ng}-\hbar\omega_a}\Big].
\tag{2}
$$

**The denominators** carry $-i0$, or $-i\hbar\Gamma_n/2$ for an intermediate level of
width $\Gamma_n$. That matters only when a pole lies inside the spectrum. There the pair
form fails, and phase or amplitude shaping can raise the transfer above the
transform-limited value, as Dudovich et al. measured [20]. This OEP excludes that regime
(section 4.3).

**The omitted $E^-E^+$ terms** reach $f$ only if the field holds two components
separated by $\Omega_f$.

**Check 1, the monochromatic limit.** It gives
$W=\frac{\pi}{32\hbar^2}E_0^4|M_{zz}|^2g(\Omega)$. Hence:

$$
\sigma_0^{(2)}=\frac{\pi\omega^2|M_{zz}|^2}{8\varepsilon_0^2c^2}
=2\pi^3\alpha^2a_0^4\omega_\mathrm{au}^2\delta,
$$

with $\delta=|M_{zz}(m_f=0)|^2$ in atomic units. This is TC-952's per-level convention,
reached independently. Here $M$ has SI units C m²/V, $\sigma_0^{(2)}$ has units
m⁴, and $g$ has units s; the au strength is converted before an SI calculation.

**Check 2, the continuum.** If $M$ does not depend on the split, (1) summed over an
energy-normalised continuum gives OEP 6 §4's identity:
- $P=\frac1{2\pi}\int K|S_2|^2d\Omega$;
- with $S_2(\Omega)=\int e^{i\Omega t}(E^+)^2dt$.

**Split variables.** Write $\omega_{a,b}=\Omega/2\pm\nu$ and $D_n=E_{ng}/\hbar-\Omega/2$.
The kernel separates into a part even in $\nu$, from the symmetric product $S^n$ of the
two dipoles, and a part odd in $\nu$, from their antisymmetric product $A^n$:

$$
\hbar M^f_{\alpha\beta}=\sum_n\Big[S^n_{\alpha\beta}\frac{2D_n}{D_n^2-\nu^2}+A^n_{\alpha\beta}\frac{2\nu}{D_n^2-\nu^2}\Big].
\tag{3}
$$

**Tensor structure from the $J_g=0$ ground state:**
- **Rank 0** (the trace of $S$) reaches $J_f=0$.
- **Rank 2** (the traceless part of $S$) reaches $J_f=2$.
- **Rank 1** ($A$) reaches $J_f=1$.

**For one field of fixed polarisation,** $\epsilon_a=\epsilon_b$, so the $A$ term
vanishes at every split: no $J=1$ level is excited.
- **Linear parallel:** $J=0$, and $J=2$ with $M=0$, at strength $\delta$ per level.
- **Co-rotating circular:** $J=0$ is 0; $J=2$ is 3/2 of the linear value (TC-928).
- **Unequal polarisations** (two beams, or a polarisation that varies across the
  spectrum) open $J=1$, at order $(\nu/D)^2$ in probability.
- **Crossed linear beams** give no $J=0$.

#### 4.2 When the pair-frequency form holds

**When it is exact.** Within the two-photon rotating-wave approximation, Eq. (1) reduces to OEP 6's pair form,
$c_f=\frac{i}{2\hbar}M^f(\Omega_f)S_2(\Omega_f)$, for arbitrary spectral amplitudes when:
- one polarisation state holds across the spectrum;
- and $M^f$ does not depend on $\nu$ over the pair spectrum's support.

For a physical kernel with finite intermediate distances $D_n$, this generally
requires an approximation; agreement for a particular pulse is not an equivalence
for all fields.

**The error.** Expand $M(\nu)=M(0)[1+\kappa\nu^2+O(\nu^4)]$:

$$
\kappa=\frac{\sum_n w_n/D_n^2}{\sum_n w_n},\qquad w_n=S^n_{zz}/D_n .
\tag{4}
$$

For a Gaussian spectrum of rms width $\sigma$ (of $|\tilde E|^2$, in photon frequency)
and group-delay dispersion $\phi''$:

$$
\langle\nu^2\rangle=\frac{\sigma^2}{1-2i\phi''\sigma^2},\qquad
\frac{P-P_0}{P_0}=2\,\mathrm{Re}(\kappa\langle\nu^2\rangle)+O(\sigma^4),\qquad
\left|\frac{P-P_0}{P_0}\right|\lesssim2|\kappa|\sigma^2 .
\tag{5}
$$

- **If every $w_n$ has one sign,** $0<\kappa\le D_\mathrm{min}^{-2}$.
- **If they cancel,** $\kappa$ is unbounded, and it can be negative: OEP 6's
  "cancellations" warning, made quantitative. A negative $\kappa$ means bandwidth lowers
  the transfer. Hydrogen's exact amplitudes show both signs: $\kappa=36.8$ hartree⁻² for
  1s–2s and $-16.3$ hartree⁻² for 1s–6s (TC-1034). The rule below estimates $|\kappa|$.
- **For real $\kappa$ and quadratic spectral phase,** chirp reduces the magnitude
  of this leading Gaussian correction. This is not a bound for arbitrary shaping,
  complex $\kappa$, cancellations of the pair amplitude, or higher-order terms.
- **For a mixed pair** (pump plus ultraviolet), $\nu$ is half the colours' separation,
  (5) does not apply, and the two-variable $M$ is required.

**Tested against curated data:**

| level | evidence | $D_\mathrm{eff}=\lvert\kappa\rvert^{-1/2}$ against $D_1$ |
|---|---|---|
| He 1s2s | Drake et al.'s variational two-photon decay spectrum [2], against Hättig et al.'s small-split CCSD value [3] (1.00073 against 1.00075) | 0.75–0.86 (a second, opposite-signed effective pole) |
| H 1s–2s, 3s, 6s | the exact hydrogen amplitude against the split (Quattropani et al., Drachman et al.; TC-1034) | 1s–2s: 0.88; 1s–3s: 2.4 (a cancelling path); 1s–6s: $\kappa<0$ |
| He 1s2s, two colours | Kumar et al.'s two-colour coupling [1] at 19.06 + 1.54 eV: $\lvert M/M_\mathrm{deg}\rvert^2=13.83$, against Drake's 14.04–14.11 | a check of the full two-variable $M$, 1.5–2 % |
| Ne 3p ¹S, ¹D | Hättig et al. at $\nu=0$ and at $\omega_2=10.49$ eV [3] | 0.89–0.99 |
| Ar 4p ¹S, ¹D | Hättig et al. [3]; Sundholm et al. [4] | 0.93–1.01 |
| Xe 6p, Kr 5p | single-path estimates from measured oscillator strengths, checked against Shekhtman et al.'s [5] and Gornik et al.'s [6] | at least 0.77, assuming the second path's sign |

Here $D_1=E_1/\hbar-\Omega_f/2$ is an angular frequency, with $E_1$ the gas's
lowest odd $J=1$ level; the table below gives $\hbar D_1$ in eV. For negative
$\kappa$, define $D_\mathrm{eff}=|\kappa|^{-1/2}$.

**The adopted empirical rule** is $D_\mathrm{eff}=0.75\,D_1$ (Q3, TC-984). The
cases above support this estimate under their stated assumptions; they do not prove a
universal bound on cancellations or the $\nu^4$ remainder. The estimated rms width
for a leading fractional error $\varepsilon$ is
$\sigma_\mathrm{max}=0.75\,D_1\sqrt{\varepsilon/2}$.

| final levels | $\hbar D_1$ (eV) | estimated shortest transform-limited intensity FWHM for 1 % (fs) |
|---|---|---|
| He 1s2s | 10.9 | 1.3 |
| Ne 3p | 7.2–7.4 | 2.0 |
| Ar 4p | 4.9–5.1 | 2.9–3.0 |
| Kr 5p | 3.9–4.3 | 3.4–3.7 |
| Xe 6p; 7p and 6p′ | 3.5–3.6; 2.9–3.0 | 4.1–4.2; 4.9–5.1 |

These estimates suggest a small split correction for smooth Gaussian pulses longer
than about 5 fs. Qualification still requires the actual spectrum, intermediate-pole
distances, cancellations and field-strength limit to be checked. A shaped pulse with
$S_2(\Omega_f)$ near zero needs an absolute amplitude error, not a relative bound.

**The continuum is different.** Its one curated nondegenerate value is McCann et al.'s
Xe 118 + 355 nm [7]. It lies past the 5d[3/2]₁ pole, where the pair form fails by
$2.7\times10^4$. It shows the size of mixed-pair errors and tests nothing at small
splits.

**Tier-3 evidence at small splits** comes from the maintainer's tRecX runs (TC-1048,
TC-1049). With a one-colour reference at $\Omega/2$, $P_\mathrm{cross}/P_\mathrm{ref} =
4|M(\omega_1,\omega_2)/M(\Omega/2,\Omega/2)|^2$ directly.
- **Helium**, two-electron haCC at $\Omega/2$ = 14, 17 and 20 eV: $\kappa>0$ everywhere,
  and $D_\mathrm{eff}$ = 1.1–2.2 $D_1$. The one-electron model gives 1.03–1.8 $D_1$.
- **Argon**, one electron at 8.5 eV: $D_\mathrm{eff}\approx1.9\,D_1$.

These runs support the 0.75 $D_1$ estimate at the sampled continuum energies. They do
not qualify the continuum's transient response. A solprop SEP must derive and validate
that response; missing published amplitudes do not preclude an explicit approximate
model, but its additional material inputs need separate qualification (Q6, TC-1057).

#### 4.3 The effective-level dynamics

**Adiabatic elimination.** Remove the odd intermediate states. This needs two
conditions:
- **spectral:** no spectral component lies within a few $\sigma$ of an odd resonance,
  and $|\nu|\ll D$ for the local pair approximation (section 4.2);
- **field strength:** $|d_{ng}E_0|/(2\hbar D_n)\ll1$ for every intermediate $n$. Its
  error is the admixture of the intermediate states, $2(d E_0/2\hbar D)^2$ in
  probability in the cited three-level example, not a general multilevel bound.
  The excited-to-intermediate couplings must also be small compared with their
  detunings; a test of $d_{ng}$ alone is insufficient. In Kaufman et al.'s three-level
  case [21] that error is 3.1 % of the transfer, against 0.2 % from the spectral estimate (5). It grows with peak intensity at
  fixed bandwidth (TC-990).

The basis is then the ground level $g$ and the even levels $f$ inside the pair spectrum:

$$
H_{gg}=-\alpha_g^R|E^+|^2,\quad H_{ff}=\hbar\Omega_f-\alpha_f^R|E^+|^2,\quad
H_{fg}=-\tfrac12\mathcal Q_f(t),\quad \mathcal Q_f=M_fE^+(t)^2 .
\tag{6}
$$

- **The Hamiltonian in (6) is Hermitian**, denoted $H_R$ below, with
  $H_{gf}=H_{fg}^*$ and $\alpha_j^R=\mathrm{Re}\,\alpha_j(\omega_c)$ at a specified
  carrier. $M_f$ and these real shifts are constant coefficients in this local model.
  Loss enters the density equation separately, once. A frequency-dependent kernel or
  polarizability is nonlocal in time: substituting a convolution into (6) does not
  establish (8) or the energy identity without a corresponding functional derivation
  (deferred by Q8).
- **Ionisation of the excited level:** $\mathrm{Im}\,\alpha_f=\varepsilon_0c\,\sigma^{PI}_f/\omega$,
  so its rate is $\Gamma^{PI}_f=2\alpha_f^I|E^+|^2/\hbar=\sigma^{PI}_f\Phi$,
  where $\alpha_f^I=\mathrm{Im}\,\alpha_f$. The ground-state one-photon loss is
  excluded here. Putting this imaginary part in $H_R$ as well would count loss twice.
- **Kumar et al.'s Eq. 7 [1] is (6),** with $\Omega_{ge}=\sqrt{2\sigma_0/\pi}\,\Phi$.
  Kaufman et al. [21] and Tóth et al. have the same form (TC-990, TC-991).
- **A naming trap:** Kaufman et al. and Tóth et al. call their "two-photon Rabi frequency" is the off-diagonal element
  itself, $|H_{fg}|/\hbar$. That is half of this OEP's $\Omega_R$ and of Kumar et al.'s
  $\Omega_{ge}$.

**The population and coherence equations** for independent loss channels, with
$\rho_{ab}=\langle c_a c_b^*\rangle$ and no direct effective coupling between excited levels:

$$
\begin{aligned}
\dot\rho_{ff}&=\hbar^{-1}\mathrm{Im}[\mathcal Q_f^*\rho_{fg}]-(\Gamma^{PI}_f+\Gamma^{rad}_f)\rho_{ff},\qquad
\dot\rho_{gg}=-\textstyle\sum_f\hbar^{-1}\mathrm{Im}[\mathcal Q_f^*\rho_{fg}],\\
\dot\rho_{fg}&=-i\big[\Omega_f-(\alpha_f^R-\alpha_g^R)|E^+|^2/\hbar\big]\rho_{fg}
+\tfrac{i}{2\hbar}\mathcal Q_f(\rho_{gg}-\rho_{ff})
-\Gamma^{coh}_f\rho_{fg}-\tfrac{i}{2\hbar}\textstyle\sum_{f'\ne f}\mathcal Q_{f'}\rho_{ff'} .
\end{aligned}
\tag{7}
$$

- **Coherence decay:** $\Gamma^{coh}_f=\tfrac12(\Gamma^{PI}_f+\Gamma^{rad}_f)+\gamma^\phi_f$.
  Half of the population decay goes to the coherence; pure (collisional) dephasing
  $\gamma^\phi_f$ is separate (OEP 6, Q7).

For $f\ne f'$, the missing equation is fixed by the same commutator:

$$
\dot\rho_{ff'}=-\frac{i}{\hbar}(H_{ff}-H_{f'f'})\rho_{ff'}
+\frac{i}{2\hbar}\left(\mathcal Q_f\rho_{gf'}-\rho_{fg}\mathcal Q_{f'}^*\right)
-\left[\frac{\Gamma_f+\Gamma_{f'}}2+\gamma_f^\phi+\gamma_{f'}^\phi\right]\rho_{ff'},
\qquad \Gamma_f=\Gamma_f^{PI}+\Gamma_f^{rad}.
$$

This dephasing convention is the Lindblad choice
$L_f=\sqrt{2\gamma_f^\phi}|f\rangle\langle f|$, with independent level noise.
Correlated collision shifts or transfer between retained levels require a specified
relaxation operator; single-line widths alone do not determine it. The sign of the
last term of (7) follows from $-i[H_R,\rho]/\hbar$ (TC-1054).

- **The coherences between excited levels** $\rho_{ff'}$ are needed when two levels lie
  inside one pair spectrum. The Xe 6p levels are 0.03–0.25 eV apart.
- **The radiative decay of np goes to ns,** not to $g$, and is tracked separately.
  He 1s2s instead decays to $g$ by two-photon emission: for times comparable with
  its lifetime, add $+\Gamma^{rad}_f\rho_{ff}$ to $\dot\rho_{gg}$. The short-pulse
  fixture may neglect this decay only with its timescale stated.
- **The star-shaped Hamiltonian is an approximation.** Elimination can also produce
  Raman couplings among the excited levels and common-continuum loss terms. Neglecting
  them needs a separate check when several levels are populated (Q9).
- **The direct continuum** has an integrated pair-spectral yield in OEP 6, not a
  qualified time-local depletion rate. If a propagation model independently supplies
  a justified ground-state loss rate $\Gamma_g$, it must add
  $-\Gamma_g\rho_{gg}$ and $-\Gamma_g\rho_{fg}/2$ consistently. OEP 8 does not
  supply that missing continuum closure (Q6).

**The polarisation** in the same local, narrowband rotating-wave model:

$$
P_c^+=N\Big[\sum_fM_f^*E^-\rho_{fg}+\alpha_g^R\rho_{gg}E^++\sum_f\alpha_f^R\rho_{ff}E^+\Big]
=-N\,\partial\langle H_R\rangle/\partial E^- ,\qquad
P_{PI}^+=iN\sum_f\alpha_f^I\rho_{ff}E^+,\quad P=2\,\mathrm{Re}(P_c^++P_{PI}^+) .
\tag{8}
$$

**Why (8) follows from the amplitudes,** and is not chosen:
1. The coupling in (6) is fixed by (1), and reproduces its constant-kernel limit
   in the weak field.
2. The two-photon term of (8) equals the direct third-order dipole: both
   $\langle\psi^{(1)}|d|\psi^{(2)}\rangle$ and $\langle\psi^{(0)}|d|\psi^{(3)}\rangle$
   contribute $\tfrac12M^*E^-\rho_{fg}$ at positive frequency.
3. The conservative part, as $-\partial H_R/\partial E^-$ at fixed $\rho$,
   obeys the energy identity below. The absorptive term has its own field-work term.

Products labelled $P^+$ must remain in the positive-frequency band under these
approximations. An arbitrary broadband product need not do so; a positive-frequency
projection and a derivation of its work balance belong to Q8's deferred extension.

A spectral factor with a chosen phase, the candidate OEP 6 §4 rejected, has none of
these properties.

#### 4.4 The energy theorem and the absence of double counting

**The local-model identity.** Let $\dot\rho=-i[H_R,\rho]/\hbar+\mathcal D\rho$,
where $\mathcal D$ contains every retained loss, return and dephasing term, and let
$U=N\operatorname{Tr}(\rho H_R)$. For a pulse whose field and interaction energy vanish
at the endpoints, integration by parts and (8) give

$$
W_c=\int E\dot P_c\,dt
=\Delta U-N\int\operatorname{Tr}(H_R\mathcal D\rho)\,dt .
$$

The unitary commutator contributes zero to $\dot U$; the remaining explicit field
derivative is $-P_c\dot E$. This proves the identity for the stated local model,
including multilevel coherences and dephasing. It does not prove the accuracy of the
elimination or of a broadband extension.

For one excited level with an outgoing population sink, no return and no pure
dephasing, initially in $g$, write $\Gamma=\Gamma^{PI}+\Gamma^{rad}$. The total
work is

$$
W=N\left\{\hbar\Omega_f\rho_{ff}(t_{\rm end})
+\int\Gamma\left[\rho_{ff}(\hbar\Omega_f-\alpha_f^R|E^+|^2)
+\operatorname{Re}(\rho_{fg}H_{gf})\right]dt\right\}+W_{PI},
\qquad
W_{PI}=-2\operatorname{Re}\int P_{PI}^+\dot E^-dt .
$$

For constant $\alpha_f^I$ this last term is
$W_{PI}=N\int\Gamma^{PI}\rho_{ff}\hbar\omega_{\rm inst}(t)dt$, where
$\omega_{\rm inst}=-\operatorname{Im}(\dot E^+/E^+)$ wherever $E^+\ne0$; use the
field-work form at a zero. For a narrowband pulse this is the absorbed third-photon
energy. It is not a general broadband photon-energy prescription. With pure dephasing,
$W_c$ also contains $2N\int\gamma_f^\phi\operatorname{Re}(\rho_{fg}H_{gf})dt$.
A multilevel or returning-decay model must use the trace identity, not the one-level
sink formula (TC-1054).

The outgoing sink term includes the dressed energy and coupling energy lost from the
retained states. For radiative decay to an excited lower state, it includes both the
emitted photon and the energy retained by that lower state; it is not all fluorescence.
A plasma model must separately track the ionisation threshold, electron kinetic energy
and these sink populations.

**Reported during drafting** with Kumar et al.'s helium coefficients (10 TW cm⁻²,
203 fs), using the one-level, zero-dephasing local model:
- relative residual $3.8\times10^{-11}$;
- counting ionisation at the bare level instead gives $2.2\times10^{-4}$.

Phase 1 must regenerate this result and demonstrate convergence in integration and
quadrature tolerances. A small residual checks consistency of the chosen model, not
its physical accuracy.

**No double counting:**
- **The linear response:** (8) includes the ground-state linear polarisation. If
  propagation already includes $N\alpha_g^R E^+$ in its linear index, add only the
  correction $P_c^+-N\alpha_g^R E^+$ and the absorptive part. Excited, ionic and
  radiative-sink populations need their own response accounting; they cannot also
  remain neutral ground-state atoms in the linear model.
- **The continuum and bound terms** have disjoint final states: pairs above $I_p$
  against pairs near $\Omega_f<I_p$. Their probabilities add.
- **(2+1) ionisation** is the $\Gamma^{PI}_f$ term: a third photon, of energy $\hbar\omega$,
  per event.
- **An approximation:** (2+1) ionisation through $f$ and direct three-photon ionisation
  from $g$ reach the same continuum, so their amplitudes can interfere. Laser-induced
  continuum structure is such interference [22]. Adding probabilities neglects it. That
  requires the direct three-photon amplitude to be small against the resonant one;
  proximity to a resonance alone does not establish it under saturation or depletion.
- **A real, constant local Kerr coefficient**,
  $P^+=N(\gamma_\mathrm{bg}/2)|E^+|^2E^+$, does no net work for a vanishing pulse. Its
  reactive part must exclude the Bloch term's adiabatic tail (section 4.6).

#### 4.5 Linear response, residues and when rate equations apply

**The weak-field monochromatic particular solution** of (7) is
$\rho_{fg}=\mathcal Q_f/[2\hbar(\Omega_f-2\omega-i\Gamma^{coh}_f)]$. With (8):

$$
\gamma^\mathrm{res}(\omega)=\sum_f\frac{R_f}{\Omega_f-2\omega-i\Gamma^{coh}_f},\qquad
R_f=\frac{|M_f|^2}{\hbar}=\frac{8\varepsilon_0^2c^2\sigma_{0,f}}{\pi\hbar\omega_f^2}.
\tag{9}
$$

Here $\omega_f=\Omega_f/2$; the particular solution assumes $\rho_{gg}=1$ and
$\rho_{ff}=0$ to leading order. In atomic units this is
$\gamma^\mathrm{res}=\sum_f\delta_f/(E_f-2\hbar\omega-i\hbar\Gamma^{coh}_f)$: **the residue
is $\delta^{TP}$ per level,** the quantity OEP 6 curated.

**Three checks:**
- $\int\mathrm{Im}\gamma\,d\Omega=\pi R_f$, which gives OEP 6's
  $\mathrm{Im}\chi=32\varepsilon_0^2c^2\sigma/(\hbar\Omega^2)$ with $\sigma=\sigma_0g_L$;
- at line centre, $\mathrm{Im}\gamma=2|M|^2/(\hbar\Gamma_\mathrm{FWHM})$, Henriksson et
  al.'s Eq. 3 [8];
- (9) is Lehmberg et al.'s Eq. 3 [9] identically.

**Polarisation:** co-rotating circular light multiplies $R$ by 3/2 for $J=2$ and by 0
for $J=0$.

**Validity of (9):** a monochromatic, perturbative response with field-independent
$\Gamma^{coh}$ and negligible depletion. Near resonance,
$\Omega_R\ll\Gamma^{coh}$ is a necessary small-coherence condition; for a level
without population return, the observation time must also be short compared with
depletion. Off resonance the relevant coherence scale includes the detuning. Here $\Omega_R=|M||E^+|^2/\hbar=\sqrt{2\sigma_0/\pi}\,\Phi$ is the
two-photon Rabi frequency.

**Rate equations** replace (7) only in one of two cases:
- $\Gamma^{coh}\gg\Omega_R$ and $\Gamma^{coh}$ exceeds the rates of change of the
  envelope and of the light shift;
- or the field is chaotic, with coherence time $\tau_c\ll\tau_p$ and $\Omega_R\tau_c\ll1$.

The rate is then $W=\sigma_0g\Phi^2$, with $g$ the line (in the first case) or the
photon-statistics-weighted pair spectrum convolved with the line (in the second).

**For coherent femtosecond pulses in a dilute gas,** natural decay alone does not
justify a rate equation: $\Gamma^{rad}$ is typically $10^7$–$10^8$ s⁻¹ for the np
levels and about 51 s⁻¹ for He 1s2s. Photoionisation can dominate the coherence decay. In Kumar et al.'s case:
- $\Omega_R^\mathrm{peak}=1.04\times10^{13}$ s⁻¹;
- $\Gamma^{coh}=5\times10^{12}$ s⁻¹;
- $1/\tau=5\times10^{12}$ s⁻¹.

So (7) is required. Only the integrated weak-field yield (OEP 6 §4) is rate-like.

#### 4.6 The relation to OEP 5's $\gamma$ (OEP 6, Q8)

**What the fit already holds.** OEP 5's fitted $\gamma_\mathrm{fit}$ contains every
pole's adiabatic tail. A code adding (7) must use the background

$$
\gamma_\mathrm{bg}=\gamma_\mathrm{fit}-\sum_f\mathrm{Re}\,R_f/(\Omega_f-2\omega),
$$

which removes the co-rotating term with the same residue the Bloch model uses. The
counter-rotating term $R_f/(\Omega_f+2\omega)$ and the other permutations stay in
$\gamma_\mathrm{bg}$.

**The size of the subtraction is not settled.** For xenon at 1064 nm:
- **with fixed residues,** the three 6p poles give 4591 au, 63 % of
  $\gamma_\mathrm{fit}=7312$ au;
- **but the residue falls off resonance:** with single-pole scaling,
  $(D_\mathrm{res}/D_\omega)^2\approx0.23$, and the physical tail is about 1071 au
  (15 %).

So a fixed-residue term must be confined to a window around $\Omega_f$, with the
matching subtraction. Otherwise $\gamma_\mathrm{bg}$ is a large, model-dependent
difference (Q2, resolved: deferred to Phase 4).

#### 4.7 The dispersive part fails one test

**The test.** Eq. (9) at 248.4 nm, 760 Torr and 300 K, where $E_f-2\hbar\omega=-396.3$
cm⁻¹ for 6p[1/2]₀, against Lehmberg et al.'s $n_2$ [9]:

| strengths used | $n_2$ ($10^{-14}$ esu) |
|---|---|
| Kröll and Bischel's 6p [10], with Henriksson et al.'s 7p [8] and Drag et al.'s 6p′ [11] | −0.73 |
| Henriksson et al.'s, all levels [8] | −0.90 |
| Lehmberg et al.'s own calculation [9] | −1.54 ± 0.30 |
| measured, by self-defocusing and by compensation with CO₂ gas [9] | **−2.06 ± 0.14 and −1.82 ± 0.08** |

**What would explain it.** The measurements imply $\delta(6p[1/2]_0)\approx1600$–1820
au, or 1698–1922 au if that level alone is counted. The values for it divide in two
(TC-997, TC-1010):
- **about 1400–1900 au:** Lehmberg et al.'s own calculation, which implies 1436 au, and
  Pindzola et al.'s Dirac-Fock length form, 1541 au [12]. So Lehmberg et al.'s
  measurement and calculation agree with each other;
- **about 600–750 au:** Kröll and Bischel's measurement, 580 ± 174 au, and Henriksson et
  al.'s calculation, 752 au. Pindzola et al.'s velocity form, 432 au, is lower still.

So the tension is between strengths, not between theory and measurement.

**What does not close the gap:**
- the non-resonant background (OEP 5's $\gamma$ adds +0.07);
- photon statistics, which can only lower Kröll and Bischel's value;
- light shifts: about 0.09 cm⁻¹ at Lehmberg et al.'s intensity, against a 396 cm⁻¹
  detuning (TC-998).

No third measurement has been found. Kligler et al. (1978), requested as one, give
only a lower bound on xenon excitation at 248.4 nm and an estimate built on placeholder
dipole moments (TC-1021).

**The discrepancy remains unresolved.** The measured strengths remain the adopted
inputs (Q1), but the same residue controls absorption and dispersion in (9). If the
strength is the source of the discrepancy, population and loss predictions change too.
Thus the absorptive coefficients are qualified by their direct measurements, with this
additional cross-check failure stated; the dispersive response is unvalidated. The
0.83–1.30 comparison below does not cover this discrepancy (TC-1054).

### 5. Primary sources and data

**In hand, curated** (OEP 6's records):
- **Level positions:** the NIST level lists [13].
- **Strengths:**
  - He: Hättig et al. [3], with Drake et al. [2] and Bondy et al. [14] as checks;
  - Xe: Kröll and Bischel [10], Drag et al. [11], Henriksson et al. [8].
- **Excited-level photoionisation:**
  - Xe 6p: Kröll and Bischel [10], Pindzola [15], Payne et al. [16];
  - He 1s2s: Hummer and Storey [17].
- **The two-colour helium amplitude:** Drake et al.'s spectral distribution [2]. Kumar
  et al.'s Table II [1] is its check.
- **Light shifts and coupling coefficients of He 1s2s:** Kumar et al. [1] and
  Kaprálová-Žďánská [18]. They are held as notes only; nothing was taken. Kumar et
  al.'s continuum light-shift term is the co-rotating one only [22]; the
  counter-rotating term may be missing, of unknown size (TC-989).
- **Measured lifetimes of the Xe 6p, 6p′ and 7p levels:** Bruce et al. [23]
  (`radiative_lifetime`), with collisional quenching by Xe and Ar
  (`collisional_quenching_rate_coefficient`); TC-1003, TC-1011.
- **Measured transition probabilities:** Wiese et al. [24] (Ar), Dzierżęga et al. [25]
  (Kr) and Sabbagh and Sadeghi (Xe), as `transition_probability` (TC-1011).
- **A measured Kr 5p[1/2]₀ light shift** at 212.55 nm: Brandi et al. [26]
  (`light_shift_coefficient`; TC-1006, TC-1011).
- **Validators:**
  - Lehmberg et al. [9] (the real part);
  - Gornik et al. [6] (Xe 6p estimates and pressure shifts);
  - Shekhtman et al. [5], Dakka et al. [19] and McCann et al. [7].
- **The near-infrared NIST lines** (TC-983, TC-985). They give:
  - the radiative widths $\Gamma^{rad}_f=\sum_kA_{fk}$ of the Ne, Ar, Kr and Xe levels,
    as derived values (He 1s2s decays by two photons, at the rate of Drake et al. [2]
    and Bondy et al. [14]);
  - the single-path amplitudes of Ar, Kr and Xe.

  **Two cautions** (TC-1010):
  - **These A values are not independent.** NIST's Ar 4s-4p and Kr 5s-5p values are
    Wiese et al.'s [24] and Dzierżęga et al.'s [25].
  - **Measured lifetimes disagree for xenon.** $1/\sum A$ from NIST is about 25 % shorter
    than every measurement for 6p[3/2]₂ and 6p[1/2]₁ (Bruce et al. [23], Bowering et
    al., Ku and Setser; TC-1025 to TC-1027), and it agrees within 3 % for 6p[3/2]₁,
    6p[5/2]₃ and 6p[5/2]₂. For the 6p′ levels, Bruce et al. and Inoue et al. measure
    32–48 ns against NIST's 23–28 ns; NIST's 6p′ lines that have an A already exceed the
    measured decay rates, so its missing lines cannot close the gap (TC-1042). For 7p,
    NIST agrees with Bruce et al. within 1 %. A measured lifetime is preferred where one
    is curated.

**To be curated or requested:**
- **The excited levels' dynamic polarisabilities** (the light shifts $\mathrm{Re}\,\alpha_f$):
  - He 1s2s from Kumar et al. [1], Table I, with the caution above;
  - Kr 5p[1/2]₀, the difference from the ground state, from Brandi et al.'s measured
    shift [26], about −16 ± 10 au;
  - for Xe 6p and Ar 4p, no dynamic value has been identified. Werbowy's static Ar 4p
    polarisabilities (TC-1007) do not apply at 185–192 nm, above the levels'
    ionisation thresholds.

  The quantity is `level_polarizability` (TC-1011), kept apart from the gas's
  `polarizability`.
- **Pressure broadening and shift coefficients** of Xe 6p (the pure dephasing
  $\gamma^\phi_f$), as `pressure_line_coefficient` (TC-1040): Gornik et al.'s shifts [6]
  (TC-1039) and Hawkins et al.'s broadening, 25 ± 15 MHz/Torr (TC-1022). Hawkins et al.
  do not say whether it is a full or half width, or on which frequency scale. Collisional
  quenching is a population loss, not this dephasing; for 6p[1/2]₀ in xenon it is about
  800 times smaller than that broadening (TC-1003). Depolarisation rates put elastic
  collisions 6–100 times above quenching for the 6p levels (TC-1027).
- **Ne, Ar and Kr strengths that resolve their conventions** (TC-952). Candidates:
  - McCann et al.'s Ar measurement [7], ambiguous by about 4.5;
  - tier-3 TDSE calculations of Ar 4p by the maintainer (tRecX). Done (TC-1049): the
    one-electron model's ¹D/¹S ratio, 0.45, favours the per-level reading (0.70) over
    one component (0.14), weakly. Its absolute strengths are 2–3 times the correlated
    calculations', so argon stays unresolved;
  - the measurements on the roadmap's papers-wanted list.
- **Argon, qualified on calculated strengths** (Q7, resolved). The convention is settled
  by xenon (TC-952), and the correlated calculations agree well.
  - **The reference** is henriksson_2005's J-resolved TDDHF strength per level: J = 2 at
    24.2, 52.5 and 38.2 au; J = 0 at 34.6 and 127.8 au.
    - Its J sums (115.0, 162.4) reproduce its own LS TDHF ¹D and ¹S (113.4, 161.5).
  - **Correlation raises both terms by the same factor.** hattig_1998b's CCSD is 1.22 times
    TDHF for ¹S (197.7 au, the same in either convention) and 1.23 times for ¹D: its
    printed single-component $S^{0f}$, 104.59 au, gives 139.5 au per level (×4/3).
    - CC2 and MCSCF (Sundholm et al. [4], as Hättig et al. reprint it) lie within 4 % of
      CCSD.
  - **The stated uncertainty** is therefore −5 % and +30 % about the TDDHF values,
    labelled "calculated only".
  - **McCann et al.'s measurement [7]** agrees in its 0.85 reading.
  - **The one-electron TDSE** is 2–3 times high and does not enter (TC-1049).
  - **Inputs** (TC-1068):
    - **Radiative widths:** Wiese et al.'s complete 4s–4p A values [24]. Sadeghi et
      al.'s lifetimes of the two-photon-excited levels (24–29 ns, at 184–189 nm)
      confirm them independently.
    - **Quenching by argon:** 0.13–0.47 × 10⁻¹⁰ cm³ s⁻¹ (TC-1074). Sadeghi et al.'s
      table values for 4p[1/2]₀ and 4p[5/2]₂ are lower limits by about 2, reconciled with
      Chang and Setser's and Nguyen and Sadeghi's state-to-state data; 4p′[1/2]₀ has
      only Sadeghi et al.'s lower limit.
    - **Photoionisation of each level at its drive:** 2.0–2.6 Mb (Petrov et al.,
      calculated). It agrees with the J-averaged values at 248 nm, and for xenon 6p it is
      0.75 of Kröll and Bischel's measurement. Baig et al.'s measured threshold values
      for three 4p′ levels are 2.8–3.1 times the same calculation. So its absolute scale
      for argon is uncertain by up to a factor of 3 (TC-1074).
    - **Not curated, and stated as unqualified:** the light shifts and the pressure
      broadening.
- **The real-part tension:** the earlier Xe $n_2$ measurements and theory that
  Lehmberg et al. cite. They are not yet requested.
- **A channel-resolved, nondegenerate continuum amplitude $M_c(\Omega,\nu)$** at small
  splits, for the continuum's time dependence (out of scope, Q6). Candidates:
  - the RPAE and R-matrix groups' methods;
  - tier-3 TDSE calculations for He and Ar: done (TC-1048, TC-1049; section 4.2).

**Uncertainty:** each published coefficient carries the level's assessment from TC-952
and its spread across sources:
- He 1s2s: 1–4 % source spread;
- Xe: 0.83–1.30 of the measured values.

The draft's light-shift sensitivity estimate (±15 % for He 1s2s) is stated separately
from source uncertainties. In its helium example this moves the transferred population
by about 10 %; Phase 1 must reproduce the sensitivity, including the pulse parameters.

### 6. Proposed interface

**Qualification describes coefficients and a regime, not all pulse dynamics.**
The Q1 and Q7 decisions stand. Q9 initially qualifies He 1s2s as the nonlinear
dynamical example and retains Xe and Ar coefficients for their perturbative response;
a coupled strong-field Xe or Ar model is not qualified by those coefficients alone.
For Xe and Ar, missing light shifts preclude an unconditional strong-field population
prediction; for Ar, missing photoionisation also
precludes a qualified (2+1) yield. Neglecting these effects requires independent evidence
that they are small for the specified pulse. Missing pressure widths preclude a
pressure-dependent line shape. He requires Phase 1's curated shifts for the strong-field
fixture. No missing coefficient is silently zero.

**A per-level response record** extends `optmat.gas.two_photon_bound`. Read-only, frozen,
identified (OEP 1):

```text
from optmat.gas.two_photon_bound import two_photon_bound_response

resp = two_photon_bound_response("He")        # TwoPhotonBoundResponse; Identified
resp.qualified                                 # ("1s.2s|1S|0",)
line = resp.line("1s.2s|1S|0")
line.omega                                     # Omega_f, rad/s (NIST)
line.j, line.reach                             # 0, "linear"
line.residue(polarization="linear")            # R_f, C m^4 V^-3 s^-1 (delta in au also)
line.residue_uncertainty                       # relative, from the assessment
line.d1, line.d_eff                            # rad/s
line.sigma_max(0.01)                           # empirical leading-order Gaussian estimate, rad/s
line.photoionisation(omega)                    # sigma^PI_f, m^2; NaN where not curated
line.radiative_width                           # Gamma^rad_f, s^-1: a measured lifetime where
                                               #   curated, else NIST A values (TC-983);
                                               #   He 1s2s its two-photon decay rate
line.kernel(omega_a, omega_b)                  # M_f(omega_a, omega_b) / M_f(Omega_f/2, Omega_f/2);
                                               #   He 1s2s only, NaN elsewhere
resp.gamma_resonant(omega, coherence_rate)     # Eq. (9), co-rotating poles; complex;
                                               #   coefficient qualification and real-part validation
                                               #   reported separately, by species and level
```

The record also carries machine-readable provenance, per-coefficient availability,
valid frequency intervals, polarisation restrictions, uncertainty type and confidence
level (or "not stated"), and the empirical nature of `d_eff` and `sigma_max`. A source
spread or argon's −5 %/+30 % model envelope is not a statistical confidence interval.
The derivation and coefficient-source digests enter its OEP 1 identity.

**Refusals:**
- Ne and Kr raise `UnavailableModelError`. The curation record
  (`two_photon_bound_levels`) stays available.
- A level that is listed but not qualified (a higher series member, $J=1$) is absent
  from `qualified`, and `line()` names the reason.
- A valid request evaluated outside available spectral support returns NaN, with an
  availability mask, as OEP 6 does. Invalid arguments (including nonfinite or negative
  coherence rates and an invalid error tolerance) raise `ArgumentError`.
- `gamma_resonant` takes a declared set of levels and a scalar or explicitly
  level-indexed coherence rate in s⁻¹; it states its polarisation and near-resonance
  support. At a zero-width pole it is unavailable. It is not a qualified extrapolation
  to the infrared or a background-subtraction recipe.

**Not supplied,** because they depend on the field or are not curated:
- the light shifts of the excited levels, until curated;
- the photon statistics;
- the band split;
- the mixed-pair policy;
- the integrator.

**Composition:**
- **With OEP 6's continuum model:** the two terms are disjoint and add (section 4.4).
- **With OEP 5's $\gamma$:** `gamma_resonant` gives the co-rotating term a code would
  subtract. No subtracted $\gamma_\mathrm{bg}$ is offered before the deferred Phase 4
  matching decision (Q2).

### 7. Implementation strategies

**Where the work belongs:** `optmat.gas`, next to the curation record that OEP 6 built.

**What optmat supplies.** The options:

| option | optmat supplies | for | against |
|---|---|---|---|
| A. Coefficients only (recommended) | per-level residues, widths, error estimates, $\sigma^{PI}$, kernel | the facts are optmat's; the field, statistics and solver are the code's (the OEP 5 and OEP 6 pattern) | each code must implement (7)–(8) |
| B. A + a public single-atom integrator | as A, plus `evolve(field)` | one reference implementation | optmat would own a field model and a solver, against its scope |
| C. A linear complex $\chi_2$ only | Eq. (9) | simple | wrong where it matters: saturation and $\Gamma^{PI}\propto I$ (section 4.5) |

**Recommendation: A.** Add a reference integrator in `tools/` (not public). It
generates the parity fixtures: the weak-field identity, Kumar et al.'s case and the
single-atom energy identity. These are checks of the material coefficients and their
normalisation. The fixtures are stored with digests so that solprop can test against
them (Q4, resolved). The propagation polarisation, plasma and energy-accounting
equations and their implementation remain solprop's responsibility (TC-1057).

**Data and vocabulary (Q10, approved; TC-1055):**
- Reuse the existing strengths and photoionisation quantities, and measured lifetimes
  (`radiative_lifetime`) where curated, otherwise NIST's sum of A (TC-983).
- Add `two_photon_decay_rate` in the existing `absorption` family: a spontaneous
  two-photon transition rate, in `1/s`, with a `label` abscissa and required
  `upper_level`, `lower_level`, `isotope` and `method` parameters. Store the printed
  rate; a lifetime is derived, not substituted into the extracted data.
- Add `two_photon_decay_spectrum` in that family: the printed rate density per unit
  $y=\hbar\omega_1/(E_u-E_l)$, in `1/s`, with an `energy_fraction` abscissa whose
  unit is `1`. Required parameters are `upper_level`, `lower_level`, `isotope`,
  `method`, `counting` (`photons` or `pairs`) and `definition` (the source equation
  and integration convention). Infinite nuclear mass is explicit in the method,
  not silently identified with a finite-mass isotope. Validate $0\le y\le1$.
- Preserve the printed half-range and raw cells. For Drake's photon-counted symmetric
  distribution, the pair-normalised full-range density is $A(y)/2$, so
  $\Gamma=\tfrac12\int_0^1A(y)dy=\int_0^{1/2}A(y)dy$. Reflection, unit conversion
  and the derived kernel belong in the documented loader/derivation. Do not store
  this squared-amplitude spectrum as `two_photon_transition_amplitude`.
- Phase 1's helium light shifts and Phase 2's further shifts, quenching and dephasing
  use `level_polarizability`,
  `light_shift_coefficient`, `collisional_quenching_rate_coefficient` (TC-1011) and
  `pressure_line_coefficient` (TC-1040).
- **The kernel's test** uses `two_photon_transition_amplitude` (TC-1040): hydrogen's
  exact signed amplitudes against the split.

**The derivation:** stored per gas as JSON, with a freshness test, like the OEP 5 and
OEP 6 records.

**The kernel of He 1s2s:** after extraction into a declared quantity, a spline of
Drake et al.'s tabulated $A(y)$, converted by
$|s(\nu)/s(0)|^2=\frac{A(y)}{A(1/2)}\frac{(1/4)^3}{[y(1-y)]^3}$. It is checked against
Kumar et al.'s two-colour value. The table gives a squared modulus, not a complex
amplitude. Publish only the on-shell ratio ($\omega_a+\omega_b=\Omega_f$) in the
sub-resonant interval where a real, continuous, nonzero amplitude permits the positive
root after fixing $M(\Omega_f/2)>0$. Do not infer signs across a zero, phases across a
pole or off-shell values. The endpoints require a demonstrated limit, not division by
$y^3(1-y)^3$. This kernel qualifies perturbative yields; a time-dependent two-colour
Hamiltonian and its polarisation are deferred by Q8. In particular,
Kumar's strong two-colour scheme also drives 1s2s–1s4s and is not a two-level test.

Drake's spectrum and the helium decay rates currently occur only in source-record
notes. Phase 1 must add extraction scripts, datasets and register entries before they
become shipped coefficients. Q10 now specifies their quantities and normalisation.
Use Drake Table II for the normalised shape and Bondy Table VI's finite-mass,
relativistically corrected helium rate, 50.93 ± 0.07 s⁻¹, for decay. Preserve the
estimated theoretical uncertainty without assigning it a statistical confidence level.
This combines sources with different approximations, explicitly; the relative kernel
shape is not promoted to Bondy's rate accuracy. Drake's half-spectrum quadrature gives
51.1457 s⁻¹ against its printed 51.3 s⁻¹; keep this check, without silently rescaling
raw data. Bondy's extracted Table VI is scrambled and must be repaired and checked
against its table image before ingestion (TC-1055).

**Cost:** milliseconds per record. No new dependencies.

### 8. Verification and validation

These are the drafting calculations reported as Phase 1 targets. Their original
residuals are not acceptance tolerances: integration, quadrature, interpolation and
finite-field errors must be separated by convergence studies. A comparison using a
paper's own input coefficients reproduces that model; it is not an independent test of
those coefficients.

| check | quantity | result for this draft | tolerance set by |
|---|---|---|---|
| per-level convention | $\sigma_0$ from (1) against TC-952's constant | $10^{-11}$ | rounding |
| weak-field limit | (7)–(8) yield against OEP 6's pair-spectral $P$ (He, 0.01 TW cm⁻², 50 fs) | 0.99973 | convergence in step and intensity separately; Stark shifts and ionisation also vanish only in the weak-field limit |
| Kumar et al.'s one-colour case | $P_e$ after the pulse (10 TW cm⁻², 9 meV, $\Delta=+14.3$ meV) | 0.276; text "about 28 %" | the paper's printed precision |
| local energy identity | relative residual, one level without pure dephasing | $3.8\times10^{-11}$ | time integration and work quadrature convergence |
| coupling | Kumar et al. against Kaprálová-Žďánská [18] | 1.9 % | their spread |
| two-colour kernel | Drake et al. against Kumar et al.'s Table II | 1.5–2 % | their spread |
| split error | quadrature against (5), He, 10 fs and 1 fs | $1.353\times10^{-4}$ both; 1.38 % against 1.35 % | the $\nu^4$ term |
| field-strength error | Kaufman et al.'s three-level model against its eliminated form [21] | $P_e$ 0.843 against 0.869 (3.1 %; 2.5 % printed) | the admixture $2(dE_0/2\hbar D)^2$ |
| Xe strengths | Henriksson et al. against Kröll and Bischel's and Drag et al.'s measurements | 0.83–1.30 | the measurements' errors |
| Xe 6p single path | against Gornik et al.'s estimates | $1.20\times10^4$ against $1.2\times10^4$ s⁻¹ | their two digits |
| real part | (9) against Lehmberg et al.'s $n_2$ | fails by 2.2–2.8 | recorded as a failure (Q1) |

**Review checks (2026-10-12; TC-1054):** the multilevel terms were compared directly
with the matrix commutator; Kumar's Table I and Lehmberg's Table 1 were checked against
the primary table images. Independent integration of Kumar's printed Table I
coefficients gives $P_e=0.27612491$ for a 202.77234 fs Gaussian with 9 meV spectral
FWHM and $\hbar(2\omega-\Omega_f)=+14.3$ meV; TC-1054 records the solver and
convergence settings. Phase 1 must additionally test Hermiticity, positivity,
trace plus sink populations, two excited levels, nonzero pure dephasing, radiative
return, invalid inputs and spectral masks. Preserve the paper's coefficients for its
reproduction fixture, and separately report results with the preferred curated inputs.

**Independent validators:**
- tier-3 TDSE runs of He 1s2s dynamics at 1–10 fs and with chirp (TC-1048, TC-1051).
  - **At 10¹² W cm⁻², at haCC's propagated resonance:** the 1s2s populations equal the
    weak-field form to ≤ 1 % for every duration and chirp (3 % at 1 fs).
  - **The fitted coupling:** |M|² = 12.31 au², 0.04 % from Hättig et al.'s δ [3].
  - **Beyond lowest order (10¹³, 10¹⁴ W cm⁻²):** the departures, down to 0.08 of lowest
    order at 10 fs, 10¹⁴, are Phase 1's Bloch-model targets;
- solprop parity on the stored fixtures.

### 9. Documentation

**The manual page** {doc}`../two_photon` gains a section *Bound two-photon resonances
in propagation*. It covers:
- the material conventions established by equations (6)–(9), with single-atom checks
  distinguished from the propagation equations owned by solprop;
- the material data optmat supplies and the response solprop constructs from them;
- material energy conventions and checks, with propagation energy bookkeeping in solprop;
- the validity limits;
- the failed real-part test;
- the refusals;
- the RDW diagnostic of section 2.1, separating continuum loss from bound excitation,
  and approved scope decisions from deferred scientific qualification.

**Also:**
- {doc}`../propagation` gains a pointer;
- docstrings carry units;
- a gallery example: He 1s2s transfer against duration for a 120 nm pulse, with the
  weak-field identity overlaid.

### 10. Compatibility and impact

- **Additive only.** No recommended value or default changes.
- **OEP 5's $\gamma$ is unchanged.** The subtraction is the code's, and is not applied
  in optmat.
- **OEP 6:** its Phase 2b and optional Phase 3 have passed to this proposal. OEP 6 is
  recorded as Final with Phases 1 and 2a implemented (TC-1053).
- **solprop:** new parity fixtures.

### 11. Alternatives considered

- **A spectral factor $H(\Omega)$ multiplying the pair spectrum, with a chosen phase.**
  Rejected by OEP 6's review, and here because it has no energy theorem and is not the
  amplitude's consequence (section 4.3).
- **Linear complex $\chi_2$ (option C).** It is correct only below saturation and with
  field-independent widths. Kumar et al.'s case violates both.
- **Rate equations with a Lorentzian line.** They fail for the coherent femtosecond
  example (section 4.5). They remain valid in the fast-dephasing or suitably chaotic
  limits stated there.
- **Using all five gases with labelled best-effort strengths.** Rejected for Kr, whose
  measurements and calculations disagree by up to 94 times. Accepted for Ar (Q7), whose
  convention is settled and whose calculations agree to about 25 %. Ne could follow
  argon's route, from calculations alone, but is not proposed here.

### 12. Phases

1. **Coefficients, helium dynamics and the weak-field handoff.** Deliverables:
   - `two_photon_bound_response` and its records for He, Xe and Ar, with separate
     coefficient, perturbative and nonlinear-dynamical qualifications (Q8–Q9);
   - the Q10 quantities and energy-fraction abscissa, documented in the vocabulary
     and data manual before ingestion; extraction of Drake's spectrum and Bondy's
     corrected helium decay rate, with the Table VI extraction repaired;
   - the radiative widths: measured lifetimes where curated (Bruce et al. for Xe),
     otherwise NIST's sum of A (TC-983), and Bondy's He 1s2s rate;
   - the He 1s2s on-shell kernel, restricted to perturbative yields;
   - He 1s2s's light shifts from Kumar et al.'s Table I, including the separate bound
     and continuum contributions and TC-989's limitation, moved here from Phase 2;
   - argon's 4p strengths: hattig_1998b's ¹De rows taken as data (`one_component`,
     ×4/3), and argon's assessment changed to "calculated", by the author's decision
     (TC-1053);
   - `gamma_resonant` for the supported perturbative response; no qualified coupled
     strong-field Xe or Ar model and no arbitrary broadband kernel;
   - the internal single-atom reference integrator in `tools/`, the local
     energy-identity tests of material conventions and stored parity fixtures;
     both the exact Kumar input reproduction and the preferred
     curated-input comparison, with their distinct purposes;
   - the tests of section 8, including decay normalisation, kernel symmetry, domain
     masks and convergence; the manual section and gallery example;
   - a documented handoff for section 2.1's continuum-first RDW investigation, using
     OEP 6's existing coefficients and integrated constraints. This is guidance for
     solprop; its continuum diagnostic does not wait for completion of this phase.
     A continuum-response SEP owns the polarisation, plasma coupling and propagation
     energy accounting; a coupled optmat OEP would cover any additional material
     inputs it requires, rather than adding a continuum phase here (TC-1057);
   - register entries and a changelog line. Extraction and public-surface changes
     require the library and release lanes in addition to the full gate.
2. **Further light shifts and collision data.** Deliverables:
   - Xe 6p's pressure coefficients from Gornik et al.;
   - additional excited-level shifts, quenching and dephasing where primary sources
     qualify them; the light-shift sensitivity stated per level.
   These additions do not by themselves qualify strong-field multilevel dynamics.
3. **More levels and improved coefficients.** Ne 3p and Kr 5p as their strengths are
   resolved, each with its own assessment. Argon's measured strengths, photoionisation,
   light shifts and pressure broadening as curated. Strong-field multilevel
   qualification requires the deferred Q9 evidence and an explicit scope revision.
4. **(Optional) The background Kerr term.** A windowed subtraction from OEP 5's
   $\gamma$, if Q1 and Q2 are settled.

General frequency-dependent/two-colour dynamics (Q8), strong-field multilevel
qualification (Q9), and the continuum's time dependence (Q6) are not silently added to
these phases. They remain on the roadmap. Approval of Q8–Q10 resolves the initial
scope and data choices, not those scientific extensions or acceptance of the whole OEP.

### 13. Open questions

Q1 to Q6 resolved by the maintainer on 2026-10-09, as recommended (TC-984); Q7 on
2026-10-12 (TC-1053).

- **Q1 (the real part; resolved).** `gamma_resonant` publishes the dispersive part,
  labelled as failing Lehmberg et al.'s xenon $n_2$ by 2.2–2.8. The qualification of the
  absorptive and dispersive parts is stated separately. The earlier xenon $n_2$
  literature is to be requested.
- **Q2 (the background; resolved).** Deferred to Phase 4. Until then optmat offers no
  subtracted $\gamma$; with fixed residues the tail is 63 % of $\gamma_\mathrm{fit}$ at
  1064 nm, with falling residues 15 %.
- **Q3 (the split rule; resolved).** $D_\mathrm{eff}=0.75\,D_1$ is adopted, labelled as a
  bound from the curated cases. The review clarifies that this is an empirical
  leading-order estimate, not a mathematical guarantee for arbitrary spectra or
  cancelling amplitudes (TC-1054). The Kr and Xe single-path cases rest on an assumed
  sign of the second path, and the label says so.
- **Q4 (the integrator; resolved).** In `tools/` only, with stored fixtures for parity;
  no public integrator.
- **Q5 (photon statistics; resolved).** The residues are for coherent light. Chaotic or
  partially coherent light is the code's responsibility; optmat states this and supplies
  nothing more.
- **Q6 (the continuum's time dependence; resolved).** Not a phase of this OEP.
  The response belongs in a solprop SEP; a coupled optmat OEP would cover any new
  material quantities, models or interfaces it requires. Published small-split
  amplitudes would provide additional material information, but are not a prerequisite
  to proposing an explicitly approximate response. Such a response must be derived and
  validated; OEP 6's integrated constraints alone do not qualify it (TC-1057).

- **Q7 (argon; resolved 2026-10-12).** Should argon be qualified without a decisive
  measurement? **Decided by the author, as recommended:**
  - argon's 4p J = 0 and J = 2 levels are qualified on henriksson_2005's J-resolved
    TDDHF strengths;
  - the uncertainty is −5 % and +30 %, from the correlated calculations, labelled
    "calculated only";
  - for absorption and populations only.

  See section 5 and TC-1053.

Q8–Q10 were resolved by the maintainer on 2026-10-12, approving the AI review's
recommendations (TC-1055). The experimental priority is recorded in TC-1056 and
section 2.1. The wider scientific work remains deferred, not declared validated.

- **Q8 (broadband dynamics; resolved).** The first dynamical response is the local
  effective model, within its spectral and field-strength limits. Publish helium's
  two-colour kernel for perturbative on-shell yields only. A general frequency-dependent
  Hamiltonian, polarisation and work balance are deferred to a later extension. This
  scope still permits short coherent pulses; width relative to an intermediate-state
  detuning is the relevant elimination condition, not width relative to the final
  state's natural linewidth. An intense overlapping pump needs a separate dressing
  check; a closed pump–VUV two-photon channel does not establish negligible dressing.
- **Q9 (several excited levels; resolved).** Qualify He 1s2s as the first nonlinear
  dynamical example, with its light shifts in Phase 1. Retain the approved Xe and Ar
  coefficients and their weak-field multilevel applications. Defer qualified coupled
  strong-field populations and ion yields until effective excited-state couplings,
  shared continuum loss and correlated dephasing are supplied or shown negligible.
  A maintainer decision about scope cannot establish that scientific condition.
- **Q10 (helium source data; resolved).** Add the two decay quantities and the
  dimensionless energy-fraction abscissa specified in section 7. Use Drake's Table II
  for the spectrum and Bondy's Table VI for the preferred helium decay rate. Preserve
  printed values, counting conventions and uncertainty meanings. Repair the scrambled
  Bondy table extraction before ingestion; derive lifetimes and kernel ratios rather
  than treating them as printed data. This is Phase 1 curation, not new atomic physics.

### 14. References

Curated as source records unless marked (the record's id in brackets).

1. S. Kumar, J. D. Piel, C. H. Greene and N. Shivaram, "Broadband femtosecond lasers
   enable efficient two-photon excitation of the ultranarrow linewidth singlet 1s2s
   state in helium," *Phys. Rev. Research* **8**, 013009 (2026). doi:10.1103/nfbf-h6hc
   [kumar_2026]
2. G. W. F. Drake, G. A. Victor and A. Dalgarno, "Two-photon decay of the singlet and
   triplet metastable states of helium-like ions," *Phys. Rev.* **180**, 25 (1969).
   doi:10.1103/PhysRev.180.25 [drake_1969]
3. C. Hättig, O. Christiansen and P. Jørgensen, "Coupled cluster response calculations
   of two-photon transition probability rate constants for helium, neon and argon,"
   *J. Chem. Phys.* **108**, 8355 (1998). doi:10.1063/1.476262 [hattig_1998b]
4. D. Sundholm, A. Rizzo and P. Jørgensen, "Multiconfiguration self-consistent field
   quadratic response calculations of the two-photon transition probability rate
   constants for argon," *J. Chem. Phys.* **101**, 4931 (1994). doi:10.1063/1.467415
   [sundholm_1994]
5. D. Shekhtman, M. A. Mustafa and N. J. Parziale, "Two-photon cross-section
   calculations for krypton in the 190–220 nm range," *Appl. Opt.* **59**, 10826
   (2020). doi:10.1364/AO.410806 [shekhtman_2020]
6. W. Gornik, S. Kindt, E. Matthias and D. Schmidt, "Two-photon excitation of xenon
   atoms and dimers in the energy region of the 5p⁵6p configuration," *J. Chem. Phys.*
   **75**, 68 (1981). doi:10.1063/1.441856 [gornik_1981]
7. M. P. McCann, C. H. Chen and M. G. Payne, "Two-photon (vacuum ultraviolet + visible)
   spectroscopy of argon, krypton, xenon, and molecular hydrogen," *J. Chem. Phys.*
   **89**, 5429 (1988). doi:10.1063/1.455595 [mccann_1988]
8. J. Henriksson, P. Norman and H. J. Aa. Jensen, "Two-photon absorption in the
   relativistic four-component Hartree-Fock approximation," *J. Chem. Phys.* **122**,
   114106 (2005). doi:10.1063/1.1869469 [henriksson_2005]
9. R. H. Lehmberg, C. J. Pawley, A. V. Deniz, M. Klapisch and Y. Leng, "Two-photon
   resonantly-enhanced negative nonlinear refractive index in xenon at 248 nm," *Opt.
   Commun.* **121**, 78 (1995). doi:10.1016/0030-4018(95)00433-9 [lehmberg_1995]
10. S. Kröll and W. K. Bischel, "Two-photon absorption and photoionization
    cross-section measurements in the 5p⁵6p configuration of xenon," *Phys. Rev. A*
    **41**, 1340 (1990). doi:10.1103/PhysRevA.41.1340 [kroll_1990]
11. C. Drag, F. Marmuse and C. Blondel, "Measurement of the two-photon excitation
    cross-section of the 6p′[3/2]₂ and 6p′[1/2]₀ levels of Xe I at the wavelengths
    224.3 and 222.6 nm," *Plasma Sources Sci. Technol.* **30**, 075026 (2021).
    doi:10.1088/1361-6595/abfbeb [drag_2021]
12. M. S. Pindzola, M. G. Payne and W. R. Garrett, "Two-photon excitation of xenon,"
    *Phys. Rev. A* **24**, 3115 (1981). doi:10.1103/PhysRevA.24.3115 [pindzola_1981b]
13. A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team, NIST Atomic Spectra
    Database (version 5.12), https://physics.nist.gov/asd. doi:10.18434/T4W30F
    [nist_asd_2026]
14. A. T. Bondy, D. C. Morton and G. W. F. Drake, "Two-photon decay rates in heliumlike
    ions: Finite-nuclear-mass effects," *Phys. Rev. A* **102**, 052807 (2020).
    doi:10.1103/PhysRevA.102.052807 [bondy_2020]
15. M. S. Pindzola, "Photoionization of excited-state xenon," *Phys. Rev. A* **23**, 201
    (1981). doi:10.1103/PhysRevA.23.201 [pindzola_1981]
16. M. G. Payne, C. H. Chen, G. S. Hurst, S. D. Kramer, W. R. Garrett and M. Pindzola,
    "Resonance ionization spectroscopy schemes for Ar, Kr and Xe," *Chem. Phys. Lett.*
    **79**, 142 (1981). doi:10.1016/0009-2614(81)85306-7 [payne_1981]
17. D. G. Hummer and P. J. Storey, "Recombination of helium-like ions — I.
    Photoionization cross-sections and total recombination and cooling coefficients
    for atomic helium," *Mon. Not. R. Astron. Soc.* **297**, 1073 (1998).
    doi:10.1046/j.1365-8711.1998.2970041073.x [hummer_1998]
18. P. R. Kaprálová-Žďánská, "Coherent superposition mediating energy flow between
    fundamental and third-harmonic modes," *Phys. Rev. A* **111**, 043105 (2025).
    doi:10.1103/PhysRevA.111.043105 [kapralova_zdanska_2025]
19. M. A. Dakka et al., "Laser-based metastable krypton generation," *Phys. Rev. Lett.*
    **121**, 093201 (2018). doi:10.1103/PhysRevLett.121.093201 [dakka_2018]
20. N. Dudovich, B. Dayan, S. M. Gallagher Faeder and Y. Silberberg, "Transform-limited
    pulses are not optimal for resonant multiphoton transitions," *Phys. Rev. Lett.*
    **86**, 47 (2001). doi:10.1103/PhysRevLett.86.47 [dudovich_2001]
21. B. Kaufman, T. Rozgonyi, P. Marquetand and T. Weinacht, "Adiabatic elimination in
    strong-field light-matter coupling," *Phys. Rev. A* **102**, 063117 (2020).
    doi:10.1103/PhysRevA.102.063117 [kaufman_2020]
22. P. L. Knight, M. A. Lauder and B. J. Dalton, "Laser-induced continuum structure,"
    *Phys. Rep.* **190**, 1 (1990). doi:10.1016/0370-1573(90)90089-K [knight_1990]
23. M. R. Bruce, W. B. Layne, C. A. Whitehead and J. W. Keto, "Radiative lifetimes and
    collisional deactivation of two-photon excited xenon in argon and xenon," *J. Chem.
    Phys.* **92**, 2917 (1990). doi:10.1063/1.457939 [bruce_1990]
24. W. L. Wiese, J. W. Brault, K. Danzmann, V. Helbig and M. Kock, "Unified set of atomic
    transition probabilities for neutral argon," *Phys. Rev. A* **39**, 2461 (1989).
    doi:10.1103/PhysRevA.39.2461 [wiese_1989]
25. K. Dzierżęga, U. Volz, G. Nave and U. Griesmann, "Accurate transition rates for the
    5p–5s transitions in Kr I," *Phys. Rev. A* **62**, 022505 (2000).
    doi:10.1103/PhysRevA.62.022505 [dzierzega_2000]
26. F. Brandi, W. Hogervorst and W. Ubachs, "High-resolution vacuum-ultraviolet and
    ultraviolet photoionization spectroscopy of krypton," *J. Phys. B* **35**, 1071
    (2002). doi:10.1088/0953-4075/35/4/329 [brandi_2002]

## Part II — Implementation record

None yet.

## Part III — History

- 2026-10-09: Drafted, from the derivation that answers OEP 6's Q9 for the bound
  levels.
- 2026-10-09: Q1 to Q6 resolved by the maintainer, as recommended (TC-984).
- 2026-10-09: Revised in Draft with the literature of TC-985 to TC-1010: the
  field-strength condition, the pole prescription, the xenon strengths divided in two,
  the measured lifetimes and the light-shift sources.
- 2026-10-09: The radiative, light-shift and quenching quantities approved (TC-1011);
  sections 5, 7 and 12 updated to use them.
- 2026-10-09: Revised with the literature of TC-1021 to TC-1041: hydrogen's exact
  kernel and the sign of $\kappa$, no third xenon strength, the xenon lifetimes against
  NIST, the pressure-line and amplitude quantities.
- 2026-10-09: The xenon 6p′ lifetimes settled by Inoue et al. (TC-1042).
- 2026-10-11: The two-photon tRecX runs added as tier-3 evidence (TC-1048, TC-1049).
- 2026-10-11: Run 5's helium 1s2s discrepancy resolved; the weak-field form and the 1s2s
  strength confirmed (TC-1051).
- 2026-10-12: Argon's 4p levels qualified on calculated strengths, by the author (Q7, TC-1053).

- 2026-10-12: John Travers confirmed his human review. GPT-6 completed the AI review
  and revised the coherence sign, local energy identity, bandwidth-error qualification,
  coefficient availability and validation contract (TC-1054). Review status is Reviewed;
  lifecycle status remains Draft. Q8–Q10 record the remaining design gaps. The review
  did not rerun the TDSE calculations or reproduce every drafting numerical claim.

- 2026-10-12: The maintainer approved the recommended resolutions of Q8–Q10
  (TC-1055). The local dynamical scope and perturbative kernel restriction are now
  explicit; helium shifts and the approved decay-data schema are in Phase 1. Broader
  dynamics remain deferred. The initial Ar/He RDW experiment and a continuum-first
  diagnostic were added to section 2.1 (TC-1056). Status remains Draft; Reviewed.
- 2026-10-12: The maintainer clarified the package boundary (TC-1057): optmat supplies
  material properties and basic data; solprop owns the propagation polarisation,
  electron accounting and energy balance. A continuum-response SEP may require a
  coupled material-data OEP. Section 2.1 distinguishes pump detuning from the ultimate
  observed cut-off and records the effective-continuum suggestion as an unqualified
  candidate. No new phase, OEP or SEP is created; status remains Draft; Reviewed.
- 2026-10-12: Argon's radiative widths, quenching and photoionisation at the drive curated (TC-1068).
- 2026-10-12: Argon's quenching reconciled and its photoionisation scale questioned (TC-1074).

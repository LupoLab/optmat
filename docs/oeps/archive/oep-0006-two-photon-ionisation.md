# OEP 6 — Two-photon ionisation and absorption of the noble gases

| | |
|---|---|
| **OEP** | 6 |
| **Title** | Two-photon ionisation and absorption of the noble gases |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 (2026-10-08); GPT-6 (review revisions, 2026-10-08) |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | John Travers (human); GPT-6 (AI) |
| **Type** | Feature |
| **Created** | 2026-10-08 |
| **Resolution** | Accepted 2026-10-08 by John Travers in the review discussion; Q9 remains a prerequisite for implementing the complex evaluator in Phase 2. Final 2026-10-12 by the author: Phases 1 and 2a implemented; Phase 2b and the optional Phase 3 transferred to OEP 8. |
| **Requires** | [OEP 1](oep-0001-propagation-interface.md): model identity, digests and the declared public surface; [OEP 5](oep-0005-hyperpolarizability-support.md): the resonance limits $E_1$, $E_2$ and the band edge |
| **Supersedes** | — |
| **Superseded-by** | [OEP 8](../oep-0008-bound-two-photon-response.md), for Phase 2b and Phase 3 only |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phase 1 merged from `feature/oep-0006-continuum`; Phase 2a (curation) merged from `feature/oep-0006-bound-levels`; Phase 2b and Phase 3 transferred to OEP 8 |
| **Manual pages** | {doc}`../../two_photon` (*The model*; *Bound two-photon excitation and the excited levels' ionisation*); {doc}`../../propagation`; {doc}`../../ionisation_rates` |

## Part I — Proposal

### 1. Abstract

Intense vacuum- and deep-ultraviolet pulses in gas-filled fibres and capillaries lose
energy to two-photon ionisation once the photon energy exceeds half the ionisation
energy, $\hbar\omega > I_p/2$. Just below that threshold, they lose it to two-photon
excitation of even-parity bound levels, whose population a further photon may ionise.
optmat holds the gases' strong-field ionisation rates, but no two-photon cross section a
propagation code can use. Evaluating those rates with the pump's carrier does not obtain a
separately calibrated two-photon response of its vacuum-ultraviolet component.

This OEP proposes a model of each noble gas's generalised two-photon ionisation cross
section $\sigma^{(2)}(\omega)$, built on the curated literature, in two steps:
- **Phase 1:** the continuum, from the two-photon threshold to the first one-photon
  resonance. Each gas's model has thresholds, labelled support, uncertainty and an
  identity.
- **Phase 2:** the bound two-photon resonances below the threshold, as the imaginary part
  of a complex two-photon response, subject to the derivation required in Q9.

The model states facts; what a propagation code does with them is its own decision, as
for OEP 5 (TC-868). It succeeds if a paired solprop proposal can build a two-photon loss
channel from it without copying numbers, and if every value is traceable to a curated
record and its stated uncertainty. It also succeeds if it says plainly where the data do
not support a model: krypton, and xenon near threshold.

### 2. Motivation and scope

**The need.** Consider a VUV pulse in a noble-gas-filled fibre or capillary at a peak
intensity of 10¹³ W cm⁻². This is an illustrative modelling case, not a measured
wavelength limit. Argon's two-photon cross section, about 1.2×10⁻⁵⁰ cm⁴ s over the
calculated off-resonant range 8.0–11.4 eV (McKenna and van der Hart [1]; TC-856), gives an
energy-loss coefficient (section 4) of
$\beta = 2N\sigma^{(2)}/\hbar\omega \approx 5\times10^{-15}$ m W⁻¹ at 1 bar and
293.15 K near 155 nm (ideal-gas density). Thus $\beta I \approx 500$ m⁻¹ at the assumed
intensity, an initial attenuation length of millimetres.

Whether two-photon ionisation limits a generated pulse's energy or wavelength is a
question a propagation code can answer only with curated cross sections and a qualified
loss model. The existing strong-field
interface takes a carrier wavelength (manual: *Ionisation rates*). Evaluating it at the
pump carrier on a multiband field does not supply the VUV cross sections proposed here.
A future solprop proposal can use the interface independently of its choice of band
split or propagation algorithm. No unaccepted downstream design is a prerequisite for
this OEP.

**In scope:**
- per gas (He, Ne, Ar, Kr, Xe) and polarisation, $\sigma^{(2)}(\omega)$ from the
  two-photon threshold to just below the first one-photon resonance, where curated data
  support it, with:
  - the single He⁺ ²S₁/₂ threshold, or both ²P₃/₂ and ²P₁/₂ thresholds for Ne–Xe;
  - labelled support regions;
  - an uncertainty;
  - a resonance mask;
- the conversions a propagation code needs (rate, loss coefficient, the relation to
  Im γ, the broadband form);
- in Phase 2, the bound even-parity two-photon levels below threshold (positions,
  strengths where curated, and the photoionisation of the excited level), as a complex
  two-photon response compatible with OEP 5's resonance limits;
- tests, docstrings, a manual page, a gallery example and the register entries.

**Out of scope:**
- molecules;
- three- and higher-photon ionisation (covered by the strong-field rates and the curated
  $\sigma^{(N)}$ datasets);
- one-photon absorption above $E_1$;
- any refit of the hyperpolarizability;
- a fitted smooth form of $\sigma^{(2)}$ (section 11);
- the propagation code's band split, broadband algorithm and policy outside the support:
  these are the paired solprop proposal's. optmat states facts.

### 3. Context

- **Manual:**
  - {doc}`../../two_photon`: the curated two-photon data compared on one footing, with the
    thresholds and what is missing (TC-902);
  - *Ionisation rates*: the strong-field models, their 350 nm lower limit and its
    `ModelRangeWarning`, *Multiphoton saturation at 248.6 nm* and *Open questions*.
  - *Third-order nonlinearity*: *Resonance limits* and *The support record*, OEP 5.
  - *The curated data*: the vocabulary.
  - *Propagation*: the surface a propagation code builds on, OEP 1.
- **Modules:**
  - `optmat.gas.ionisation_rates` (`ionisation_rate_model`);
  - `optmat.gas.nonlinear` (`HyperpolarizabilitySupport`: $E_1$, $E_2$, the ionisation
    energy, `band_edge`);
  - `optmat.identity` (`Identified`);
  - `optmat.vocabulary`.
- **Data:**
  - `multiphoton_ionization_cross_section` with `photons = 2`;
  - the two-photon quantities added for this work (TC-872):
    - `two_photon_transition_strength`;
    - `two_photon_excitation_cross_section`;
    - `excited_state_photoionization_cross_section`;
    - `multiphoton_ionization_polarization_ratio`;
  - `nist_asd_2026`: the neutral ionisation energies, the ions' ²P₁/₂ levels, and the
    neutral levels of He to Xe as returned (TC-443, TC-862).
- **Register:**
  - TC-856: McKenna and van der Hart's Ne and Ar cross sections, and the flux convention
    $\Gamma = \sigma^{(N)}\Phi^N$ shared with Uiterwaal et al.
  - TC-853: multiphoton saturation at 248.6 nm.
  - TC-862, TC-865 and TC-868: the resonance limits, and the rule that optmat states
    facts.
  - TC-872 to TC-902: the vocabulary, the curation of the two-photon literature, and its
    comparison on one footing (the manual page *Two-photon ionisation and absorption:
    the data*, TC-902).
- **Prior OEPs:** OEP 1 (identity, digests, the declared surface); OEP 5 (the support
  record and its band edge $\min(E_1, E_2/2)$).

### 4. Theory

SI throughout, in the dilute-gas limit $n \simeq 1$, without local-field corrections.
$\omega$ is a photon's angular frequency, and $\Omega = \omega_a +
\omega_b$ a pair's sum. The fields are real, $E(t) = E^{+}(t) + E^{-}(t)$, with $E^{+}$
the positive-frequency part ($E^{-} = E^{+*}$). A monochromatic wave $E_0\cos\omega t$
has $E^{+} = (E_0/2)e^{-i\omega t}$, intensity $I = \tfrac12\varepsilon_0 c E_0^2$ and
photon flux $\Phi = I/\hbar\omega$. The analytic signal common in propagation codes is
$2E^{+}$.

**The cross section and the rate.** The generalised two-photon ionisation cross section
is defined by the per-atom rate

$$
W = \sigma^{(2)}\,\Phi^2 ,
$$

the convention of Uiterwaal et al. [2] and McKenna and van der Hart [1] (their eq. 3;
TC-856). With $\Phi$ in cm⁻² s⁻¹, $\sigma^{(2)}$ is in cm⁴ s; 1 cm⁴ s = 10⁻⁸ m⁴ s. No
factor 2! enters a cross section so defined.
- **Other conventions:** McCown et al. [3] print a coefficient $\alpha$ in cm⁴ W⁻¹
  defined by $dn/dt = \alpha[\mathrm{Xe}]I^2/\hbar\omega$ (their Eq. 9), so
  $\sigma^{(2)} = \alpha\hbar\omega$ (TC-884). The bound-state literature's "two-photon
  cross sections" are frequency-integrated (TC-872). Each is converted once, by its
  record's stated definition.
- **Intensity dependence:** a non-perturbative calculation's cross section is
  $\Gamma/\Phi^2$ at one stated intensity (10¹² W cm⁻² for [1]). It approximates
  the lowest-order value only where higher-order and field-dressing effects are
  negligible. The stored model records that intensity and the evidence for its
  perturbative interpretation; an intensity-independent interface does not establish
  validity at arbitrary intensity.

**Energy loss.** Each event removes two photons. A quasi-monochromatic wave in a gas of
number density $N$ therefore loses intensity as

$$
\frac{dI}{dz} = -2\hbar\omega\,N W = -\beta I^2,\qquad \beta = \frac{2N\sigma^{(2)}}{\hbar\omega},
$$

the form of the solids' `two_photon_absorption_coefficient`; the gas's $\beta$ is derived,
never stored. Here $N$ is the remaining neutral ground-state density. For an ionic final
channel $j$ of threshold energy $I_{p,j}$, the electron carries
$2\hbar\omega-I_{p,j}$; the upper fine-structure channel also leaves energy in the ion.
A total cross section without channel branching cannot determine that energy partition. A
propagation code that removes the field's energy through a loss current must remove
$2\hbar\omega$ per event, not $I_p$.

**The relation to Im γ.** optmat's $\gamma$ is Shelton and Rice's Taylor-series
hyperpolarizability [4]. In their eq. 3, the induced dipole at the sum frequency is
$(1/6)K^{(3)}\gamma E(\nu_1)E(\nu_2)E(\nu_3)$. For the intensity-dependent index
$(-\nu;\nu,\nu,-\nu)$, $K^{(3)} = 3/4$ ([4], Table I). For a linearly polarised
$E_0\cos\omega t$, the dipole at $\omega$ is therefore $\mu_\omega = \gamma E_0^3/8$, with
$\gamma = \gamma_{zzzz}(-\omega;\omega,\omega,-\omega)$. The time-averaged power absorbed
per atom is $\tfrac12\omega\,\mathrm{Im}\,\mu_\omega\,E_0 =
(\omega/16)\,\mathrm{Im}\gamma\,E_0^4$. At two photons per event,

$$
W = \frac{\mathrm{Im}\,\gamma\;E_0^4}{32\hbar},\qquad
\sigma^{(2)} = \frac{\hbar\omega^2\,\mathrm{Im}\,\gamma}{8\,\varepsilon_0^2c^2}
\quad(\gamma\ \text{in C m}^4\,\text{V}^{-3}).
$$

In the perturbation-series convention, $\gamma_P = \gamma/6$, the first reads
$W = (3/16)\,\mathrm{Im}\gamma_P E_0^4/\hbar$. A factor of 6 separates the two
conventions; optmat uses the Taylor series throughout. This identifies the absorptive
contribution from the two-photon channel in the perturbative, monochromatic, linearly
polarised limit. It is not a formula for circular polarisation using the same scalar
$\gamma_{zzzz}$, nor a claim that all imaginary third-order response at arbitrary
frequencies is two-photon ionisation.

**The thresholds.** Direct two-photon ionisation of the ground state needs
$2\hbar\omega \ge I_p$. For Ne, Ar, Kr and Xe the ion has two fine-structure levels,
²P₃/₂ (the ground level) and ²P₁/₂. Between the two thresholds the ²P₁/₂ Rydberg series
autoionise, so $\sigma^{(2)}$ there is a dense series of resonances. Helium has one
threshold. From `nist_asd_2026` (the ionisation energy, plus the ion's ²P₁/₂ level,
halved):

| | $I_p/2$ (lowest ionic channel) | upper ²P₁/₂ threshold | first one-photon resonance $E_1$ | wavelength at OEP 5's $E_2/2$ |
|---|---|---|---|---|
| He | 12.294 eV, 100.85 nm | — | 21.218 eV, 58.43 nm | 120.28 nm |
| Ne | 10.782 eV, 114.99 nm | 10.831 eV, 114.48 nm | 16.671 eV, 74.37 nm | 133.49 nm |
| Ar | 7.880 eV, 157.34 nm | 7.969 eV, 155.59 nm | 11.624 eV, 106.67 nm | 189.36 nm |
| Kr | 7.000 eV, 177.13 nm | 7.333 eV, 169.08 nm | 10.032 eV, 123.58 nm | 216.67 nm |
| Xe | 6.065 eV, 204.43 nm | 6.718 eV, 184.55 nm | 8.437 eV, 146.96 nm | 256.02 nm |

The table is computed from the record's printed values with $hc$ = 1.239841984×10⁻⁴
eV cm; the implementation generates and tests it. Some consequences:
- Argon at the F₂ laser lines (about 157.6 nm) is just below its threshold.
- Krypton at 193 nm and xenon at 248 nm are below theirs.
- Xenon at 193 nm lies between its two thresholds.

Phase 1 covers the window $[I_p/2, E_1)$: at and beyond $E_1$ intermediate one-photon
resonances enter, so this OEP does not qualify that region. This is not an assertion of
continuous one-photon absorption at every energy above the first discrete line.

**Energy axes of the curated calculations.** The calculations are of three kinds:
- **LS-coupled calculations** [1, 5] have one threshold and resonances at their own
  positions (e.g. argon's 4s ¹P one-photon resonance), which differ from the measured,
  J-resolved levels by some tenths of an eV.
- **Some calculations shift their axis.** Saenz and Lambropoulos [6] place their
  thresholds at the experimental values after a 0.69 eV error in the ground-state energy.
  Karamatskou and Santra [7] shift by about 0.3 eV to match another calculation.
- **Relativistic and MQDT calculations** [8, 9] place the fine-structure thresholds
  correctly.

This OEP does not re-shift a source's axis. It takes each reference
curve as published and labels as *unsupported detail* the energies where the source's
resonance positions or thresholds differ from NIST's by more than the curve's structure
allows.

**Broadband fields.** For a weak coherent pulse, fixed linear polarisation and negligible
ground-state depletion, the lowest-order amplitude to a continuum state of energy
$\hbar\Omega$ above the ground state is an integral over the pairs of spectral components $\omega_a + \omega_b = \Omega$. Let the two-photon matrix element
vary slowly with the split $\omega_a - \omega_b$ at fixed $\Omega$, separately in each
final channel. Staying far from intermediate resonances is necessary but does not by
itself establish this approximation for a broad band: cancellations and channel
interference also matter. The degenerate cross section contains no dependence on the
frequency split, and cannot test that approximation (Q9). The amplitude is then the
degenerate one times the Fourier transform of $E^{+}(t)^2$. The ionisation probability of a pulse follows as

$$
P = \frac{1}{2\pi}\int_0^\infty d\Omega\;
\frac{16\,\varepsilon_0^2c^2\,\sigma^{(2)}(\Omega/2)}{\hbar^2\,\Omega^2}\;
\bigl|S_2(\Omega)\bigr|^2,\qquad S_2(\Omega) = \int dt\;e^{i\Omega t}\,E^{+}(t)^2 .
$$

In the monochromatic limit, $|S_2|^2 \to (E_0^4/16)\,2\pi T\,\delta(\Omega-2\omega)$ over a
duration $T$, which recovers $P = T\sigma^{(2)}\Phi^2$. The onset of two-photon ionisation
in this approximation depends on the spectrum of $E^{+2}$. A strictly band-limited
field on $[\omega_{\min},\omega_{\max}]$ has pair support inside
$[2\omega_{\min},2\omega_{\max}]$. The appreciable yield depends on spectral shape
and phase; there is no universal half-bandwidth shift. For a fluctuating source the
ensemble average of $|S_2|^2$ is required, not the squared mean spectrum or intensity.
The integral is a perturbative yield, not a saturated probability when it approaches 1.

The energy removed per event is $\hbar\Omega$, drawn from the pairs that make it. Three
treatments follow:
1. **The centroid:** $\sigma^{(2)}$ at the band's centroid frequency, the narrowband
   limit. It requires little variation of the response across the band, including
   distance from thresholds and all resonances; distance from $E_1$ alone is insufficient.
2. **The pair spectral density:** the integral above. Exact in the stated approximation,
   and needs $S_2$ on each step.
3. **A spectral factor, as an integral identity only:** define
   $K(\Omega)=16\varepsilon_0^2c^2\sigma^{(2)}(\Omega/2)/(\hbar^2\Omega^2)$
   on the qualified positive-frequency domain. Any chosen extension and factor
   $|H|^2=K$ gives $P=\int |h*E^{+2}|^2dt$ by Parseval. This does not determine a
   physical time-resolved ionisation rate: its time dependence changes with the
   undetermined phase of $H$. Parseval establishes neither causality nor the existence
   of a minimum-phase completion with these exact spectral gaps. A causal approximation
   would require additional assumptions outside the measured band. Phase 1 offers neither
   a causal filter nor a loss current derived from it.

The corresponding energy removed per initially neutral atom, within the same
approximation, is

$$
\Delta\mathcal E = \frac{1}{2\pi}\int_0^\infty
\hbar\Omega\,K(\Omega)|S_2(\Omega)|^2\,d\Omega.
$$

A positive integrated yield alone does not specify a current that removes this energy
from the contributing frequency pairs. optmat documents these integrated constraints;
a propagation algorithm must establish its own time dependence and energy conservation.
Only pair spectra within the reference's qualified domain can be evaluated without an
additional model. Multiplication by $K=0$ in a data gap would silently discard unknown
loss and is not permitted. The scalar formulas are for fixed linear polarisation;
circular monochromatic data do not determine a general polarisation tensor.

**Below the threshold.** A ground-state atom below $I_p/2$ can still be excited by two
photons, to the even-parity levels:
- $J = 0, 2$ for a scalar (one linear) field;
- also $J = 1$ for general polarisation at unequal frequencies (OEP 5, section 4).

The lowest is OEP 5's $E_2$. Their series converges on the ionisation threshold. A
further photon may ionise a populated level, with instantaneous rate
$W_{n\to i}=\sigma_{n\to i}(\omega)I/(\hbar\omega)$ where that cross section is
known. Whether it does so during a pulse depends on the remaining fluence, branching,
and competing relaxation, not just its peak intensity. The measured Xe 6p cross
sections [10] and the Kr population model [21] concern particular levels and wavelengths;
they do not establish prompt ionisation of every bound level in every gas.

Two-photon excitation and subsequent ionisation are separate events: the first removes
two photons and creates an excited neutral, the second removes a third and creates an
electron. A model that sets the *direct continuum* cross section to zero below $I_p/2$
is correct for that channel, but leaves bound absorption and (2+1) ionisation unmodelled.
Phase 2 must distinguish these populations and energy transfers. Field dressing,
stimulated emission and coherent population dynamics may require more than rate equations.

The strength density of high Rydberg members may approach the continuum in an appropriate
channel-resolved, energy-averaged limit. That must be tested from data, not imposed by
extrapolating the few curated low-lying levels. Finite pulse resolution may leave lines
unresolved; it does not turn bound states into an ionisation continuum.

**A complex two-photon response (Phase 2).** Pair denominators in the sum over states
([4], eq. 14; OEP 5, section 4) motivate a schematic positive-frequency resonant part

$$
\chi_2^{\mathrm{res}}(\Omega) =
\sum_n \frac{S_n}{\Omega_{ng}-\Omega-i\Gamma_n/2}
+ \int_{I_p/\hbar}^\infty
\frac{s(\Omega')\,d\Omega'}{\Omega'-\Omega-i0}.
$$

This is a candidate representation, not yet a defined susceptibility or evaluator.
If its units and absorptive part are chosen to match the Taylor-series $\gamma$ on the
linearly polarised degenerate diagonal, section 4 requires

$$
\operatorname{Im}\chi_2^{\mathrm{res}}(\Omega)
= \frac{32\varepsilon_0^2c^2}{\hbar\Omega^2}
\sigma^{(2)}(\Omega/2),\qquad
\pi s(\Omega)=\operatorname{Im}\chi_2^{\mathrm{res}}(\Omega)
$$

for the unbroadened continuum term. Thus $s$ scales with $\sigma^{(2)}/\Omega^2$,
not with $\sigma^{(2)}\Omega^2$. This fixes an absorptive normalisation on that
slice, not the full multifrequency response. Residues, tensor factors, other frequency
permutations and negative-frequency terms still need a derivation (Q9). A finite
continuum window cannot fix the dispersive real part without missing spectral weight
and background terms; a one-variable dispersion relation for the full degenerate Kerr
response must not be assumed from causality alone.

Nor is $\chi_2^{\mathrm{res}}$ the spectral factor $H$: absorption is proportional
to $\operatorname{Im}\chi_2^{\mathrm{res}}$, whereas the integral identity needs
$|H|^2$. For one isolated pole with real strength $S$,
$\operatorname{Im}\chi=S(\Gamma/2)/(\Delta^2+(\Gamma/2)^2)$ but
$|\chi|^2=S^2/(\Delta^2+(\Gamma/2)^2)$; even there their ratio depends on the
strength and width. Several poles also generate interference terms in $|\chi|^2$.
The cross section does not supply the continuum amplitudes or their phases.

Physical population decay contributes half its rate to coherence decay; pure dephasing
is separate. In the displayed pole, $\Gamma/2$ is the coherence decay rate, in s⁻¹.
Pulse bandwidth belongs to the driving field and its statistics, not to an intrinsic
level width. Photoionisation can add a field-dependent decay rate, but then a fixed
linear convolution is insufficient in general. Phase 2 must reconcile these definitions
against [10, 20, 21] before choosing a response model (Q7, Q9). Its real part, if ever
offered, also needs a subtraction or matching prescription to avoid counting the
existing Kerr response twice (Q8).

### 5. Primary sources and data

Curated before this OEP (TC-872 to TC-902) and compared on one footing on the manual
page *Two-photon ionisation and absorption: the data*. All are in hand as source
records.

**Two-photon ionisation, $\sigma^{(2)}$:**

| Gas | Calculations (vector figures unless stated) | Measurements | Notes |
|---|---|---|---|
| He | Saenz & Lambropoulos [6] (partial ¹S, ¹D; L and V); Mihelič & Horvat [11]; Nikolopoulos & Lambropoulos [12] | Sato et al. [13], FEL: one point in the window (61.4 nm), one on the 1s2p resonance, two above $E_1$ | three calculations, within about 15 % (TC-902) |
| Ne | McKenna & van der Hart [1] | — | McGuire [14], circular and $R$ (raster, check only) |
| Ar | McKenna & van der Hart [1]; Petrov et al. [5] (L and V) | — | within 10 % with core polarisation; Starace & Jiang [15] (raster) bracket them; McGuire [14] |
| Kr | L'Huillier, Tang & Lambropoulos [9], MQDT (raster); McGuire [14] (raster) | — | **no vector or tabulated source** |
| Xe | Lagutin et al. [8], linear and circular, from the ²P₁/₂ threshold; Karamatskou & Santra [7], above the window | McCown et al. [3], 193 nm (between the thresholds) | L'Huillier & Wendin [16] and L'Huillier et al. [9] (raster) are 0.3–0.6 of Lagutin's |

Under data rule 8, every curve read from a figure is `use = "check"`. A curve read from a
vector drawing may become the model's reference only with the maintainer's explicit
approval, recorded in the register; a raster curve never may (Q1).

**Bound two-photon levels and their ionisation:**
- **Positions:** `nist_asd_2026`, from the neutral level lists already held as returned
  (Phase 2 selects the even $J = 0, 2$ series).
- **Strengths:**
  - Hättig, Christiansen & Jørgensen [17] (He, Ne, Ar; $\delta^{TP}$);
  - Henriksson et al. [18] (He to Xe, relativistic);
  - Shekhtman et al. [19] (Kr, 192–217 nm);
  - Kröll & Bischel [10] (Xe 6p, measured);
  - Drag et al. [20] (Xe 6p′, measured);
  - Dakka et al. [21] (Kr 5p, measured).
- **Excited-level photoionisation:** Kröll & Bischel [10]; Dakka et al. [21];
  Geohegan et al. [22].

**Validators** (held out of any reference curve):
- the other calculations of each gas;
- Sato's He and McCown's Xe measurements;
- Uiterwaal et al. [2] at 248.6 nm, for the convention ($\Gamma = \sigma^{(N)}\Phi^N$),
  not for $\sigma^{(2)}$.

**Missing, and requested:**
- independent near-threshold Kr calculations and any Kr measurement;
- the relativistic Ne–Xe calculations of Fink & Johnson (not curated);
- the 193 nm Xe measurements of Bischel et al. and Hodges et al. (not curated);
- any measured Ar $\sigma^{(2)}$;
- the bound strengths of Ne 3p, Ar 4p and He 1sns/nd above 1s2s.

**Uncertainty and model spread.** Preserve three distinct kinds of evidence:
- the source's own accuracy or uncertainty, including its stated level and limitations;
- digitising uncertainty, including energy calibration where it affects a steep curve;
- spread between comparable independent calculations over an explicitly named overlap.

Gauge variants and repeated drawings of one calculation are diagnostics, not independent
replicates. Raster curves may inform comparisons but never become the reference. A spread
is an empirical envelope, not a confidence interval or a bound on the true answer. Report
its lower and upper ratios to the reference, the sources included and any exclusion's
reason. Do not use a maximum with digitising error to erase a larger stated source
uncertainty (e.g. Lagutin's gauge-based accuracy estimate, TC-880).

For Ne the reference's accuracy remains *single source, unstated*; a digitising estimate
is still available and is not the model's total uncertainty. An absent numerical
uncertainty is represented explicitly, never as zero (Q6).

### 6. Proposed interface

A new module, `optmat.gas.two_photon`. Phase 1:

| Name | Meaning | Units |
|---|---|---|
| `two_photon_ionisation_model(gas, polarization="linear")` | the model of one gas | |
| `.cross_section(omega)` | $\sigma^{(2)}$ at photon angular frequency `omega`; exactly 0 below the lower threshold | m⁴ s |
| `.pair_cross_section(Omega)` | the same, indexed by the pair's sum ($\sigma^{(2)}(\Omega/2)$), for the broadband forms | m⁴ s |
| `.uncertainty(omega)` | immutable record of source uncertainty, digitising uncertainty and model spread, with their levels, availability and provenance | relative values: 1 |
| `.rate(omega, intensity)` | $W = \sigma^{(2)}(I/\hbar\omega)^2$ | s⁻¹ |
| `.loss_coefficient(omega, number_density)` | $\beta = 2N\sigma^{(2)}/\hbar\omega$ | m W⁻¹ |
| `.thresholds` | one He⁺ ²S₁/₂ threshold, or two ²P thresholds for Ne–Xe, with final-state labels, uncertainties and NIST datasets | rad s⁻¹ |
| `.support` | the labelled regions (supported, unsupported detail, unavailable, below threshold), each with its datasets, reasons and register entries | rad s⁻¹ |
| `.reference` | the reference dataset(s) and the validators | |
| `.identity`, `.digest` | OEP 1 | |

**Warnings and refusals:**
- **Finite domain:** reference coverage is the intersection of the source's actual
  intervals with $[I_p/2,E_1)$. No endpoint clamping, extrapolation, interpolation across
  a missing segment, or manufactured value between the physical and calculated
  thresholds. Missing data is *unavailable*, not a physical zero.
- **Outside coverage:** positive frequencies below the lower threshold return zero for
  direct continuum ionisation only. At or above the threshold, a request outside the
  reference coverage (including $E_1$ itself) raises a documented `OptmatError` subclass.
  A warning alone cannot define a value where no continuation is specified.
- **Unsupported detail:** a finite source value may be returned with
  `ModelRangeWarning` where the source exists but has an unqualified resonance or axis
  mismatch. Its region includes a reason and source intensity. A clipped or missing
  pole segment has no reconstructible value and remains unavailable. Region masks are
  metadata, not a new intensity-independent resolution of a resonance.
- **Polarisation:** circular is available only where an approved vector reference gives
  it (Xe). Other requests raise an `OptmatError` subclass; no universal linear/circular
  conversion is inferred from a raster ratio (Q4).
- **Krypton:** the constructor raises a missing-model error until an eligible reference
  is curated and approved (Q2); its thresholds remain documented.
- **Inputs and arrays:** frequencies must be finite and strictly positive; intensity and
  density finite and nonnegative, in SI. Arrays broadcast consistently with the existing
  public surface. An array containing an unavailable point raises rather than filling
  with zero or NaN. Metadata inspection remains possible before evaluation; the zero
  continuum value below threshold carries no claim of zero total absorption.

A sketch (fenced as `text` until it runs; values to be generated):

```text
>>> from optmat.gas.two_photon import two_photon_ionisation_model
>>> ar = two_photon_ionisation_model("Ar")
>>> [round(2 * math.pi * c / w * 1e9, 2) for w in ar.thresholds.omega]
[157.34, 155.59]
>>> ar.reference.datasets  # illustrative candidate, subject to Q1 approval
('mckenna_2004/fig5_sigma2_Ar',)
>>> ar.loss_coefficient(omega_155nm, number_density_1bar)  # m/W
...
```

**Composition:**
- **With the strong-field rates:** no double counting. The PPT/ADK models describe the
  pump in their qualified regime; this model describes degenerate two-photon continuum
  ionisation. Applying both to the same component risks overlapping contributions, and
  a band split is a downstream modelling choice. It does not prove completeness:
  $(E_p^+ + E_v^+)^2$ contains $2E_p^+E_v^+$ as well as the two same-band terms. Mixed
  pairs, interference and field dressing require nondegenerate amplitudes or a separately
  justified neglect; the degenerate curve alone does not supply them (Q9).
- **With OEP 5:** this model's lower threshold $I_p/2$ lies above OEP 5's band edge
  $\min(E_1, E_2/2)$. Between them lie the bound two-photon levels (Phase 2). Phase 1's
  support labels that interval *below threshold, bound resonances not modelled*, and its
  value there is 0.

### 7. Implementation strategies

**Representation.**
- *(A) The curated reference curves, interpolated* (recommended for Phase 1). For each
  gas and polarisation, a selection file (`data/gas/two_photon.toml`) names:
  - the reference dataset, with the reason;
  - the validators;
  - the support regions with their labels;
  - the uncertainty per region.

  The model interpolates the reference's vertices, log-linearly in $\sigma$ against photon
  energy, within a continuous positive-valued segment. Zeros, clipped poles and
  disconnected paths must be identified before taking logarithms: no arbitrary positive
  floor and no interpolation across a clipped resonance. If only partial channels are
  available (He), interpolate each on its own grid, sum the required disjoint channels
  on their common coverage, and never count a total alongside its partials. Duplicate
  energy coordinates with distinct values need an explicit resolution, not sorting into
  an arbitrary branch. Nothing is fitted; the source grid and every mask remain traceable.
- *(B) A smooth near-threshold form with fitted parameters.* This is a new fitting
  methodology needing its own justification. It gains little where the curves are flat
  off resonance (Ne and Ar are within 10 % of a constant over most of the window [1]),
  and it would smooth away the structure between the fine-structure thresholds. Not
  proposed (section 11).
- *(C) A pole-plus-background complex response* $\chi_2(\Omega)$ (section 4). A candidate for
  Phase 2's resonant response, subject to Q7–Q9. Phase 1 supplies continuum absorption
  evidence for that work, not a uniquely determined causal amplitude or complex response.

**Where the numbers live.** As for OEP 5: the selection file is reviewed data, and a
stored derivation (`data/gas/models/two_photon/<gas>.json`, written by a tool and checked
by a freshness test) holds the reference grid, thresholds, regions and uncertainties. The
model reads only the stored derivation, so loading is fast. The derivation fingerprints
its selected source records, grids, thresholds, validators and selection file; a freshness
test detects stale inputs. OEP 1 identity also includes evaluator versions and immutable
fields, so a change to interpolation or conversion changes the digest even with unchanged
data. Provenance and region diagnostics are included explicitly; a bare file hash does
not replace OEP 1's identity contract.

**Vocabulary and data.** No new quantity in Phase 1: TC-872 added what the curation
needed. Phase 2 may need a level-width quantity (Q7).

**Cost.** Small JSON files; no new dependency.

### 8. Verification and validation

- **Conventions:**
  - McCown's coefficient against $\sigma^{(2)}$, by its record's definition;
  - the cm⁴ s to m⁴ s conversion;
  - $\beta$ against a direct integration of $dI/dz = -2\hbar\omega NW$;
  - the Im γ relation, against an effective pair transition with a specified coupling
    and decay, or a three-level electric-dipole model with an intermediate odd state;
    a centrosymmetric two-level electric-dipole atom alone cannot supply this process;
  - the broadband formula's monochromatic limit, numerical Parseval identity and
    energy-weighted integral on a weak Gaussian pulse; these check normalisation, not a
    unique time-resolved rate. Two-frequency pulses test convolution factors explicitly;
    physical nondegenerate validation remains Q9.
- **Thresholds:** reproduce `nist_asd_2026`'s printed values and the table of section 4.
- **The reference reproduces its record:** the stored input vertices exactly, and their
  evaluated values within floating-point rounding on every available segment; every
  printed value the record checks (e.g. Ar within 8 % of 1.2×10⁻⁵⁰ cm⁴ s from 8.0 to
  11.4 eV; Ne within 10 % of 1.65×10⁻⁵¹ cm⁴ s from 11 to 16 eV; TC-856).
- **Validators:**
  - **Ar:** McKenna against Petrov, length and velocity.
  - **He:** Saenz against Mihelič and Nikolopoulos, and against Sato's measured value at
    61.4 nm (20.2 eV, inside the window).
  - **Xe:** Lagutin against independent calculations only over their common coverage.
    McCown's randomly polarised, finite-linewidth 193 nm measurement lies below Lagutin's
    grid. Keep it as context and a future validator of that missing region, with its
    polarisation and spectrum accounted for; do not extrapolate Lagutin to pass it.
  - **Ne:** the single source, against McGuire's circular curve and $R$ (raster, a check
    only).

  Each comparable ratio and its overlap are tabulated and recorded in the register;
  noncomparisons state the missing coverage or experimental information. Sato's He point
  is a conditional experimental check: the assumed duration, intensity dependence and
  uncorrected photon statistics preclude using its printed error as a full acceptance
  tolerance (Q5). The reported model-spread envelope contains every calculation selected
  for that comparison, and keeps source uncertainty separate.
- **Domain and interpolation:** scalar and broadcast inputs; exact endpoints; unavailable
  intervals and array failures; positivity; zero below threshold only; pole masks;
  reconstruction of disjoint partial channels; no extrapolation or arbitrary log floor.
  Numerical interpolation tolerances follow floating-point error; comparison tolerances
  follow the named sources and digitising resolution, not an arbitrary common percentage.
- **Identity and freshness:** golden digests; a change to the selection, the reference or
  a validator changes the digest; the stored derivation equals a fresh one; the installed
  wheel contains the files.
- **Later validators:**
  - time-dependent Schrödinger-equation runs (tRecX) of argon or helium at a few photon
    energies, held as tier-3 validators below the published calculations (manual:
    *Ionisation rates*);
  - stationary or single-active-electron calculations a paired solprop proposal
    documents.
- **Parity with solprop:** the paired proposal reads the model through its adapter and
  records its digest.

### 9. Documentation

Extend the existing `docs/two_photon.md`, *Two-photon ionisation and absorption*, with:
- the conventions and derivations of section 4;
- the per-gas thresholds;
- the curated data and their comparison (TC-902);
- the model, its support and uncertainty;
- what is missing.

Links from *Ionisation rates* (*Open questions*), *Third-order nonlinearity* (the support
record) and *Propagation* (the surface table).

A gallery example shows the loss coefficient against wavelength for argon and helium at
1 bar. A prescribed-pulse attenuation example
uses $I(z)=I(0)/(1+\beta I(0)z)$ at fixed neutral density, stating the path length,
temperature and pulse profile. It illustrates two-photon loss, not a prediction of
resonant dispersive-wave generation or an energy ceiling without a propagation model:
- argon near 155 nm (above both thresholds);
- helium near 97 nm. Helium's threshold is at 100.85 nm, so 105 nm is below it, near the
  1s4s/4d two-photon resonances: a Phase 2 example.

### 10. Compatibility and impact

- New module and data. No recommended value, default, stored fit or evaluated quantity of
  any existing model changes.
- The public-API snapshot gains the new names; the release test changes.
- The strong-field rate models are unchanged; their manual page gains the
  no-double-counting note.
- No solprop parity fixture changes until the paired proposal uses the model.

### 11. Alternatives considered

- **A fitted near-threshold form** (option B). A new fitting methodology, small gain, and
  it would blur the fine-structure region.
- **One constant per gas** (e.g. argon's 1.2×10⁻⁵⁰ cm⁴ s). Simple, and within 10 % over
  most of the window for Ne and Ar. However, it hides the structure near $E_1$ and between
  the thresholds, and has no basis for Kr and Xe.
- **Re-shifting LS-coupled curves onto NIST's levels.** Not proposed: one shift cannot
  resolve channel splitting or independently misplaced resonances. A future model could
  explicitly transform an unchanged source grid with a justified, registered mapping;
  that would be a model decision, not an erratum to the printed source.
- **Extending the strong-field rate to short wavelengths.** Evaluating a pump-carrier
  rate on a multiband field does not recover a calibrated VUV two-photon cross section.
  Changing its carrier is not validation either; its existing 350 nm qualification
  limit stays.
- **Using raster-read curves (Kr) as the model.** Rejected by data rule 8; Q2 decides
  what optmat offers for krypton.

### 12. Phases

1. **The continuum model.** Deliverables:
   - the selection file and the stored derivation;
   - `optmat.gas.two_photon` with the interface of section 6;
   - the thresholds;
   - the tests of section 8;
   - the manual page and the gallery example;
   - the register entries, including the maintainer's rule-8 approval of each reference
     curve;
   - a changelog line.

   It verifies section 8 except solprop parity.
2. **Bound two-photon resonances.** Deliverables:
   - the even-parity level series from the held NIST level lists;
   - the strengths and widths from the curated bound data, with their conventions
     reconciled ($\delta^{TP}$ against $\sigma_0^{(2)}$);
   - a reviewed derivation and interface for the response, including Q7–Q9, before
     implementation of a complex evaluator;
   - level-specific absorption and photoionisation inputs for downstream populations;
     an ionisation fraction additionally requires a pulse and a dynamical model;
   - distinct bound-channel support; Phase 1's direct continuum still vanishes below
     threshold. Missing higher series members remain explicitly unmodelled.

   It checks the line strengths against every comparable measured value (Kröll, Drag,
   Dakka), preserving photon-statistics and frequency-integration conventions. It tests
   continuity of the strength density only if enough high-lying levels are available;
   otherwise it records that limit as unverified. Q9 is a blocker for accepting the
   complex-response deliverable, not permission to invent a closure during implementation.
3. **(Optional) Consistency with γ.** The real part of $\chi_2$ against OEP 5's fitted
   $\gamma$ below $E_2/2$, and whether a two-photon-resonant Kerr term should be offered.
   Only if Phase 2 shows the data support it.

### 13. Open questions

- **Q1 (the references; Phase 1).** Which vector curve is each gas's reference, by the
  maintainer's rule-8 approval:
  - He: Saenz and Lambropoulos's ¹S + ¹D, or Mihelič and Horvat's;
  - Ne: McKenna and van der Hart (the only one);
  - Ar: McKenna and van der Hart, or Petrov et al.;
  - Xe: Lagutin et al.

  Proposed references must be identified by dataset and gauge, and each approval must
  precede its use in an evaluated model. Review status is not that approval.
- **Q2 (krypton; resolved recommendation).** No evaluated continuum model until an
  eligible reference is curated and approved. The raster data remain comparisons. No
  interpolated Ar–Xe constant or raster-based fallback is offered.
- **Q3 (between the fine-structure thresholds; Phase 1).** With the present candidate
  references, Xe is unavailable between its lower threshold and Lagutin's first vertex;
  Kr has no model. LS-coupled Ne and Ar curves do not resolve ionic fine structure, and
  any covered values in that region need an explicit *unsupported detail* label. Any
  newly proposed curve or averaging scheme requires a separate documented decision.
- **Q4 (circular polarisation; resolved recommendation).** Only an approved Xe circular
  reference is offered. McGuire's raster ratio may be checked, but supplies no model
  factor for the others.
- **Q5 (Sato's measurement; resolved review finding).** A conditional experimental check,
  not a quantitative lowest-order calibration from its printed errors alone: the pulse
  duration is assumed, the cross section depends on intensity, and the paper does not
  correct photon statistics (TC-876). A quantitative comparison requires the experimental
  intensity and coherence information; no correction factor is guessed.
- **Q6 (uncertainty; resolved recommendation).** Publish source and digitising uncertainty
  separately from the nonstatistical model-spread envelope (section 5). An unstated
  source uncertainty stays unstated even where a spread is available.
- **Q7 (Phase 2 widths; resolve in its derivation).** Prefer curated intrinsic lifetimes,
  dephasing where known, and level-specific photoionisation cross sections. The driving
  pulse determines bandwidth and field-induced decay in the dynamical calculation.
  Specify whether every quoted width is in photon or pair frequency, angular or ordinary
  frequency, and full or half width; any new curated quantity needs vocabulary first.
  Dakka's submitted manuscript also has an unresolved factor in its photoionisation
  rate definition (TC-892); resolve it before using that rate as an absolute validator.
- **Q8 (Phase 2 real part; optional Phase 3).** Do not add a dispersive term to OEP 5's
  fitted $\gamma$ until the missing spectral weight, other response terms and a
  double-counting subtraction are justified. Defer provision of the real part if the
  data cannot support these.
- **Q9 (broadband response; acceptance blocker for the complex evaluator).** Derive or
  obtain the nondegenerate, channel- and polarisation-resolved amplitudes needed to
  qualify the pair-frequency approximation. Define the units and residues of the bound
  response, the population and coherence equations, and their energy balance. Demonstrate
  why the proposed dynamical response follows from those amplitudes; a positive loss
  spectrum or a choice of spectral-factor phase is insufficient. The current data and
  this review establish the integrated normalisation, not that closure. Deferred to
  Phase 2: before implementing its complex evaluator, supply and review that theory.
  If it cannot be established, restrict Phase 2 to curation and report the evaluator
  as deferred to a later proposal; acceptance of this OEP does not resolve Q9.

### 14. References

All curated as source records (the record's id in brackets).

1. C. McKenna and H. W. van der Hart, "Multiphoton ionization cross sections of neon and
   argon," *J. Phys. B* **37**, 457 (2004). doi:10.1088/0953-4075/37/2/013
   (`mckenna_2004`).
2. C. J. G. J. Uiterwaal et al., "Generalized multiphoton-ionization cross sections of
   the rare gases for 500-fs, 248.6-nm pulses," *Phys. Rev. A* **57**, 392 (1998).
   doi:10.1103/PhysRevA.57.392 (`uiterwaal_1998`).
3. A. W. McCown, M. N. Ediger and J. G. Eden, "Two-photon ionization of xenon at 193
   nm," *Phys. Rev. A* **26**, 3318 (1982). doi:10.1103/PhysRevA.26.3318
   (`mccown_1982`).
4. D. P. Shelton and J. E. Rice, "Measurements and calculations of the
   hyperpolarizabilities of atoms and small molecules in the gas phase," *Chem. Rev.*
   **94**, 3 (1994), eq. 3, Table I and eq. 14. doi:10.1021/cr00025a001
   (`shelton_1994`).
5. I. D. Petrov, B. M. Lagutin, V. L. Sukhorukov, A. Knie and A. Ehresmann,
   "Correlation and polarization effects in two-photon photoionization of Ar," *Phys.
   Rev. A* **93**, 033408 (2016). doi:10.1103/PhysRevA.93.033408 (`petrov_2016`).
6. A. Saenz and P. Lambropoulos, "Theoretical two-, three- and four-photon ionization
   cross sections of helium in the XUV range," *J. Phys. B* **32**, 5629 (1999).
   doi:10.1088/0953-4075/32/23/316 (`saenz_1999`).
7. A. Karamatskou and R. Santra, "Time-dependent configuration-interaction-singles
   calculation of the 5p-subshell two-photon ionization cross section in xenon," *Phys.
   Rev. A* **95**, 013415 (2017). doi:10.1103/PhysRevA.95.013415 (`karamatskou_2017`;
   the preprint was read).
8. B. M. Lagutin, I. D. Petrov, V. L. Sukhorukov, Ph. V. Demekhin, A. Knie and A.
   Ehresmann, "Relativistic, correlation, and polarization effects in two-photon
   photoionization of Xe," *Phys. Rev. A* **95**, 063414 (2017).
   doi:10.1103/PhysRevA.95.063414 (`lagutin_2017`).
9. A. L'Huillier, X. Tang and P. Lambropoulos, "Multiphoton ionization of rare gases
   using multichannel-quantum-defect theory," *Phys. Rev. A* **39**, 1112 (1989).
   doi:10.1103/PhysRevA.39.1112 (`lhuillier_1989`).
10. S. Kröll and W. K. Bischel, "Two-photon absorption and photoionization cross-section
    measurements in the 5p⁵6p configuration of xenon," *Phys. Rev. A* **41**, 1340
    (1990). doi:10.1103/PhysRevA.41.1340 (`kroll_1990`).
11. A. Mihelič and M. Horvat, "Calculation of multiphoton ionization amplitudes and
    cross sections of few-electron atoms," *Phys. Rev. A* **103**, 043108 (2021).
    doi:10.1103/PhysRevA.103.043108 (`mihelic_2021`).
12. L. A. A. Nikolopoulos and P. Lambropoulos, "Multichannel theory of two-photon single
    and double ionization of helium," *J. Phys. B* **34**, 545 (2001).
    doi:10.1088/0953-4075/34/4/304 (`nikolopoulos_2001`).
13. T. Sato et al., "Determination of the absolute two-photon ionization cross section
    of He by an XUV free electron laser," *J. Phys. B* **44**, 161001 (2011).
    doi:10.1088/0953-4075/44/16/161001 (`sato_2011`).
14. E. J. McGuire, "Two- and three-photon ionization in the noble gases," *Phys. Rev. A*
    **24**, 835 (1981). doi:10.1103/PhysRevA.24.835 (`mcguire_1981`).
15. A. F. Starace and T.-F. Jiang, "Transition-matrix theory for two-photon ionization of
    rare-gas atoms and isoelectronic ions with application to argon," *Phys. Rev. A*
    **36**, 1705 (1987). doi:10.1103/PhysRevA.36.1705 (`starace_1987`).
16. A. L'Huillier and G. Wendin, "Linear and nonlinear screening effects in two-photon
    ionization of xenon," *Phys. Rev. A* **36**, 4747 (1987).
    doi:10.1103/PhysRevA.36.4747 (`lhuillier_1987`).
17. C. Hättig, O. Christiansen and P. Jørgensen, "Coupled cluster response calculations
    of two-photon transition probability rate constants for helium, neon and argon,"
    *J. Chem. Phys.* **108**, 8355 (1998). doi:10.1063/1.476262 (`hattig_1998b`).
18. J. Henriksson, P. Norman and H. J. Aa. Jensen, "Two-photon absorption in the
    relativistic four-component Hartree-Fock approximation," *J. Chem. Phys.* **122**,
    114106 (2005). doi:10.1063/1.1869469 (`henriksson_2005`).
19. D. Shekhtman, M. A. Mustafa and N. J. Parziale, "Two-photon cross-section
    calculations for krypton in the 190-220 nm range," *Appl. Opt.* **59**, 10826
    (2020). doi:10.1364/AO.410806 (`shekhtman_2020`).
20. C. Drag, F. Marmuse and C. Blondel, "Measurement of the two-photon excitation
    cross-section of the 6p′[3/2]₂ and 6p′[1/2]₀ levels of Xe I at the wavelengths 224.3
    and 222.6 nm," *Plasma Sources Sci. Technol.* **30**, 075026 (2021).
    doi:10.1088/1361-6595/abfbeb (`drag_2021`; the accepted manuscript was read).
21. M. A. Dakka et al., "Laser-based metastable krypton generation," *Phys. Rev. Lett.*
    **121**, 093201 (2018). doi:10.1103/PhysRevLett.121.093201 (`dakka_2018`; the
    submitted manuscript was read).
22. D. B. Geohegan, A. W. McCown and J. G. Eden, "Resonantly enhanced three-photon
    ionization of krypton," *Phys. Rev. A* **33**, 269 (1986).
    doi:10.1103/PhysRevA.33.269 (`geohegan_1986`).
23. A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team, NIST Atomic Spectra Database
    (version 5.12). doi:10.18434/T4W30F (`nist_asd_2026`).

## Part II — Implementation record

### Phase 1 — The continuum model (2026-10-09)

- **Scope:** Phase 1, 2026-10-09, Opus 5.5 (AI) for John Travers, branch
  `feature/oep-0006-continuum`.
- **Delivered:**
  - `src/optmat/gas/two_photon.py`: `two_photon_ionisation_model`,
    `two_photon_ionisation_models`, `TwoPhotonIonisationModel` and its records
    (`TwoPhotonThreshold`, `TwoPhotonRegion`, `TwoPhotonReference`, `TwoPhotonSpread`,
    `TwoPhotonUncertainty`).
  - `src/optmat/gas/_two_photon_derivation.py`: the derivation.
  - `src/optmat/data/gas/two_photon.toml`: the selection.
  - `src/optmat/data/gas/models/two_photon/*.json`: five stored models (He, Ne, Ar
    linear; Xe linear and circular), written by `tools/fit/two_photon_record.py`.
  - `ArgumentError` and `UnavailableModelError` in `optmat.exceptions`.
  - Tests and goldens: `tests/test_two_photon.py`; the golden values and digests in
    `tests/data/golden_models.json`; the public-API snapshot.
  - Documentation:
    - the manual page {doc}`../../two_photon` (*The model*), with entries on
      {doc}`../../propagation` and {doc}`../../ionisation_rates`;
    - the gallery example *Two-photon ionisation of the noble gases*.
- **Deviations from the proposal:**
  - **Evaluation never raises on a frequency grid** (section 6 raised on an array with an
    unavailable point). It returns NaN there, with `available` as the mask, for scalars
    too: the author's decision after acceptance (TC-903, TC-943). It raises only for
    arguments outside their domain.
  - **`uncertainty` is a property, not a function of frequency:** the spreads are over
    the supported band, each with its overlap and number of points.
  - **Added:** `region` (the labels), `available` (the mask) and the registry
    `two_photon_ionisation_models()`.
- **Decisions:**
  - **The references** (Q1): approved by the author (TC-943). He: Mihelič and Horvat;
    Ne and Ar: McKenna and van der Hart; Xe: Lagutin et al., linear and circular.
  - **Krypton:** no model (Q2).
  - **The staircase rule:** vertices at one frequency become their geometric mean.
  - **The resonance wing:** 0.6 eV below $E_1$ for every gas, chosen where the
    independent calculations depart from the references (TC-944).
  - **Status:** Ne is `untested` (a single source); the others `validated`.
- **Register entries:** TC-943, TC-944. The data work before them is TC-872 to TC-942.
- **Verification:** run on the branch's tree:
  - `uv run poe check` (the full gate);
  - `PAPER_LIBRARY=… uv run poe test-library`;
  - `uv run poe test-release`.

  The tests reproduce:
  - every stored vertex from its dataset;
  - mckenna_2004's printed plateaus (Ar within 8 % of 1.2×10⁻⁵⁰ cm⁴ s from 8.0 to 11.0
    eV; Ne within 11 % of 1.65×10⁻⁵¹ from 11 to 16 eV);
  - the thresholds of section 4 to 0.01 nm;
  - $\beta$ against the energy-loss law and a fourth-order integration of it;
  - the region semantics on a grid;
  - the manual's table.

  No refit or fit change.
- **Sources checked:**
  - **Against their papers:** every reference's record was checked by its curation (TC-856,
    TC-874, TC-880), and the thresholds against `nist_asd_2026`'s printed values.
  - **Not independently re-audited in this phase:** the validators' records.
- **Limitations:**
  - **Gases:** no krypton model.
  - **Xenon** is unavailable between its two thresholds (Q3), where its only measurement
    (McCown et al.) lies.
  - **Neon** rests on one calculation.
  - **Circular polarisation** is xenon's only.
  - **Not modelled:** the bound two-photon resonances below $I_p/2$ (Phase 2), and a
    time-resolved loss for a broadband pulse (a documented gap, TC-903; Q9).
- **Open questions:** Q7–Q9 remain with Phase 2; Q8 with the optional Phase 3. Q1–Q6
  are settled or carried as above.
- **Handoff:**
  - **Phase 2** needs Q9's theory before any complex evaluator.
  - **The bound strengths** must first settle the convention of a J = 2 level's
    strength (one tensor component or the level, a factor 5; TC-942), and the line-shape
    factors of the cm⁴ s literature (TC-929).
  - **The paired solprop proposal** reads the model through an adapter, records its
    digest, and splits the bands itself.

### Phase 2a — The bound levels, curation only (2026-10-09)

- **Scope:** Phase 2 restricted to curation by the author (TC-945); 2026-10-09, Opus 5.5
  (AI) for John Travers, branch `feature/oep-0006-bound-levels`.
- **Delivered:**
  - `src/optmat/gas/two_photon_bound.py`: `two_photon_bound_levels`,
    `TwoPhotonBoundLevels`, `BoundLevel`, `BoundStrength`, `BoundPhotoionisation`,
    `sigma0_from_delta`.
  - `src/optmat/gas/_two_photon_bound_derivation.py`: the derivation.
  - `src/optmat/data/gas/two_photon_bound.toml`: the selection, every assignment with
    its reason.
  - Five stored records, written by `tools/fit/two_photon_bound_record.py`.
  - Tests and goldens: `tests/test_two_photon_bound.py`; golden values and digests.
  - Documentation: the manual's section *Bound two-photon excitation and the excited
    levels' ionisation*.
- **Deviations:**
  - **No complex evaluator** (Phase 2b, which waits on Q9's theory; TC-945), as OEP 6's
    own fallback provides.
  - **The record is a module of its own,** {mod}`optmat.gas.two_photon_bound`, beside the
    continuum model.
- **Decisions:**
  - **The convention** (TC-952): $\delta^{TP}$ per level, summed over $M_J$, parallel
    linear, equal photons, a.u.; $\sigma_0^{(2)} = 2\pi^3\alpha^2a_0^4\omega_\mathrm{au}^2
    \delta^{TP}$, normalised in the two-photon angular frequency, coherent light. A
    single tensor component of a $J = 2$ level is a fifth of it.
  - **Measured values** are held corrected for their stated photon statistics.
  - **Placement:** each level at NIST's energy.
  - **The rate constant K:** $\sigma_0 = K/4t_0$ for one beam of equal photons.
  - **TC-929 decided:** Pindzola's coefficient converts as $\alpha(\hbar\omega)^2$.
  - **Assessments:** He and Xe checked; Ne, Ar and Kr unresolved.
- **Register entries:** TC-952, TC-953 (the round's decisions TC-945).
- **Verification:**
  - **Gates:** the full gate, `poe test-library` and `poe test-release`.
  - **Xenon:** Henriksson et al.'s strengths reproduce Kröll and Bischel's and Drag et
    al.'s measurements to 0.83–1.30. Drag et al.'s $J = 0/J = 2$ ratio (1.38) matches
    per level (1.41), not per component (7.0).
  - **Helium:** the coupled-cluster 1s2s strength matches Drake's exact two-photon decay
    rate and two Rabi couplings within 1–4 %.
- **Sources checked:** each strength's defining equation against its paper (TC-952).
  Not re-audited: the photoionisation records beyond their own curation entries.
- **Limitations:**
  - **Krypton:** its measurements and calculations disagree by up to 94 times.
  - **Argon:** its only measurement is ambiguous by about 4.5.
  - **Neon:** calculations only.
  - **Strengths:** most levels have positions only; the higher series members have no
    curated strength.
  - **Theory placement:** the TDDHF assignments in xenon's 6p′/7p cluster are the
    curators'.
- **Open questions:**
  - Q7: Dakka et al.'s rate factor remains open, held as a check.
  - Q9: the evaluator's theory.
  - Krypton's and argon's measurements: to the roadmap's papers wanted.
- **Handoff:** Phase 2b needs Q9's derivation of the nondegenerate, channel-resolved
  amplitudes and of the population dynamics. It then reads these records: the level
  positions, the per-level strengths and the excited levels' photoionisation.

### Closure — Phases 2b and 3 transferred to OEP 8 (2026-10-12)

- **Scope:** the author closes OEP 6 as Final; 2026-10-12, Opus 5.5 (AI) for John
  Travers.
- **Delivered:** Phase 1 (the continuum model, TC-943, TC-944; optmat 0.5.0) and Phase 2a
  (the bound-level records, TC-952, TC-953), as reported above.
- **Transferred:**
  - **Phase 2b (the complex evaluator) and the optional Phase 3 (the real part against
    OEP 5's γ)** pass to [OEP 8](../oep-0008-bound-two-photon-response.md), as Q9's
    fallback provides.
  - **OEP 8 answers Q9 for the bound levels** and takes up Q7 and Q8 there.
- **Evidence gathered after Phase 2a (TC-985 to TC-1052):**
  - the literature of OEP 8, in five rounds;
  - five quantities for radiative and collisional data, and two more for pressure-line
    coefficients and bound two-photon amplitudes;
  - the maintainer's tRecX runs. haCC's helium cross section is 0.90–1.16 of the Phase 1
    model. The two-colour runs support the pair-frequency rule for the helium and argon
    continua (TC-1049).
- **Left open, on the roadmap:**
  - Q2, krypton (no tabulated or vector-drawn source);
  - Q3, xenon between its thresholds;
  - the continuum's time-resolved loss (TC-903), for a new OEP (OEP 8, Q6);
  - Q7, Dakka et al.'s rate factor.
- **Status:** Final.

## Part III — History

- 2026-10-08: Drafted.
- 2026-10-08: GPT-6 completed the AI review and revised Part I; findings and remaining
  acceptance question Q9 are recorded below.

- 2026-10-08: John Travers confirmed his completed human review and accepted the OEP
  in the review discussion. Review status is Reviewed. Q1 remains with Phase 1; Q7 and
  Q9 remain with Phase 2, with Q9 required before implementing the complex evaluator;
  Q8 remains with optional Phase 3. Part I is frozen on acceptance.

- 2026-10-09: The author approved the Phase 1 references (Q1) and the evaluation on
  grids (TC-943). Phase 1 implemented (TC-944); status In progress.

- 2026-10-09: The author restricted Phase 2 to curation (2a); the complex evaluator
  (2b) waits on Q9's derivation, OEP 6's own fallback. Xenon between its thresholds
  stays unavailable (Q3). The five further krypton leads were curated; none holds
  krypton data, so Q2 stands (TC-945 to TC-951). Released in optmat 0.5.0.

- 2026-10-09: Phase 2a (curation of the bound levels) implemented (TC-952, TC-953).
- 2026-10-09: The complex evaluator (Phase 2b) and the optional Phase 3 proposed in
  draft [OEP 8](../oep-0008-bound-two-photon-response.md), OEP 6's fallback for Q9.

### AI review, 2026-10-08

GPT-6 reviewed all of Part I against OEP 0. The review is complete; maintainer acceptance
and rule-8 reference approvals are separate decisions. Revisions address:

- **Broadband response:** retain the conditional integrated yield, correct the continuum
  spectral weight to $\sigma^{(2)}/\Omega^2$, remove the identification of susceptibility
  with a spectral factor, and require a separate derivation for time dependence (Q9).
- **Physical interpretation:** distinguish pulse bandwidth from level decay, bound
  excitation from (2+1) ionisation, and ionic-channel energy from electron kinetic energy;
  state the dilute-gas and perturbative limits and the missing mixed-frequency response.
- **Coverage and evidence:** remove the direct Lagutin–McCown validation outside shared
  coverage; preserve Sato's experimental limitations and source uncertainty; specify
  unavailable intervals, partial-channel sums, pole masks and evaluator identity.
- **Scope:** reuse the existing manual page, remove the unmerged-branch reference, and
  make the complex evaluator's remaining theory an explicit acceptance question.

**Primary-source checks:** inspected Shelton and Rice [4], pp. 3–4 (Taylor-series
normalisation and Table I); McKenna and van der Hart [1], Eq. 3 and the Ne/Ar
cross-section discussion; Lagutin et al. [8], p. 7 and Fig. 4 (coverage and the authors'
extrapolation statement); McCown et al. [3], p. 3322, Table I and Eqs. 9–10; Sato et al.
[13], p. 3, Table 2 and Eqs. 3–5. Thresholds were recomputed from the curated NIST
records and checked against the archived NIST responses. The other cited papers were
not independently re-audited in this review; their existing source records and register
entries remain the evidence. No reference curve was promoted, source value changed,
nondegenerate amplitude recovered, or Phase 2 model validated by this review.

**Verification of this revision:** the proposal-only gate (`git diff --check`, codespell
on the changed OEP and roadmap, proposal/register tests, and `poe docs`); additionally
the NIST-response and two-photon-page tests (90 passed, one expected process-OEP skip
across the four files). Sphinx ran all 53 gallery examples successfully. Recomputing the
illustrative argon coefficient gives 4.627×10⁻¹⁵ m W⁻¹, or 463 m⁻¹ at the assumed
intensity. The Gaussian Parseval check agrees to 3.6×10⁻¹⁶ relative and the monochromatic
flux and Im γ rates agree to floating-point precision. These are normalisation checks,
not a validation of the proposed broadband physical approximation. The full package,
extraction, release and refit gates were not run: this revision changes only the OEP
and its roadmap summary.
- 2026-10-12: Final, by the author. Phase 2b and Phase 3 transferred to OEP 8. Archived.

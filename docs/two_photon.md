# Two-photon ionisation and absorption

The two-photon processes of the noble gases He, Ne, Ar, Kr and Xe:
- two-photon ionisation from the ground state, above half the ionisation energy;
- two-photon excitation of bound even-parity levels below it;
- the photoionisation of the levels so excited.

The first has a model, {mod}`optmat.gas.two_photon` (OEP 6, Phase 1; *The model*,
below). The companion chapters describe {doc}`user-guide/two-photon-data` and
{doc}`user-guide/two-photon-excitation`. The reference curves are not fitted. Bound two-photon excitation has data but no model yet
(OEP 6, Phase 2).

The tables are written by `tools/two_photon_tables.py --update-docs` from the records,
and `tests/test_two_photon_page.py` checks that the page holds its output.

## The model

```python
import numpy as np

from optmat.constants import BAR, K_B, C
from optmat.gas.two_photon import two_photon_ionisation_model

argon = two_photon_ionisation_model("Ar")  # "He", "Ne", "Ar", "Xe"; Xe also "circular"
wavelength = np.array([160e-9, 150e-9, 130e-9, 100e-9])  # m, one photon
omega = 2 * np.pi * C / wavelength  # rad/s
argon.cross_section(omega)  # m^4 s: 0 below threshold, NaN where unavailable
argon.available(omega)  # where it has a value
argon.region(omega)  # below_threshold, supported, unsupported_detail, unavailable
argon.rate(omega, 1e17)  # 1/s per atom at 1e17 W/m^2
argon.loss_coefficient(omega, BAR / (K_B * 293.15))  # beta, m/W: dI/dz = -beta I^2
argon.pair_cross_section(2 * omega)  # indexed by a photon pair's summed frequency
argon.thresholds, argon.e1  # the 2P3/2 and 2P1/2 thresholds, and E1
argon.reference, argon.uncertainty, argon.identity  # what it is and how well known
```

This grid samples below-threshold, supported and unavailable regions. Evaluating
an `unsupported_detail` region also returns a value, with `ModelRangeWarning` as
described below.

Each model is a curated reference curve approved by the author under data rule 8
(TC-943). It is interpolated linearly in $\log\sigma^{(2)}$ between the curve's own
drawn vertices, and never fitted or extrapolated. Where a vector drawing's staircase puts
several vertices at one energy, they are one point at their geometric mean; partial
cross sections are summed on their common coverage (TC-944).

| Gas | Reference | Supported (nm) | Unsupported detail (nm) | Other calculations / reference | Status |
|---|---|---|---|---|---|
| He | Mihelič and Horvat 2021, ¹S + ¹D, velocity | 100.80–60.13 | 60.13–58.47 | 0.75–1.19 (Saenz; Nikolopoulos; van der Hart) | validated |
| Ne | McKenna and van der Hart 2004 | 114.48–77.15 | 114.80–114.48; 77.15–74.39 | none: single source | untested |
| Ar | McKenna and van der Hart 2004 | 155.59–112.47 | 156.71–155.59; 112.47–106.68 | 0.91–0.99 (Petrov et al., length and velocity) | validated |
| Kr | none (OEP 6, Q2) | — | — | — | — |
| Xe, linear | Lagutin et al. 2017, length | 184.47–158.21 | 158.21–146.96 | 0.65–0.87 (Fink and Johnson, raster) | validated |
| Xe, circular | Lagutin et al. 2017, length | 184.47–158.21 | 158.21–146.97 | 0.62–0.71 (Fink and Johnson, raster) | validated |

**The regions** tile the frequency axis:
- **below_threshold:** below $I_p/2$, the direct continuum cross section is exactly 0.
  Bound two-photon excitation and (2+1) ionisation there are not modelled (OEP 6, Phase
  2).
- **supported:** the reference's values.
- **unsupported_detail:** the reference's values, with
  {class}`~optmat.exceptions.ModelRangeWarning`. These are the 0.6 eV below $E_1$, where
  the first one-photon resonance's wing is placed differently by each calculation. For
  an LS-coupled reference (Ne, Ar) they are also the band between the ion's two
  thresholds.
- **unavailable:** NaN. This is between the threshold and the reference's first vertex
  (xenon: its two thresholds; Lagutin et al. compute nothing between them, OEP 6, Q3),
  and at and above $E_1$.

**Evaluation never raises on a frequency grid.** It returns NaN where the model is
unavailable and gives the mask (`available`); what to do there is the propagation
code's choice (the author's decision, TC-903, TC-943). It raises only for frequencies
that are not finite and positive, and for intensities or densities that are not finite
and non-negative. Asking for krypton, or for circular polarisation of a gas other than
xenon, raises {class}`~optmat.exceptions.UnavailableModelError` with the reason.

**The uncertainty** is kept in three kinds that are never merged (OEP 6, Q6):
- **the source's own accuracy:** Lagutin et al.'s 7.7 %; unstated for the others, and
  then held as None, not zero;
- **the digitising uncertainty at a drawn vertex:** 0.5 % for Mihelič and Horvat, 2 %
  for the others;
- **the spread of the other calculations over the supported band:** an empirical
  envelope, not a confidence interval. A gauge variant of the reference is a diagnostic
  and does not enter it.

**Composition.**
- **With the strong-field rates:** the model is the degenerate two-photon continuum. The
  strong-field rates ({doc}`ionisation_rates`) describe a pump at their carrier, and
  applying both to the same spectral component counts it twice. A band split is the
  propagation code's.
- **Mixed pairs:** a pump photon and a vacuum-ultraviolet one need nondegenerate
  amplitudes that no curated source gives (OEP 6, Q9).
- **With OEP 5's support record:** the lower threshold $I_p/2$ lies above OEP 5's band
  edge $\min(E_1, E_2/2)$. Between them lie the bound two-photon levels.

The gallery example *Two-photon ionisation of the noble gases* plots the models and the
loss of a prescribed field.

## Conventions

**The generalised two-photon ionisation cross section** $\sigma^{(2)}$ is defined by
the rate per atom,

$$
W = \sigma^{(2)}\Phi^2,\qquad \Phi = I/\hbar\omega ,
$$

with $\Phi$ the photon flux in cm⁻² s⁻¹, so $\sigma^{(2)}$ is in cm⁴ s
(1 cm⁴ s = 10⁻⁸ m⁴ s). This is the convention of Uiterwaal et al. and of McKenna and
van der Hart (TC-856), and of every calculation curated here, each record stating its
own definition. No factor 2! enters it.
- **Exception, McCown et al.:** they print a coefficient $\alpha$ in cm⁴ W⁻¹ defined by
  $dn/dt = \alpha[\mathrm{Xe}]I^2/\hbar\omega$ (their Eq. 9), so
  $\sigma^{(2)} = \alpha\hbar\omega$ (TC-884).
- **Intensity dependence:** a non-perturbative calculation (R-matrix Floquet) or a
  measurement gives $\Gamma/\Phi^2$ at a stated intensity. Off resonance it equals the
  lowest-order cross section.

Each event removes two photons. A quasi-monochromatic wave in a gas of number density
$N$ therefore loses intensity as

$$
\frac{dI}{dz} = -2\hbar\omega N W = -\beta I^2,\qquad \beta = \frac{2N\sigma^{(2)}}{\hbar\omega}.
$$

Of the $2\hbar\omega$ absorbed, $I_p$ ionises the atom and the rest is the electron's
kinetic energy.

**Bound two-photon excitation** is curated in two conventions (TC-872):
- **`two_photon_transition_strength`:** the rotationally averaged strength
  $\delta^{TP}$ in atomic units (Hättig et al.; Henriksson et al.), with the rotational
  average's convention and the polarisation stated.
- **`two_photon_excitation_cross_section`:** the frequency-integrated cross section
  $\sigma_0^{(2)}$ in cm⁴ or m⁴. The rate is $W = \sigma_0^{(2)} g\,\Phi^2$, with $g$
  the line profile normalised in angular frequency (every curated source: TC-889 to
  TC-892), and with the photon statistics $G^{(2)}(0)$ of the measurement stated.

A cross section in cm⁴ s at a line's peak depends on the line's width, and is not
curated as data. The ionisation of the excited level is
`excited_state_photoionization_cross_section`, in cm².

**Energy axes.** Every curve is held against the photon energy as its source drew it,
and no axis is shifted:
- **LS-coupled calculations** (McKenna and van der Hart; Petrov et al.) have one
  threshold, with resonances at their own computed positions.
- **Saenz and Lambropoulos** scaled their photon-energy axis so that the energy above
  threshold matches experiment (their ground state is 0.69 eV off; TC-873).
- **Karamatskou and Santra** shifted their Fig. 3 by 0.33 eV (TC-881).
- **L'Huillier, Tang and Lambropoulos** plot the total two-photon energy, held as the
  photon wavenumber, half of it (TC-882).

Near a resonance, curves are compared only at their own resonance positions, never at
NIST's.

## Thresholds

Direct two-photon ionisation of the ground state needs $2\hbar\omega \ge I_p$. For Ne,
Ar, Kr and Xe the ion has two fine-structure levels, ²P₃/₂ and ²P₁/₂. Between the two
thresholds, the Rydberg series converging on ²P₁/₂ autoionise.

Above the first one-photon resonance $E_1$ a single photon is resonantly absorbed; the
window of two-photon ionisation without one-photon resonance is therefore
$[I_p/2, E_1)$. Below $I_p/2$, down to OEP 5's $E_2/2$ (the lowest two-photon level),
lie the bound two-photon resonances.

From `nist_asd_2026`, the ion's ²P₁/₂ level added to the ionisation energy for the upper
threshold; photon energy in eV, vacuum wavelength in nm:

<!-- two-photon-thresholds:start -->
| Gas | I_p/2, ²P₃/₂ (eV (nm)) | ²P₁/₂ threshold (eV (nm)) | E₁ (eV (nm)) | E₂/2 (eV (nm)) |
|---|---|---|---|---|
| He | 12.2937 (100.85) | — | 21.2180 (58.43) | 10.3079 (120.28) |
| Ne | 10.7823 (114.99) | 10.8307 (114.48) | 16.6708 (74.37) | 9.2879 (133.49) |
| Ar | 7.8798 (157.34) | 7.9686 (155.59) | 11.6236 (106.67) | 6.5474 (189.36) |
| Kr | 6.9998 (177.13) | 7.3327 (169.08) | 10.0324 (123.58) | 5.7223 (216.67) |
| Xe | 6.0649 (204.43) | 6.7181 (184.55) | 8.4365 (146.96) | 4.8428 (256.02) |
<!-- two-photon-thresholds:end -->

Some consequences:
- Argon at the F₂ laser lines (157.52 and 157.63 nm; TC-900) is just below its
  threshold.
- Krypton at 193 nm and xenon at 248 nm are below theirs. Their ionisation there is
  two-photon resonant (2+1) ionisation (Geohegan et al., Kröll and Bischel).
- Xenon at 193 nm, McCown et al.'s measurement, lies between its two thresholds, in the
  autoionising region.

<span id="the-curated-two-photon-ionisation-cross-sections"></span>
<span id="helium"></span>
<span id="neon"></span>
<span id="argon"></span>
<span id="krypton"></span>
<span id="xenon"></span>

See {doc}`/user-guide/two-photon-data`.



<span id="bound-two-photon-excitation-and-the-excited-levels-ionisation"></span>

See {doc}`/user-guide/two-photon-excitation`.

## What optmat does not supply

**A time-resolved two-photon loss.** The cross section, and the integrated yield and
energy of a broadband pulse in the pair-spectral approximation, are OEP 6's (section 4).
How that loss is distributed in time within a pulse is not determined by
$\sigma^{(2)}$. It needs the nondegenerate two-photon amplitudes (OEP 6, Q9), which no
curated source supplies. It is a documented gap, for a new OEP (OEP 8, Q6) or for a
propagation code's own model (TC-903).

The continuum model and the bound-level records are what
[OEP 6](oeps/archive/oep-0006-two-photon-ionisation.md) proposed and implemented: its
theory, decisions and phase reports are the record. The bound resonances' response for
propagation is [OEP 8](oeps/oep-0008-bound-two-photon-response.md)'s.

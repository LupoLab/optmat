# OEP 3 — The molecular gases between their infrared bands

| | |
|---|---|
| **OEP** | 3 |
| **Title** | The molecular gases between their infrared bands |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | GPT-6 (AI); John Travers (human) |
| **Type** | Feature |
| **Created** | 2026-09-30 |
| **Resolution** | 2026-09-30: accepted by the author, John Travers, after the human and AI reviews; merged directly to `main` |
| **Requires** | [OEP 1](oep-0001-propagation-interface.md): the qualified intervals on the public surface, and identities |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phases 1, 2 and 3 and the open questions merged from `feature/oep-0003-windows-data`, `feature/oep-0003-intervals`, `feature/oep-0003-refits` and `feature/oep-0003-final` |
| **Manual pages** | {doc}`../../data` (*Infrared bands (HITRAN)*); {doc}`../../gas_models` (*The molecules' qualified intervals*); {doc}`../../models` (the air models' intervals) |

## Part I — Proposal

### 1. Abstract

The fitted models of CO, CO2, CH4 and N2O are qualified today from the ultraviolet to
1.4 µm (CO2, CH4) or 2 µm (CO, N2O). Their vibrational bands, from HITRAN, enter as one
pole per band. This OEP extends those fits through the near and mid infrared, into
every **window** between the bands where a pole model can be qualified. A pulse pumped at
2 µm and broadened to 4 µm in CO, or in CO2 outside its 2.7 µm band, is then described
by a qualified model.

It has four parts:
- **A stated target.** Between the bands the model is the gas's **background**
  refractivity at a stated spectral resolution. It is not the resolved structure of the
  weak lines that lie in every window.
- **Several qualified intervals per model.** A model can declare several intervals, and
  warns between them with the reason for each gap.
- **A band representation accurate at the window edges.** One pole per band is replaced
  by a Gaussian quadrature rule of a few nodes per band, built from HITRAN's lines. The
  pure rotational band of a polar molecule enters as its own band system.
- **Refitted parents.** The parent fits keep every band up to the HITRAN tables' limit.
  In each window they state a **representation bound**, measured against the complete,
  smoothed line sum, beside the fit's statistical uncertainty.

Success is judged in each window against a stated tolerance: the representation bound
on $n - 1$ and on the group-velocity dispersion, and the curated mid-infrared
measurements.

The feasibility calculation (§8) shows where this reaches.
- **CO:** its windows qualify almost entirely.
- **CO2:** most of 2.2-2.6 and 2.9-4.2 µm qualify, and 5.5-6.9 µm.
- **N2O and CH4:** their mid-infrared windows mostly do not. Weak bands inside them
  dominate the local dispersion.
- **The CO2-laser lines in pure CO2:** they lie inside CO2's own 9.4 and 10.4 µm bands,
  and are not qualified. In mixtures with a little CO2 they can be (§13 Q13).

SF6 is deferred to its own treatment (§13 Q14). OEP 2 (the compact models) then follows
the extended parents.

### 2. Motivation and scope

**The use.**
- **Propagation in gas-filled fibres and cells beyond 2 µm.** Mid-infrared pumps and
  supercontinua span several µm, and need $n(\lambda)$ and its dispersion across the
  whole spectrum. Where a molecule absorbs, or where its model is not accurate, the
  propagation code needs to know that it is outside the model.
- **The CO2-laser lines at 9-11 µm** in gases and mixtures that contain CO2, air among
  them (§13 Q13). In pure CO2 these lines lie inside CO2's own laser bands, and are out
  of this OEP's reach (§8).
- **A publication** of the gas models that covers the infrared windows (OEP 2's printed
  forms), which requires the windows to be qualified first.

**What limits the models today.**
- **The declared ranges end at 1.4-2 µm.** For CO2 and CH4 they end at 1.4 µm, for CO
  and N2O at 2 µm; only SF6 is declared to 5 µm. The atoms and the infrared-inactive
  molecules (N2, O2, H2, D2) are declared to the static limit.
- **The bands are cut.** Each fit keeps its HITRAN bands only up to its
  `infrared_max_wavenumber`, 4500-7000 cm⁻¹, although the tables reach 11000 cm⁻¹.
- **One pole per band is not accurate at mid-infrared window frequencies.**
  - Compared with the sum over HITRAN's lines outside a window, its error in the
    group-velocity dispersion there is 0.4 to 20 times the molecule's electronic GVD
    (medians over the window; §8).
  - Its error in $n - 1$ reaches $3\times10^{-5}$ to $7\times10^{-2}$ of $n - 1$ away from
    the window edges.
  - It was built for the optical side, where its error is $6\times10^{-6}$ of a band's
    contribution (TC-723).
- **Every window contains weak lines.** An undamped line sum diverges at each of them,
  so the resolved refractivity has narrow features of any size near every weak line. A
  pole model cannot follow them, and a propagation code for broadband pulses does not
  need it to (§4). CO2's 9.4 µm and 10.4 µm bands, whose lines lie in the 8-12 µm
  region, are the largest such bands.

**Where the windows might be.** The molecules are transparent between their bands. A
screening by each band's absorption averaged over its line span, below 0.1 m⁻¹ at 1 atm
and 296 K, gives these candidates in 0.8-12.5 µm:

| Gas | Candidate windows by the band-average screening (µm) |
|---|---|
| CO | 0.8-2.29, 2.52-4.34, 5.39-12.5 |
| CO2 | 0.8-1.94, 2.11-2.65, 2.93-4.17, 4.65-12.2 |
| CH4 | 0.8-1.31, 1.34-1.63, 1.70-2.08, 2.73-2.89, 4.77-5.41, 11.8-12.5 |
| N2O | 0.8-1.98, 2.00-2.10, 2.21-2.25, 2.34-2.56, 2.74-2.81, 3.08-3.51, 4.85-5.17, 5.59-7.34, 9.23-12.5 |

This screening is not the qualification.
- It spreads a band's strength over its whole span, and treats each band alone.
- The qualified windows follow from §4's criteria, spectral absorption and the
  representation bound, and §8 shows they are narrower.
- CO2's 4.70-5.18 µm absorbs above 0.1 m⁻¹ at the stated resolution although every band
  there passes the screening.

A 2-4 µm spectrum in CO crosses only the 2.3 µm overtone band. In CO2 it crosses the
2.7 µm band. In CH4 and N2O it crosses bands over most of that range. The OEP qualifies
the model where the criteria hold, and says plainly where it is not valid, and why.

**In scope:**
- CO, CO2, CH4 and N2O, whose HITRAN line lists give the lines;
- their pure rotational bands, which matter in the mid infrared (§13 Q12);
- qualified intervals for every model: a single interval for the models that do not
  change;
- the band quadrature rules, the windows and the representation bound, as derived data
  with their provenance;
- the refits of the molecular parents, their validators in the infrared and the manual.

**Out of scope:**
- **Inside the bands, and the resolved line structure in the windows.** The resonant,
  absorbing index: line shapes, line mixing, far wings, the complex index, and the
  narrow dispersion of each weak line in a window. The models are backgrounds at a stated
  resolution (§4), and warn inside the bands. A complex infrared index, which could also
  qualify pure CO2 at its laser lines, would be its own proposal.
- **SF6.** Its infrared input is a measured, broadened cross-section, which needs a
  different reference (§13 Q14). Its model and single interval are unchanged here.
- **Temperature.** HITRAN's intensities refer to 296 K. Hot bands, the CO2 laser bands
  among them, change strongly with temperature. The polarizabilities of CO2 and N2O
  refer to 293.15 K (TC-730). CO's and CH4's levels rest on 297.5 K data, and they
  declare no reference temperature. The difference from 296 K is stated, not modelled.
- **Pressure beyond the reference state.** The windows are defined at 1 atm and 296 K.
  - The absorption at the stated resolution scales with density while the collisional
    widths stay well below the resolution.
  - The peak absorption of a line does not: it is independent of pressure in the
    Lorentz regime.
  - Far-wing and collision-induced absorption grow with the square of the density, and
    are not included.

  The absorption map in the manual lets a user judge their own path and pressure (§9).
- **The compact models** (OEP 2), ozone (a Kramers-Kronig model with measured
  absorption), air and water vapour.

### 3. Context

**What optmat does today.**
- {mod}`optmat.gas.infrared` reads the band totals of
  `data/external/hitran/<gas>_bands.csv`. Per band these are the number of lines, the sum
  of the intensities $S$ and of $S/\tilde\nu^2$, the intensity-weighted wavenumber and
  the line span (`min_cm`, `max_cm`). Each band becomes one pole, keeping its strength
  and static contribution. Pure rotational bands are left out (TC-156).
- The tables are written by `tools/hitran/band_strengths.py` from the complete HITRAN2024
  line lists between 0 and 11000 cm⁻¹, fetched through the HITRAN API into a local cache
  and never redistributed. Only derived totals are committed
  (`data/external/hitran/PROVENANCE.md`).
- A band is a group of lines with the same isotopologue and vibrational labels. CH4's
  unassigned lines have an empty upper label, so they form 7 groups that span up to
  900-11000 cm⁻¹ and cross every window. They carry 1.3 % of CH4's infrared strength,
  in 132,832 lines.
- Each molecular fit configuration's `[model] infrared_max_wavenumber` keeps the bands
  below it as fixed trailing terms (`docs/fitting.md`). The values are CO 5000, CO2 7000,
  CH4 and N2O 4500, SF6 2000 (its cross-section windows) and O3 2500 cm⁻¹.
- `Model.wavelength_range` holds one interval, and `_check_range` warns
  (`ModelRangeWarning`) outside it.
- The declared ranges of the molecular fits:

  | Fit | Declared range |
  |---|---|
  | `fit/CO` | 167-2000 nm |
  | `fit/CO2` | 250-1400 nm |
  | `fit/CH4` | 240-1400 nm |
  | `fit/N2O` | 320-2000 nm |
  | `fit/SF6` | 230-5000 nm |

- CO2's configuration description and `docs/accuracy.md` also call it valid "between
  the bands, 8-12 µm", on the strength of Marchetti's 10 µm values. That window has
  never been part of its declared range.
- SF6's cached cross-sections show weak absorption at 1.5-5 µm, inside its declared
  range. Its peak cross-sections are $1.6\times10^{-22}$ to $1.6\times10^{-21}$ cm², or
  0.4-3.9 m⁻¹ at 1 atm. Its model keeps only the bands below 2000 cm⁻¹.

**Curated infrared data.** For these molecules, only CO2 has refractivity data beyond
2 µm:
- Marchetti & Simili's CO2-laser values at 9.2-10.7 µm (2007), which `fit/CO2` fits;
- their 2006 value at 10.57 µm, curated but neither fitted nor a validator today;
- Simmons's and Willetts & Harris's values at 10.6 µm, which are validators.

Old, Gentili & Peck's CO2 points at 1695 and 1817 nm are held out of the fit today
because one pole per band does not represent the edge of the 2 µm bands there (TC-148).
For SF6, Kumar & Meath's infrared oscillator strengths are curated.

**Register entries.**
- TC-104: HITRAN's strengths against Shelton's phase-matching ratios.
- TC-146 and TC-148: the molecular fits with their bands.
- TC-156: the pure rotational bands are left out; their static contribution is the
  orientation term.
- TC-720 to TC-723: Voronin & Zheltikov, Mathar, and Tang's bounds with their
  application to the per-band pole.
- TC-726: CO2's electronic form.
- TC-730: the reference temperature.

**Related OEPs.**
- **OEP 1**, the propagation-code interface (Draft), declares the public surface that
  propagation codes rely on. It declares `Model.intervals`, with one interval per model
  and `wavelength_range` as the envelope (its §4.1, Q30). It also covers:
  - the intervals and the representation bound in the identity (its §4.8);
  - the infrared terms' 296 K as metadata (its §4.9);
  - the infrared terms as one override input (its §6.3).

  This OEP supplies the several intervals of the molecules.
- **OEP 4**, the index of gas mixtures (Draft), and OEP 1's air model take a mixture's
  validity as the intersection of its constituents' intervals. This OEP states when a
  mixture may widen them (§13 Q13).
- **OEP 2**, the compact models (Draft), derives few-pole models from the parents. It
  is written to follow the parents' qualified intervals, and is sequenced after this
  OEP.

**Roadmap.** No informal item covers this. OEP 2's former Q3 (a band's second moment)
is taken up here (§4, §11).

### 4. Theory

**The line sum.** A line $\ell$ of HITRAN intensity $S_\ell$ (cm/molecule at 296 K,
including the lower-state population and stimulated emission) is an oscillator of
strength $f_\ell = S_\ell/(\pi r_e)$ at $\tilde\nu_\ell$. Its contribution to the
polarizability, in atomic units, is

$$
\alpha_\ell(\omega) = \frac{f_\ell}{E_\ell^2 - \omega^2}
= \frac{f_\ell}{2E_\ell}\left[\frac{1}{E_\ell - \omega} + \frac{1}{E_\ell + \omega}\right].
$$

The undamped sum diverges at every line. Near a weak line in a window its resolved
dispersion is narrow but large. On CO2's 4.65-12.2 µm candidate window the median GVD of
the resolved sum differs from the 8-node model by about $2\times10^{6}$ times the
electronic GVD (§8).

**The target: the background at a stated resolution.** This OEP's model is the gas's
refractivity smoothed by a normalized Gaussian of standard deviation $\sigma$ in
frequency. The resonant term of each line becomes

$$
\left[\frac{1}{E_\ell - \omega}\right]_\sigma
= \frac{\sqrt{2}}{\sigma}\, F\!\left(\frac{E_\ell - \omega}{\sqrt{2}\,\sigma}\right),
$$

with $F$ Dawson's function. The antiresonant term is smooth and is left as it is.
- **What it means.** Smoothing in frequency multiplies the response in time by a
  Gaussian gate of standard deviation $1/\sigma$. For $\sigma = 1$ cm⁻¹ that is 5.3 ps.
  The target is the gas's prompt response, and leaves out the long free-induction decay
  of the narrow lines. This is the dispersion a pulse sees whose spectrum spans many
  lines.
- **Far from the lines.** Where no line lies within several $\sigma$, the smoothed and
  undamped sums agree to a relative $(\sigma/\Delta)^2$, with $\Delta$ the detuning.
- **Its imaginary part** is the absorption at the same resolution,
  $a_\sigma(\tilde\nu) = N\sum_\ell S_\ell\, g_\sigma(\tilde\nu - \tilde\nu_\ell)$. This
  is HITRAN's definition with the line profile replaced by the Gaussian $g_\sigma$. The
  two are the real and imaginary parts of one smoothed susceptibility.
- **Collisional widths.** The reference neglects them. That is accurate while they are
  well below $\sigma$, below about 3 bar for self-broadened CO2 at 0.1 cm⁻¹ per atm.
  Above that, the physical broadening smooths further; Phase 1 assesses its effect.

**The reference.** The reference of this OEP is that smoothed sum over **every** line of
the line lists: the vibrational lines inside and outside each window, and the pure
rotational lines. It is evaluated in the windows. Its derivatives are $D_0 = \alpha$ and
$D_m = \partial^m(\omega\alpha)/\partial\omega^m$ for $m = 2, 3$, as in OEP 2. These
give the group-velocity dispersion ($m = 2$) and the third-order dispersion ($m = 3$) in
the dilute limit.

**One band as a measure.** Write $x = E^2$. A band $b$ is a positive measure with
weights $w_\ell = f_\ell/x_\ell$ at the nodes $x_\ell$. Its moments are
$\mu_k = \sum_\ell w_\ell x_\ell^{k}$, which are $M_{2k-2}$ in the notation of OEP 2:
$\mu_0 = M_{-2}$ is the static contribution and $\mu_1 = M_0$ the optical-side limit.
The band's contribution is the Stieltjes transform of that measure,

$$
\alpha_b(\omega) = \sum_\ell \frac{w_\ell\,x_\ell}{x_\ell - s}, \qquad s = \omega^2 .
$$

**The $n$-node rule.** The $n$-point Gaussian quadrature of this measure has positive
weights $W_j$ at $n$ nodes $X_j$ inside the band's line span. It reproduces
$\mu_0, \ldots, \mu_{2n-1}$, that is $M_{-2}, M_0, \ldots, M_{4n-4}$. Each node
becomes a pole of strength $W_j X_j$ at $E = X_j^{1/2}$:
- $n = 1$ is optmat's present pole per band (Tang's construction, TC-723);
- $n = 2$ is the exact two-point rule of OEP 2's former Q3, which needs $M_2$ and $M_4$;
- as $n$ grows, the rule converges to the line sum at every frequency outside the span.

The nodes and weights are the eigenvalues and first eigenvector components of the
$n\times n$ Jacobi matrix of the measure, computed from the lines by the Lanczos process
(Golub & Welsch [9], Eq. 2.6). A measure with fewer than $n$ distinct nodes gives its
exact rule with fewer nodes.

Gaussian quadrature of a Stieltjes transform is its Padé approximant from the same
moments; Tang & Poon [5] compare their moment bounds with such approximants. Where the
moments are known to all orders, the error falls rapidly with $n$ for frequencies well
outside the span, and slowly at the span's edge. §8 measures both.

**Band systems.** A measure is a group of lines with the same isotopologue and labels,
split at every window, so that no measure's span crosses a window. CH4's unassigned
lines are split the same way.

**Weak bands inside a window.** A band whose lines lie inside a window is treated in one
of two ways.
- **Folded.** Its measure moves to two poles at the nearest rule nodes below and above
  the window. Their strengths follow from the band's $M_0$ and $M_{-2}$, exactly and
  positively, because the band's $\bar\omega = (M_0/M_{-2})^{1/2}$ lies between them. The
  fold keeps the molecule's $M_0$ and $M_{-2}$. What it misplaces is part of the
  representation error, measured against the complete reference.
- **Kept.** It becomes its own band system with its own rule. Its poles lie inside its
  span, so the span, with its trimmed edges, is excluded, and the window splits there.

Phase 1 chooses per band, for the widest qualified span at the tolerance. A band may be
folded only where the fold's error is within the tolerance. Otherwise it is kept or its
neighbourhood is excluded. §8 shows that the choice decides CO2's mid-infrared windows.

**The rotational band.** A polar molecule's pure rotational lines lie below every
window: CO's below 300 cm⁻¹, N2O's below 85 cm⁻¹. At window frequencies they contribute
$-\sum_\ell f_\ell/\omega^2$ to a relative $(E_\mathrm{rot}/\omega)^2$. That is
$-77$ ppm of CO's electronic polarizability at 10 µm, and $-34$ % of its electronic GVD.
- They enter as one more band system, with its own rule, far below every window, and
  marked as rotational.
- **The static convention is unchanged.** The rotational system's static limit is the
  orientation term $\mu^2/3kT$ (TC-156), not a polarizability. The static
  polarizability, the static validators and the electric validators leave the
  rotational system out. Every refractive evaluation in a qualified interval includes it.
- Below the lowest interval the model warns, as it does anywhere outside its intervals.

**The representation error.** In a window the parent's infrared terms differ from the
reference by the signed, deterministic residual
$\delta_m(\omega) = D_m[\text{rules}] - D_m[\text{reference}]$. The parent stores a
**bound** $B_m(\omega) \ge |\delta_m(\omega)|$, established on a grid and checked
between its points (§7). It is not a standard deviation.
- It is stated beside the fit's statistical sigma, and not merged into it.
- A stated error interval adds the two linearly, $k\,\sigma_\mathrm{stat} + B_m$, since
  one is a bound and the other a standard deviation.

**Windows.** A window is an interval in which all of these hold:
1. **Transparency.** The smoothed total absorption $a_\sigma$, summed over every band,
   is below a threshold at the reference state (§7).
2. **Representation.** $B_0$ and $B_2$ are below the tolerance (§7) everywhere in it.
3. **Width.** It is at least a minimum width wide.

The first condition is physics. The second is the model's qualification, and it sets
how close to a band a window reaches. A gap between two windows records which condition
failed: absorption, representation or width. The windows are derived data, recomputed
when the line lists, the rules or the settings change.

**Qualified intervals.** A model's validity is a set of disjoint, ascending intervals of
vacuum wavelength. For the atoms and the infrared-inactive molecules it is one
interval, as now. For the molecules it is the optical interval, split where a
near-infrared band absorbs or cannot be represented, followed by the infrared windows.

**Conventions.** SI at every interface, vacuum wavelengths, the refractive
polarizability of the conventions page. The density conversion is unchanged. The line
list and the band rules refer to 296 K. The resolution $\sigma$ is part of the model's
stated meaning in the windows.

### 5. Primary sources and data

- **HITRAN2024** [6]. The line lists of CO, CO2, CH4 and N2O between 0 and 11000 cm⁻¹ are
  already fetched by `tools/hitran/band_strengths.py`. No new query is needed; a
  re-query under the same policy updates them. Only derived data are committed:
  - per band system, the $n$-node rule (nodes and weights), with the band's existing
    totals, and whether it is folded, kept or rotational;
  - per gas, the windows, the reasons for each gap and the representation bound on its
    grid;
  - each band's smoothed absorption at the reference state (Q13).
- **In hand, as validators (not fitted today):**
  - Old, Gentili & Peck's CO2 values at 1695 and 1817 nm [7], held out today for the
    reason this OEP removes (TC-148);
  - Simmons's [10] and Willetts & Harris's [11] 10.6 µm CO2 values (validators today).
- **In hand, fitted:** Marchetti & Simili's 2007 CO2-laser refractivity [8], which
  `fit/CO2` fits. Their 2006 value at 10.57 µm [12] is curated but not used in `fit/CO2`.
  - Both come from the same group, and probably the same apparatus. Phase 1 assesses
    their shared systematic errors before the 2006 value is used as a validator.
  - Neither is then independent validation of the 2007 fit.
  - Every one of these values lies in CO2's 9.4 and 10.4 µm bands, outside pure CO2's
    qualified intervals (§8). The representation error there is 1 % of their stated
    uncertainty, so the 2007 values stay fitted, with that error stated in their
    comparison.
- **To be sought** in Phase 1, by a literature search, and curated before they are used:
  - refractivity or dispersion measurements of CO, CO2, CH4 and N2O between 2 and
    12 µm;
  - static and dynamic vibrational polarizabilities from theory, as checks of $M_{-2}$
    and $M_0$.
  None are cited here from memory.
- **Methods:** Tang [4] and Tang & Poon [5] (curated); Mathar's line-by-line refractivity
  [2] (curated, method only); Golub & Welsch's quadrature algorithm [9]. The
  eigenvector-weight relation (Eq. 2.6) was checked against the paper in GPT-6's review.
  A source record is curated in Phase 1.
- **No new vocabulary quantity** is needed. The band rules and windows are derived
  tables with provenance, like the band totals.
- **Uncertainty.** In each window, two numbers are stated:
  - the fit's statistical sigma, unchanged in meaning;
  - the representation bound, measured against the complete smoothed reference.

  The limits of the line sum itself are stated as limitations, not quantified: line
  shapes, line mixing, far wings, missing weak lines and HITRAN's intensity
  uncertainties.

### 6. Proposed interface

- **Qualified intervals** (declared by OEP 1, §4.1).
  - `Model.intervals`: a tuple of `(low, high)` vacuum wavelengths in m, disjoint and
    ascending. `Model.wavelength_range` is its envelope,
    `(intervals[0][0], intervals[-1][1])`.
  - OEP 1 gives every model one interval, its present range. This OEP gives the
    molecular fits several.
  - A wavelength outside every interval warns (`ModelRangeWarning`). Between two
    intervals the message states the gap's recorded reason: an absorption band, a
    region the representation does not reach, or a gap narrower than the minimum width.
- **Fit configurations.**
  - `range_nm` remains the optical interval.
  - `[model] windows = true` adds the gas's derived windows. The derivation's
    resolution, thresholds, tolerance and treatment of each band are in a separate
    file, as OEP 2 does, so that the parent's digest changes only when the parent's
    inputs do.
  - `[model] band_nodes`, the largest number of nodes per band, replaces
    `infrared_max_wavenumber`. Every band of the tables is kept, and the rotational
    system of a polar molecule.
- **Stored fits.** `intervals_nm` beside `range_nm`, the gaps' reasons, and the
  representation bound on its grid, with the digest of the derived tables used.
- **The representation bound.**
  `FittedModel.representation_bound_au(omega, order)`, for `order` 0, 2 or 3.
  - It returns $B_m(\omega) \ge |D_m[\text{model}] - D_m[\text{reference}]|$ (§4), in
    atomic units: those of $\alpha$ for order 0, and those of $\alpha$ times
    $E_h^{1-m}$ for order $m$.
  - It is defined in every qualified interval, including the optical one, where it is
    computed like any other and not assumed to be zero.
  - Outside every qualified interval it returns NaN and warns (`ModelRangeWarning`).
  - The statistical methods, `polarizability_sigma_au` and `dispersion_sigma_relative`,
    keep their meaning.
- **Derived data.**
  - `data/external/hitran/<gas>_band_rules.csv`: per band system, its nodes (cm⁻¹) and
    weights, its treatment, with the band key of `<gas>_bands.csv`.
  - `data/external/hitran/<gas>_windows.csv`: the windows, each gap's reason, the
    criterion values at their edges, and the representation bound on its grid.
  - Each band's smoothed absorption at the reference state, beside its rule, so that a
    mixture's windows can be derived (Q13).
  - All are written by `tools/hitran/`, from the cached line lists, with the digests of
    those inputs (§8).
- **Example.**

  ```text
  from optmat.gas.models.registry import model
  co2 = model("fit/CO2")
  co2.intervals            # ((2.5e-07, 1.94e-06), (2.21e-06, 2.59e-06), ...)
  co2.refractivity(3.5e-6, state)            # in a qualified window: no warning
  co2.refractivity(2.7e-6, state)            # in the 2.7 um band: ModelRangeWarning
  co2.representation_bound_au(omega, 2)      # the bound on the GVD's representation error
  ```

  The intervals shown are illustrative; the derived windows decide them.

### 7. Implementation strategies

**Where.** `optmat.gas.infrared` (the rules, the band systems and the folding),
`optmat.gas.models.base` (`intervals`), `optmat.gas.fitting` and
`optmat.gas.models.fitted` (configuration and records), and `tools/hitran/`
(derivation). There are no core changes and no new dependencies; Dawson's function is
in SciPy.

**The representation. Recommended: an adaptive Gaussian rule per band system.**
- For each system, the smallest number of nodes $n \le$ `band_nodes` (recommended 8)
  that keeps its share of the error below a tenth of the tolerance at every point of
  the adjacent windows.
- Weak systems far from any window keep $n = 1$.
- The nodes lie inside the system's line span by construction. They never lie inside a
  window, since a kept in-window band splits the window.
- The rule reproduces the system's $M_{-2}$ and $M_0$ for every $n$. So the static and
  optical-side limits of the parent are those of the full line list, whatever $n$ is.
- Each in-window band is folded or kept (§4), chosen for the widest qualified span.

**The resolution. Recommended:** $\sigma = 1$ cm⁻¹ (a 5.3 ps gate). Its full width at
half maximum, 2.4 cm⁻¹, spans more than one line spacing of CO2 and N2O. Its collisional
limit is a few bar. §8 gives the results at
0.3 and 3 cm⁻¹ beside it.

**The tolerance. Recommended** (§13 Q2): in every window,
- the representation bound of $\alpha$ is below 0.1 of the fit's statistical sigma of
  $\alpha$. In the windows that is a relative $4\times10^{-6}$ for CO, $1\times10^{-6}$
  for CO2, $4\times10^{-6}$ for CH4 and $1.5\times10^{-5}$ for N2O;
- that of the GVD is below 1 % of the larger of the reference's $|D_2|$ and the
  electronic $|D_2|$ at the same wavelength. In the mid-infrared candidate windows the
  median total GVD is 45 to 100,000 times the electronic GVD. The electronic term is a
  floor where the total crosses zero.

This replaces the draft's 1 % of the electronic GVD. Once the reference includes every
line, the weak bands' real contributions alone exceed that bar by orders of magnitude,
although they are 1 % or less of the total GVD.

**The windows. Recommended** (§13 Q1, Q11):
- the transparency criterion is the smoothed total absorption $a_\sigma$ below 0.1 m⁻¹,
  for the pure gas at 101325 Pa and 296 K;
- each window's record also states the largest $a_\sigma$ in it and the largest
  single-line Lorentz peak, self-broadened at 1 atm, for users of narrowband light;
- the edges are then trimmed to where the representation tolerance holds;
- windows narrower than 50 cm⁻¹ are dropped;
- each gap's reason is recorded.

**The bound between grid points.** The representation bound is a maximum over a grid.
- The grid spacing is at most $\sigma/3$ and at most a tenth of the distance to the
  nearest model pole.
- The bound is stored per cell of the grid.
- An independent grid, offset by half a spacing and twice as dense, checks it: every
  sample must lie within its cell's bound, or the cell is refined and its bound
  raised.
- The smoothed reference is an entire function, and the model's poles lie outside the
  windows. So refining the grid converges.

**The fits.**
- Every band of the tables (to 11000 cm⁻¹) is kept, as rules or folded. The rotational
  system is added for CO and N2O. CH4's is about 3 % of its tolerance on $\alpha$, and
  CO2 has none.
- The refits take the new fixed terms in dependency order, as for any refit:
  - the molecules only;
  - their partners in Shelton's ratios are N2, which has no bands, and so unchanged.
- The fitted data are unchanged at first. The validators gain Old's 1695 and 1817 nm
  points and the infrared data found in Phase 1.
- Whether Old's two points are then fitted is a decision for the register, after the
  comparison.

**Cost.**
- The rules at 8 nodes for every band would give CO2 about 35,000 fixed poles, against
  4296 now; the adaptive rule needs far fewer, and is measured in Phase 1.
- The fixed terms are evaluated once per observation and cached, so the fits' cost grows
  little.
- Evaluating a model at a wavelength costs one sum over its poles.
- The line-by-line reference is computed by the tool, not in the tests.

**Option B, one pole per band with narrower windows.** Qualify only the parts of each
window where one pole is accurate. At the recommended tolerance that leaves almost
nothing in the mid infrared (§8). Rejected.

**Option C, the line sum at run time** (Mathar [2]). Exact, but hundreds of thousands of
lines per evaluation, and the lists may not be redistributed. Rejected for the model.
Kept as the reference.

**Option D, two poles per band from $M_2$ and $M_4$ only** (the former OEP 2 Q3 as
proposed). This is the $n = 2$ case of the recommended rule. It is not enough alone: in
the mid-infrared windows its GVD error is 0.2-44 % of the electronic GVD (medians, §8).

**Option E, the resolved refractivity as the target** ($\sigma \to 0$). No pole model
approaches it in any window with lines (§8). Qualification would need an excluded
neighbourhood around each of tens of thousands of weak lines. Rejected.

**Option F, damped oscillators for the weak bands.** A complex pole per band, or a
smoothed band shape, could follow a weak band's background inside its own span, and
would qualify pure CO2 at its laser lines. It is a complex-index model. Deferred to a
later proposal.

### 8. Verification and validation

- **Feasibility against the lines outside each window** (computed 2026-09-30 from the
  cached HITRAN2024 line lists and the stored fits' electronic terms). This measures how
  fast the rules converge.
  - Method: lines outside each window, grouped into bands (and split at the window),
    with each band replaced by its $n$-node rule. Compared with the undamped line sum
    over the same lines, on 1500 geometric points per window.
  - "GVD" is the error of the total $D_2$ in units of the electronic $|D_2|$. "Core"
    excludes 2 % of the window at each edge.

  | Gas, window | Nodes per band | $n - 1$, core max | GVD, median | GVD < 1 %: window (core) |
  |---|---|---|---|---|
  | CO, 2.52-4.34 µm | 1 | $2\times10^{-3}$ | 0.73 | 2 % (2 %) |
  | | 2 | $7\times10^{-5}$ | $6\times10^{-3}$ | 55 % (59 %) |
  | | 4 | $3\times10^{-7}$ | $8\times10^{-7}$ | 88 % (91 %) |
  | | 8 | $4\times10^{-11}$ | $2\times10^{-13}$ | 98 % (100 %) |
  | CO2, 2.11-2.65 µm | 1 | $1\times10^{-4}$ | 0.43 | 3 % (4 %) |
  | | 4 | $1\times10^{-9}$ | $2\times10^{-7}$ | 88 % (99 %) |
  | | 8 | $8\times10^{-15}$ | $3\times10^{-13}$ | 97 % (100 %) |
  | CO2, 2.93-4.17 µm | 1 | $7\times10^{-2}$ | 20 | 0 % (0 %) |
  | | 2 | $2\times10^{-3}$ | 0.10 | 25 % (28 %) |
  | | 4 | $7\times10^{-6}$ | $6\times10^{-6}$ | 80 % (85 %) |
  | | 8 | $3\times10^{-10}$ | $2\times10^{-12}$ | 94 % (100 %) |
  | CO2, 4.65-12.2 µm | 1 | $1\times10^{-3}$ | 9.4 | 0 % (0 %) |
  | | 4 | $2\times10^{-8}$ | $7\times10^{-7}$ | 91 % (95 %) |
  | | 8 | $1\times10^{-11}$ | $7\times10^{-13}$ | 98 % (100 %) |
  | N2O, 3.08-3.51 µm | 1 | $6\times10^{-5}$ | 3.6 | 0 % (0 %) |
  | | 4 | $5\times10^{-10}$ | $1\times10^{-4}$ | 77 % (100 %) |
  | | 8 | $1\times10^{-14}$ | $1\times10^{-9}$ | 93 % (100 %) |
  | CH4, 4.77-5.41 µm | 1 | $5\times10^{-5}$ | 2.5 | 0 % (0 %) |
  | | 4 | $6\times10^{-8}$ | $4\times10^{-2}$ | 34 % (49 %) |
  | | 8 | $4\times10^{-10}$ | $5\times10^{-4}$ | 71 % (98 %) |
  | CH4, 1.70-2.08 µm | 1 | $3\times10^{-5}$ | 0.42 | 2 % (3 %) |
  | | 4 | $4\times10^{-6}$ | $1\times10^{-2}$ | 51 % (56 %) |
  | | 8 | $3\times10^{-7}$ | $2\times10^{-4}$ | 75 % (84 %) |

  With 8 nodes per band, the rules converge to the outside lines over the cores of the
  windows. One pole per band, today's model, fails almost everywhere in the mid-infrared
  windows. This table leaves the in-window lines out of both sides, so it does not
  qualify any window.
- **Feasibility against the complete reference** (computed 2026-09-30, after GPT-6's
  review). The reference now includes every line, smoothed at $\sigma$ (§4).
  - **The resolved sum.** On CO2's 4.65-12.2 µm, with every in-window band folded and
    8 nodes per band, the median GVD error is:
    - $2.0\times10^{6}$ electronic GVDs at $\sigma = 0$, which reproduces the review;
    - $6.9\times10^{3}$ at 0.3 cm⁻¹;
    - 56 at 1 cm⁻¹;
    - 43 at 3 cm⁻¹.

    The median error in $n - 1$ is $8\times10^{-6}$ at every $\sigma$. Above 1 cm⁻¹ the
    error comes from the weak bands' smooth envelopes, not their lines.
  - **The qualified parts.** At $\sigma = 1$ cm⁻¹, with 8 nodes per band and the
    tolerance of §7: $\alpha$ within 0.1 statistical sigma, and the GVD within 1 % of
    the larger of the total and the electronic GVD. "Kept" bands are the in-window bands
    above $10^{-7}$ of the molecule's $M_0$; the others are folded. The rotational
    system is included for CO and N2O, as a 4-node rule.

    | Gas, candidate window | In-window bands | Treatment | Passing the tolerance |
    |---|---|---|---|
    | CO, 2.52-4.34 µm | 7 ($7\times10^{-10}$ of $M_0$) | all folded | 98 % (core 100 %) |
    | CO, 5.39-12.5 µm | 2 ($2\times10^{-13}$) | all folded | 100 % |
    | CO2, 2.11-2.65 µm | 404 | all folded | 60 % (core 71 %) |
    | | | 2 kept | 74 % of its 2.21-2.59 µm part |
    | CO2, 2.93-4.17 µm | 396 | all folded | 75 % (core 84 %) |
    | | | 5 kept | 98 % of its 3.21-3.73 µm part |
    | CO2, 4.65-12.2 µm | 681 ($7.6\times10^{-5}$) | all folded | 9 % |
    | | | 19 kept | all of 5.54-6.95 µm; none of 8.41-8.95 µm |
    | N2O, 3.08-3.51 µm | 223 | all folded | 10 % |
    | N2O, 9.23-12.5 µm | 80 | all folded | 0 % |
    | CH4, 4.77-5.41 µm | 11 | all folded | 0 % (core 1 %) |

    Keeping N2O's and CH4's in-window bands above $10^{-7}$ of $M_0$ leaves no part
    wider than 50 cm⁻¹.
  - **What this shows.**
    - CO's windows qualify almost entirely.
    - CO2's short-wave windows qualify over most of their span. Its 4.65-12.2 µm
      candidate qualifies over 5.54-6.95 µm.
    - In N2O's and CH4's mid-infrared candidates, weak bands dominate the local
      dispersion. There the median GVD error is 30-95 % of the total GVD. A pole model
      does not qualify them at this tolerance (Option F).
  - **CO2's laser bands.** The 9.4 and 10.4 µm bands span 854-1110 cm⁻¹ (9.0-11.7 µm).
    - Kept, their poles fill that span.
    - Folded, the error at 10.6 µm is $4.2\times10^{-5}$ of $n - 1$ (4.3 statistical
      sigma) and 54 % of the total GVD, or 18 % at $\sigma = 3$ cm⁻¹.
    - Pure CO2 is therefore not qualified at its laser lines.
    - The error is 1 % of Marchetti's stated uncertainty (§5).
    - At 400 µmol/mol in air the folded error scales with the mole fraction (§13 Q13).
  - **The rotational band.** At 10 µm CO's contributes $-77$ ppm of its electronic
    polarizability and $-34$ % of its electronic GVD. N2O's contributes $-26$ ppm and
    $-12$ %. That is 0.15 % of the total GVD in the 5-12 µm windows. On $\alpha$ it is
    20 times CO's tolerance and 1.7 times N2O's. With it included as a 4-node system,
    CO's 5.39-12.5 µm error stays below $10^{-12}$ in $n - 1$.
  - **Absorption** ($a_\sigma$, 1 atm, 296 K, pure gas).
    - CO2's 4.70-5.18 µm exceeds 0.1 m⁻¹, reaching 1.8 m⁻¹ at 4.81 µm. The largest
      single-line Lorentz peak there is 1.7 m⁻¹, at 4.859 µm; 135 lines exceed
      0.1 m⁻¹.
    - N2O reaches 0.11 m⁻¹ at 3.16 µm, at $\sigma = 1$ cm⁻¹ only, and 0.13 m⁻¹ at
      10.53 µm, over 10.43-10.85 µm.
    - CH4's 1.70-2.08 µm candidate exceeds 0.1 m⁻¹ over 1.70-1.80 µm, reaching
      1.15 m⁻¹.
    - The other candidates computed (CO's, CO2's 2.11-2.65 and 2.93-4.17 µm, CH4's
      4.77-5.41 µm) stay below 0.05 m⁻¹.
  - These calculations use 600-1500 geometric points per interval, not the grids of §7,
    and the draft's rule for keeping bands. Phase 1 recomputes every window with the
    retained tool and the §7 grids, and chooses fold or keep per band.
- **Invariants and synthetic tests** (no HITRAN access needed):
  - every rule has positive weights and nodes inside its system's span;
  - it reproduces $\mu_0, \ldots, \mu_{2n-1}$ of its measure to a relative $10^{-12}$;
  - on small synthetic measures, it reproduces the direct line sum at points outside the
    support;
  - the synthetic cases cover duplicate nodes, fewer distinct lines than $n$ (an exact
    rule with fewer nodes), a measure split at a window (the parts sum to the whole), and
    folding (the molecule's $M_0$ and $M_{-2}$ kept, both poles outside the window);
  - independent samples between the qualification grid's points lie within the stored
    bound;
  - the rules of a gas sum to the band table's totals;
  - no pole lies inside a qualified interval (an analytic check, independent of any
    grid);
  - the intervals are disjoint and ascending, and every gap carries a reason;
  - a model with one interval behaves exactly as before.
- **Provenance and freshness.** The derived tables carry:
  - the SHA-256 digest of each raw input, the cached line list and its query record;
  - the derivation's version and the tool's commit;
  - every setting: $\sigma$, the thresholds, the tolerance, the node cap, the grid
    rules and each band's treatment;
  - each window's exact endpoints and its grid's definition.

  The feasibility calculation is retained in `tools/hitran/`, and reruns from the cache.
  The stored fits carry the derived tables' digest, and the tests fail if either is
  stale. The tests cannot reach HITRAN, so they check the chain from the tables onwards.
- **The reference.** The representation bound is recomputed from the line lists by the
  tool and recorded in the register.
- **Validators:**
  - Old's CO2 values at 1695 and 1817 nm;
  - the 10.6 µm CO2 values of Simmons and of Willetts & Harris;
  - Marchetti's 2006 value, after the assessment of §5;
  - the infrared data found in Phase 1.

  Each is compared with the parent within its own uncertainty and the parent's
  statistical sigma plus representation bound. Marchetti's 2007 values are reproduction
  checks.
- **Unchanged:** every non-molecular fit, and SF6. A test compares the refitted
  molecules at 250-1400 nm with the present fits and reports the change. It is expected
  to be a small fraction of the statistical sigma, from the added weak bands. A change
  larger than one sigma is a finding for the register.
- **Parity:** solprop does not consume these ranges today. If its interface reads
  `wavelength_range`, the envelope semantics are to be checked with OEP 1.

### 9. Documentation

- `docs/data.md`, *Infrared bands (HITRAN)*: the band rules, the folding, the rotational
  system and the windows, with their provenance.
- `docs/fitting.md`: the representation, its bound, the resolution, and the windows in a
  fit.
- `docs/gas_models.md` and `docs/accuracy.md`:
  - each molecule's qualified intervals, and the reason for each gap;
  - an absorption map of 0.8-12.5 µm with the windows marked, and the peak single-line
    absorption, so a user can judge their own path length, pressure and bandwidth;
  - that the model is the background at a stated resolution, not the resolved structure
    of the weak lines;
  - the representation bound beside the statistical sigma, and how to combine them;
  - CO2's 8-12 µm wording replaced by the qualified windows, and the statement that
    pure CO2 is not qualified at its laser lines.
- `docs/conventions.md` or `docs/models.md`: qualified intervals, the warning between
  them, and the rotational system's place in the static convention.
- A gallery example: CO2's refractivity and GVD from 0.25 to 12 µm, with the bands and
  unqualified regions shaded and the representation bound.

### 10. Compatibility and impact

- **Models with one interval** (the atoms, N2, O2, H2, D2, SF6, the published formulas
  and the comparison entries) are unchanged.
- **The molecular fits are refitted.** Their optical-interval predictions change by the
  weak bands now kept and by the rotational systems (reported in Phase 3). Their stored
  fits and robustness records are regenerated, with register entries.
- **New warnings.** A molecular model evaluated in one of its absorption bands inside
  its old range now warns. CH4 near 1.31-1.34 µm is one case, if the criteria keep that
  band. That is intended.
- **`wavelength_range` is the envelope.** OEP 1 declares it so, and its value does not
  change until this OEP's refits. For the molecules, code that tests one wavelength
  against it would then accept wavelengths inside a band. Such code uses `intervals` or
  the warning, as OEP 1's page says.
- **Air and mixtures.** OEP 1's full-composition air model and OEP 4's mixtures with
  molecules inherit the constituents' gaps through the intersection of intervals, until
  the checks of §13 Q13 widen them.
- **Fixtures.** solprop's parity fixtures for CO, CO2, CH4 and N2O are re-recorded after
  Phase 3 (OEP 1, §6.11).
- **OEP 2** follows the parents: its compact models take the parents' qualified
  intervals, with their poles kept outside every one of them. Its weights use the
  parents' statistical sigma. The representation bound is stated beside it, not
  combined in quadrature. The rotational system is carried unchanged, outside the
  reduction.
- **Marchetti's 2007 values** stay fitted, at wavelengths outside pure CO2's qualified
  intervals, with the representation error stated in their comparison.
- No register entries are superseded. TC-148's exclusion of Old's two points is
  revisited in Phase 3. TC-156's static convention is kept; the rotational system
  extends it (§4).

### 11. Alternatives considered

- **Two poles per band from fitted wavelengths** (Voronin & Zheltikov [1]). Their pairs
  are fitted to the optical side of air, not to the lines at the window edges, and
  constrain only the band strength. Superseded here by rules computed from the lines.
- **Each band's $M_2$ alone** (the former OEP 2 Q3). It gives two-sided bounds, or a
  two-point rule with one node fixed at a band edge (a Radau rule). The full two-point
  rule needs $M_4$ as well. It is option D at best, and does not reach the tolerance in
  the mid infrared.
- **Sub-band binning.** Split each band into wavenumber bins, one pole per bin. It
  converges, but more slowly than a Gaussian rule with the same number of poles, and
  puts poles at bin centres rather than where the band's structure needs them.
- **Mathar's window polynomials** [3]. Valid inside a window only, without limits or
  physical meaning, and per mixture. Rejected, as in OEP 2.
- **A Lorentzian smoothing** (a complex frequency shift). It is exact for a causal
  function, but its kernel's tails move a strong band's absorption into its
  neighbouring windows as $1/\Delta^2$. Rejected for the transparency criterion. The
  Gaussian's tails fall off fast enough.
- **Dropping the rotational band.** Its omission is 20 times CO's tolerance on
  $\alpha$ in the 5-12 µm windows, and 1.7 times N2O's. Rejected.
- **Leaving the models at their present ranges.** This is the status quo, which fails
  the mid-infrared use.

### 12. Phases

1. **The representation and the windows (derived data, no model change).**
   - `tools/hitran/`:
     - the band systems, split at windows;
     - the adaptive rules and the rotational systems;
     - the fold-or-keep choice per in-window band;
     - the smoothed reference and absorption;
     - the windows with their gaps' reasons;
     - the representation bound with its independent check;
     - the retained feasibility calculation.
   - `<gas>_band_rules.csv` and `<gas>_windows.csv` with their provenance and raw-input
     digests.
   - The collisional-width assessment of $\sigma$, and the Lorentz-peak and far-wing
     absorption at the reference state.
   - The literature search for infrared refractivity and dispersion of the four
     molecules, and the curation of what it finds. The Golub & Welsch source record.
     The assessment of Marchetti's 2006 and 2007 values.
   - Register entries:
     - the method, the resolution and the window criteria;
     - the findings of §2 and §8: one pole per band in the windows, the resolved sum,
       the weak bands' envelopes, the rotational band, CO2's laser bands, and CH4's
       unassigned lines;
     - SF6's weak absorption at 1.5-5 µm inside its declared range;
     - each source curated.
   - Settle Q1, Q2, Q3 and Q6 with the recomputed windows.
   - Verify: the invariants and synthetic tests, and the gate.
2. **Several intervals per model** (after OEP 1's Phase 1, which declares `intervals`
   with one per model).
   - The warning's message with the gap's reason.
   - The configuration keys and stored-record fields; `representation_bound_au`, and
     its place in the identity.
   - The rotational system's exclusion from static and electric comparisons.
   - Every existing model keeps its single interval, behaviour unchanged.
   - Docs: the conventions or models page.
   - Verify: every stored model's warnings unchanged, and the gate.
3. **The refitted molecules.**
   - The molecular fits with every band as rules, the rotational systems, the windows,
     and the refits and robustness runs in dependency order.
   - The validators in the infrared, and the decision on Old's two points.
   - The mixture rule of Q13, with OEP 1's air models and OEP 4.
   - The manual pages of §9 and the gallery example.
   - Register entries for each fit decision and change.
   - Verify: the gate, `test-fits`, the optical-interval comparison with the present
     fits, and the infrared validators.

OEP 2 is implemented after Phase 3.

### 13. Open questions

1. **The window criterion.** The band-average absorption over the line span is crude.
   It spreads a band's strength over its whole span, while individual lines peak far
   above that average.
   - Alternatives: a bin-averaged absorption (for example over 10 cm⁻¹), or a
     line-shape absorption at a reference pressure.
   - Also open: the threshold, and the reference state (1 atm and 296 K, or the
     models' 293.15 K).
   - **Recommended:** the band average at 101325 Pa and 296 K, with a threshold of
     0.1 m⁻¹, trimmed by the representation tolerance. Record the criterion's value at
     each window edge, so a user with a longer path can see how close they are.
   - To be settled in Phase 1, when the absorption map is computed.
   - *Revised (2026-09-30) by Q11:* the criterion is the smoothed total absorption
     $a_\sigma$ below 0.1 m⁻¹ at 101325 Pa and 296 K. Each window's record states the
     largest $a_\sigma$ and the largest single-line Lorentz peak. The band average
     remains a screening. **Phase 1**, with the recomputed windows.
2. **The tolerance in the windows:** $10^{-6}$ in $n - 1$, and 1 % of the electronic
   GVD. Is the GVD bar right for propagation? Should $\beta_3$ be bounded too?
   - **Recommended:** as stated, with $\beta_3$ reported but not bounded. Near a band
     it is dominated by the band.
   - *Revised (2026-09-30):* 0.1 of the statistical sigma in $\alpha$, which is
     $1\times10^{-6}$ for CO2 and $4\times10^{-6}$ to $1.5\times10^{-5}$ for the others.
     For the GVD, 1 % of the larger of the total and the electronic $|D_2|$ (§7). Against
     the complete reference, the electronic bar fails on real weak-band contributions of
     1 % or less of the total GVD (§8). $\beta_3$ is reported, not bounded.
     **Phase 1**, confirmed with the recomputed windows.
3. **Weak bands inside a window: fold, drop, or keep.**
   - **Recommended:** fold them into the nearest band system outside, which keeps $M_0$
     and $M_{-2}$, and add their misplaced local contribution to the representation
     uncertainty.
   - Dropping them loses those moments. Keeping them keeps a spurious pole in the
     window.
   - *Revised (2026-09-30) by Q10:* fold or keep per band, for the widest qualified
     span. A fold goes to two poles, one on each side of the window, which keeps $M_0$
     and $M_{-2}$ exactly. A kept band splits the window at its span. The error is
     measured against the complete reference, not bounded by the band's own lines
     (§4, §8). **Phase 1.**
4. **How the representation uncertainty is exposed.**
   - **Recommended:** a separate method, with the statistical methods unchanged, so a
     user or OEP 2 can combine the two as they need. The manual states both.
   - The alternative is a combined sigma, which would change the meaning of the
     existing methods.
   - *Resolved (2026-09-30) by Q15:* a separate method,
     `representation_bound_au`, returning a bound, with a linear combination rule.
5. **`wavelength_range` and OEP 1.**
   - **Recommended:** `intervals` becomes the declared validity, and
     `wavelength_range` the envelope, documented as such.
   - The alternative is keeping `wavelength_range` as the principal (optical) interval.
     That is safer for code that reads it, but less honest in the infrared.
   - *Resolved (2026-09-30)* in OEP 1 (§4.1, Q30), as recommended. OEP 1 declares the
     tuple with one interval per model, so this OEP adds values, not a change of
     surface.
6. **Fixed or adaptive nodes.**
   - **Recommended:** adaptive, up to 8, recorded per band.
   - The alternative is a fixed $n$ for every band: simpler, but it multiplies CO2's
     poles by $n$.
   - *Disposition (2026-09-30):* as recommended. **Phase 1** measures the poles needed.
7. **SF6.** Its cross-sections give each band's measure directly, at 0.112 cm⁻¹
   resolution, and its range already reaches 5 µm.
   - **Recommended:** include it in Phase 1 with the same rules from the measured
     spectrum. Its windows beyond 5 µm are to be qualified if the rules pass.
   - *Revised (2026-09-30) by Q14:* SF6 is deferred to its own treatment. Its model is
     unchanged. **Phase 1** records the weak absorption found inside its declared range.
8. **Temperature.** The hot bands (CO2's laser bands among them) change strongly with
   temperature, and HITRAN refers to 296 K.
   - **Recommended:** out of scope here, stated as a limitation in each window. A later
     change can scale the intensities with temperature from HITRAN's lower-state
     energies.
   - OEP 1 records the infrared terms' 296 K beside the reference temperature (its §4.9),
     so a consumer can see the difference.
   - *Disposition (2026-09-30):* out of scope, as recommended. **Phase 3** states it in
     the manual.
9. **The validity of mixtures and air.** A trace constituent's band absorbs in
   proportion to its mole fraction. The intersection of the constituents' intervals,
   OEP 1's and OEP 4's rule for now, is therefore conservative: air at 400 µmol/mol of
   CO2 would warn across the whole 4.3 µm band.
   - **Recommended:** store each band's absorption at the reference state with the
     windows (Phase 1). A mixture's windows then follow from its constituents' bands
     weighted by their mole fractions, with the same criterion as a pure gas.
   - Offer that as the rule for OEP 1's air models and OEP 4's mixtures in Phase 3.
   - *Revised (2026-09-30) by Q13:* absorption alone cannot widen a mixture's
     intervals. **Phase 1** stores the absorption and bounds, and **Phase 3** the rule.

The following questions record GPT-6's review of all of Part I on 2026-09-30. The
recommendation is to **revise before acceptance**. Q10-Q13 are major scientific or
validity issues; Q14-Q16 concern the reference, interface and verification design;
Q17 records textual corrections. They require resolution or an explicit disposition
by the author before acceptance. The recommendations in Q1-Q9 remain the draft's
recommendations, not decisions made by this review.

10. **Weak in-window lines and the reference (§4, §8; Q3).** The reference retains
    every non-rotational line, while folding moves in-window contributions to poles
    outside the window. A reference line with positive strength has
    $\alpha_\ell = f_\ell/(E_\ell^2 - \omega^2)$, which diverges as $\omega$ approaches
    $E_\ell$. A model with no pole there remains finite. Their difference therefore
    has no finite bound over an interval containing that line, however small its
    strength; excluding the exact line centre does not bound its neighbourhood.
    The same problem applies to the dispersion derivatives.

    The feasibility calculation removes the in-window lines from both sides. An
    independent calculation from the cached HITRAN inputs reproduced the strong
    convergence for the remaining lines of CO and CO2. For the quoted CO2 window
    4.65-12.2 µm, however, 69,849 vibrational lines lie strictly inside it. With eight
    nodes per band and 1500 geometric wavelength samples, the median GVD error is
    about $10^{-12}$ of the electronic GVD against the truncated reference, but
    about $2\times10^6$ when the omitted lines are restored to the reference. This
    diagnoses the proposed undamped reference, not the physical GVD peaks, and is
    not a calculation of the proposed folding algorithm.

    Decide whether the target is resolved refractivity or an explicitly smoothed
    background. Resolved qualification needs excluded resonance neighbourhoods and
    a physically regularized reference where necessary. A background model needs a
    stated averaging operation and qualification against that observable. Recompute
    feasibility against the complete chosen reference before claiming the folding
    error is bounded.

    **Disposition (2026-09-30): resolved in the text.**
    - The target is the smoothed background: the line sum convolved with a Gaussian of
      stated standard deviation $\sigma$ in frequency (recommended 1 cm⁻¹, a 5.3 ps
      gate), evaluated with Dawson's function (§4). The resolved refractivity is
      rejected as a target (Option E).
    - The reference includes every line, in and outside the windows, and the rotational
      lines.
    - The feasibility is recomputed against it (§8). The review's $2\times10^{6}$ is
      reproduced at $\sigma = 0$. It falls to 56 electronic GVDs at 1 cm⁻¹ with every
      in-window band folded. The remaining error is the weak bands' smooth envelopes.
    - So folding is not always adequate. Each in-window band is folded only where the
      fold is within the tolerance, and is otherwise kept, which splits the window.
    - §8 gives the qualified parts that follow. CO2's laser region, and most of N2O's
      and CH4's mid-infrared candidates, do not qualify.
11. **Transparency requires a spectral criterion (§4, §7; Q1).** A band's absorption
    averaged over its entire line span does not bound its local peaks. Applying the
    threshold to each band separately also misses the sum of overlapping bands.
    The cached CO2 band containing the line at 2059.781613 cm⁻¹ (4.854884 µm)
    averages 0.07435 m⁻¹ at the proposed ideal-gas reference density, below the
    0.1 m⁻¹ threshold. That line alone gives a Lorentz-profile peak estimate of
    1.676 m⁻¹ at 1 atm and 296 K, using its intensity $2.017\times10^{-22}$ cm/molecule
    and self-broadened half-width 0.095 cm⁻¹/atm. It lies inside the quoted window.

    The calculation uses HITRAN's
    [intensity, linewidth and absorption definitions](https://hitran.org/docs/definitions-and-units/):
    $a(\tilde\nu) = N\sum_\ell S_\ell\phi_\ell(\tilde\nu)$, with consistent units.
    The higher-pressure statement in §2 also needs qualification: pressure changes
    linewidths as well as density, so the spectral absorption does not simply scale
    with density. Define total absorption at a stated spectral resolution, pressure,
    temperature and broadening composition. A band average can screen candidate
    windows, but does not establish their physical qualification. Line shapes could
    be used in the qualification tool without adding a public complex-index model.

    **Disposition (2026-09-30): resolved in the text.**
    - The criterion is the total absorption summed over every band, at the model's
      resolution: $a_\sigma = N\sum_\ell S_\ell g_\sigma$. It is HITRAN's definition
      with a Gaussian profile of the stated $\sigma$, for the pure gas at 101325 Pa and
      296 K, below 0.1 m⁻¹ (§4, §7).
    - Pressure and broadening composition do not enter while the collisional widths are
      below $\sigma$, which is the stated limit.
    - Each window also records the largest single-line Lorentz peak with HITRAN's
      self-broadened width, for users of narrowband light. Phase 1 adds the far-wing
      absorption at the reference state.
    - The finding is confirmed: CO2's 4.70-5.18 µm exceeds 0.1 m⁻¹, reaching 1.8 m⁻¹
      at 4.81 µm, with a 1.7 m⁻¹ Lorentz peak at 4.859 µm (§8).
    - §2's pressure statement is corrected: the smoothed absorption scales with
      density, the Lorentz peak does not, and the wings grow with its square.
    - The band average is kept only as a screening.
12. **Reassess the rotational exclusion in the extended range (§3, §4).** TC-156
    separates static orientation from polarizability and finds the rotational terms
    negligible at optical frequencies. That does not establish their negligibility
    in the mid infrared. Using the same band classification as
    `optmat.gas.infrared`, summing the cached pure rotational lines with the oscillator
    formula, and comparing with the stored fits' electronic terms at 10 µm gives:

    | Gas | Rotational polarizability / electronic polarizability | Rotational GVD / electronic GVD |
    |---|---|---|
    | CO | -77.3 ppm | -33.5 % |
    | N2O | -25.8 ppm | -11.7 % |

    These contributions are outside the proposed representation budget because the
    reference excludes them. Assess the dynamic rotational response explicitly:
    include it with a clear treatment of the static convention, or quantify its
    omission and restrict qualification accordingly. This review does not choose a
    replacement for the static convention.

    **Disposition (2026-09-30): resolved in the text; included.**
    - The rotational lines are part of the reference.
    - A polar molecule's pure rotational band enters the parent as its own band system,
      marked rotational (§4). As a 4-node rule it keeps CO's 5.39-12.5 µm error below
      $10^{-12}$ in $n - 1$ (§8).
    - Its omission would be 20 times CO's tolerance on $\alpha$ and 1.7 times N2O's,
      although only 0.15 % of the total GVD.
    - **The static convention is unchanged:** the static polarizability, the static
      validators and the electric validators leave the rotational system out, since
      its static limit is the orientation term (TC-156).
    - CH4's rotational band is about 3 % of its tolerance on $\alpha$; CO2 has none.
13. **Mixture qualification must also check the representation (Q9; OEP 4).** A
    dilute constituent may absorb weakly enough for the mixture to be transparent
    where its pure-gas model is unqualified. That model can still have quadrature
    poles there. Multiplying a pole by a small, nonzero mole fraction does not remove
    its divergence. Expanding intervals using absorption alone could admit an
    artificial singularity, especially when OEP 4 describes the replacement rule as
    changing warnings without changing values.

    Expanded mixture intervals need absorption, absence-of-poles and representation
    error checks for the mixture. Keep the intersection rule unless those checks
    succeed; some expansions require a different representation. Reconcile Q9 with
    OEPs 1 and 4.

    **Disposition (2026-09-30): resolved in the text, and reconciled with OEPs 1 and
    4.**
    - A mixture's intervals are the intersection of its constituents' intervals. They
      widen only where all three checks pass for the mixture at its mole fractions:
      - its smoothed absorption, weighted by mole fraction, is below the threshold;
      - no constituent's pole lies in the wider interval;
      - each constituent's representation error, times its mole fraction, is within
        the mixture's tolerance.
    - Phase 1 stores each band's absorption and each constituent's bound beyond its own
      windows, so the checks can be made.
    - Widening the intervals changes where the law warns, not its values.
    - Crossing a trace constituent's kept band, such as CO2's laser bands in air, needs
      a different representation of that constituent, with the band folded. That is a
      different model, with its own identity. It changes values, and is Phase 3 work
      under the same checks.
14. **SF6 needs a distinct reference and error assessment (§5; Q7).** The input is a
    measured, broadened absorption spectrum, not discrete transition strengths. The
    primary cached file is at 298.1 K, 760 Torr, with N2 broadening and 0.112 cm⁻¹
    instrumental resolution. It has nonzero cross-sections through substantial
    portions of the regions between strong bands. Quadrature is appropriate away
    from the measure's support; treating its samples as undamped lines inside weakly
    absorbing support introduces artificial poles.

    Specify a principal-value/Kramers-Kronig reference where necessary, and assess
    resolution, baseline and spectral truncation errors. Carry the actual measurement
    state; HITRAN's
    [cross-section documentation](https://hitran.org/docs/cross-sections-definitions/)
    limits transferring these measurements between thermodynamic conditions. Either
    define this treatment separately or defer SF6.

    **Disposition (2026-09-30): deferred.**
    - SF6 leaves this OEP's scope (§2). Its model, representation and single interval
      are unchanged, and a later proposal takes up the Kramers-Kronig reference.
    - A related finding does not wait for that proposal. SF6's cross-sections show
      weak absorption at 1.5-5 µm, 0.4-3.9 m⁻¹ at 1 atm, inside its declared
      230-5000 nm range (§3). Phase 1 records it in the register, with whether it is
      real absorption or baseline, and whether the declared range should narrow.
15. **Define the representation-error contract (§4, §6; Q4).** The definition
    $D_m = \partial^m(\omega\alpha)/\partial\omega^m$ makes $D_0 = \omega\alpha$,
    whereas the interface identifies order zero with $\alpha$. OEP 2 explicitly
    defines $D_0 = \alpha$ separately; use a consistent definition here.

    A signed deterministic residual, an error bound and a standard deviation are
    different quantities. The derivation supplies a residual, but
    `representation_sigma_au` suggests a standard deviation; OEP 2 already combines
    this quantity with statistical sigma in quadrature. Define the returned quantity,
    its units for every supported order, and the justified combination rule. Define
    behaviour outside qualified intervals as well.

    Interpolating sampled residuals does not establish a bound between samples.
    Specify independent interpolation checks or conservative bounds, and demonstrate
    the tolerance that justifies returning zero in the optical interval after adding
    near-infrared bands.

    **Disposition (2026-09-30): resolved in the text; the interface in Phase 2, with
    OEP 2.**
    - $D_0 = \alpha$ and $D_m = \partial^m(\omega\alpha)/\partial\omega^m$ for $m = 2, 3$,
      as in OEP 2 (§4).
    - The method is `representation_bound_au(omega, order)`. It returns a bound on the
      signed residual's magnitude, in the atomic units of $\alpha$ times $E_h^{1-m}$
      (§6).
    - It is combined with the statistical sigma linearly, not in quadrature. OEP 2's
      weights are corrected to the statistical sigma alone.
    - Outside every qualified interval it returns NaN and warns.
    - The bound between samples is specified: grid spacing rules, a bound per cell, and
      an independent offset grid that must lie within it (§7).
    - No interval returns zero by assumption. The optical interval's bound is computed
      with the near-infrared and rotational bands, like any other.
16. **Identify the full derivation and test higher moments (§8, §12).** Hashing the
    existing band table and settings does not identify the raw line distribution
    needed for higher-order quadrature. Different distributions can share its
    recorded low-order summaries while giving different rules and window errors.
    Record raw-input digests, the derivation version and precise settings, and retain
    a reproducible feasibility calculation with exact endpoints and grid definitions.

    The proposed invariants check only $M_{-2}$ and $M_0$, although an $n$-node rule
    claims $2n$ exact moments. Add higher-moment and direct-response checks using
    small synthetic measures that require no HITRAN access. Cover duplicate nodes,
    fewer than $n$ distinct lines, splitting and folding, and independent samples
    between qualification-grid points.

    **Disposition (2026-09-30): resolved in the text; Phase 1.** The derived tables
    carry the raw inputs' digests, the derivation's version, the tool's commit, every
    setting, and each window's endpoints and grid. The feasibility calculation is
    retained in `tools/hitran/`. The invariants include $\mu_0 \ldots \mu_{2n-1}$, and
    synthetic measures cover each listed case (§8).
17. **Correct the inventory and domain wording (§2, §5, §6; Q1-Q9).** The statement
    that CO is transparent throughout 2-4 µm conflicts with §2's own exclusion of
    2.29-2.52 µm. The current CO2 fit configuration includes
    `marchetti_2007/experiment_CO2`, but not the 2006 CO2 datum; correct its role and
    assess shared apparatus and systematic errors before treating it as independent
    validation. The 293.15 K statement overgeneralizes TC-730: CO and CH4 declare no
    such reference temperature, and SF6's fitted level rests mainly on 373 K data.

    Internal gaps can result from representation failure or the minimum window width,
    as well as absorption. Preserve the exclusion reason rather than labelling every
    gap an absorption band. Give each unresolved decision in Q1-Q9 an explicit phase
    assignment before Part I is frozen.

    **Disposition (2026-09-30): resolved in the text.**
    - A 2-4 µm spectrum in CO now crosses the 2.3 µm overtone band (§2).
    - Marchetti's 2006 value is described as curated but unused, and is assessed with
      the 2007 values before it validates anything (§3, §5).
    - The reference temperatures are stated per gas (§2).
    - Each gap records its reason, and the warning states it (§4, §6).
    - Q1-Q9 each carry a disposition and a phase: Q1, Q2, Q3 and Q6 in Phase 1; Q4 and
      Q5 resolved; Q7 deferred with SF6; Q8 in Phase 3; Q9 in Phases 1 and 3.

18. **A coarser resolution for the optical ranges.** Should the optical ranges be
    qualified against a second, coarser σ rather than not enforced? Open; raised in
    Phase 1, for the roadmap.
19. **Old's points.** Old's 1695 and 1817 nm points and CO2's 1.4-2.3 µm lie in a
    representation gap (TC-752). Fit them, compare them, or narrow the gap by keeping
    bands there? Open; raised in Phase 1, for Phase 3 and the author.
20. **The infrared validators.** The candidate infrared data of TC-756 are to be
    requested and curated. Open; raised in Phase 1.
21. **Accepting the refits (Phase 3).** Do the refits, with CH4 to 1313 nm and N2O to
    1950 nm, replace the present fits? *Decided (2026-10-01):* yes (TC-797).
22. **Marchetti's values.** CO2's last window ends at 8.94 µm, so its fitted
    9.2-10.7 µm values lie outside its intervals. Keep fitting them, or fold the laser
    bands so the window can qualify? *Decided (2026-10-01):* validators; the bands
    stay kept (TC-797).
23. **A second iteration of the windows.** Rederive the windows from the refitted
    parents' sigma? *Decided (2026-10-01):* no (TC-797).

### 14. References

1. A. A. Voronin and A. M. Zheltikov, "The generalized Sellmeier equation for air,"
   *Sci. Rep.* **7**, 46111 (2017). doi:10.1038/srep46111 (`voronin_2017`)
2. R. J. Mathar, "Calculated refractivity of water vapor and moist air in the atmospheric
   window at 10 µm," *Appl. Opt.* **43**, 928 (2004). doi:10.1364/AO.43.000928
   (`mathar_2004`)
3. R. J. Mathar, "Refractive index of humid air in the infrared: model fits," *J. Opt. A:
   Pure Appl. Opt.* **9**, 470 (2007). doi:10.1088/1464-4258/9/5/008 (`mathar_2007`)
4. K. T. Tang, "One-term Sellmeier formula for dispersion of dilute gases," *J. Opt. Soc.
   Am.* **62**, 644 (1972). doi:10.1364/JOSA.62.000644 (`tang_1972`)
5. K. T. Tang and K. K. Poon, *J. Opt. Soc. Am.* **64**, 1582 (1974).
   doi:10.1364/JOSA.64.001582 (`tang_1974`)
6. I. E. Gordon, L. S. Rothman, R. J. Hargreaves et al., "The HITRAN2024 molecular
   spectroscopic database," *J. Quant. Spectrosc. Radiat. Transfer* **353**, 109807
   (2026). doi:10.1016/j.jqsrt.2026.109807. A data service, not a curated source record;
   the provenance of the derived tables is `data/external/hitran/PROVENANCE.md`.
7. J. G. Old, K. L. Gentili and E. R. Peck, *J. Opt. Soc. Am.* **61**, 89 (1971)
   (`old_1971`)
8. S. Marchetti and R. Simili, *Infrared Phys. Technol.* **50**, 78 (2007).
   doi:10.1016/j.infrared.2006.07.003 (`marchetti_2007`)
9. G. H. Golub and J. H. Welsch, "Calculation of Gauss quadrature rules," *Math. Comp.*
   **23**, 221 (1969). **Not curated**; a numerical method. Its Eq. 2.6 was checked
   against the paper in GPT-6's review (2026-09-30); a source record is curated in
   Phase 1.
10. A. C. Simmons, *Opt. Commun.* **25**, 211 (1978). doi:10.1016/0030-4018(78)90309-7
    (`simmons_1978`)
11. D. V. Willetts and M. R. Harris, *J. Phys. D: Appl. Phys.* **15**, 51 (1982).
    doi:10.1088/0022-3727/15/1/007 (`willetts_1982`)
12. S. Marchetti and R. Simili, *Infrared Phys. Technol.* **47**, 263 (2006).
    doi:10.1016/j.infrared.2004.12.005 (`marchetti_2006`)

## Part II — Implementation record

### Phase 1 — The representation and the windows

- **Scope:** Phase 1 of §12, 2026-09-30, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0003-windows-data`, merged into `main`. Derived data only: no model or
  fit changes.
- **Delivered:**
  - `tools/hitran/windows.py`: the derivation of §4 and §7.
    - The line loader.
    - The smoothed reference. The far lines are FFT convolutions on a 0.005 cm⁻¹ grid; the lines within
      30 cm⁻¹ and everything above the line lists are summed directly.
    - The smoothed absorption and the Lorentz far wings.
    - The band systems, split at the windows, with the rotational system.
    - The adaptive rules, and the fold or keep choice per window.
    - The windows with their gaps' causes, and the bound checked between its grid points.
    - The collisional-width assessment.
    - The tables' provenance: the tool's commit, the raw inputs' SHA-256 and their query record, the parent fit's digest, and
      every setting.
  - `tools/hitran/feasibility.py`: §8's two tables, recomputed from the cache.
  - For CO, CO2, CH4 and N2O, in `src/optmat/data/external/hitran/`: `<gas>_band_rules.csv`,
    `<gas>_windows.csv`, `<gas>_bound.csv` and `<gas>_absorption.csv`, and the
    derived tables' section of `PROVENANCE.md`.
  - `optmat.gas.infrared`:
    - `gauss_rule` (Golub and Welsch's Eq. 2.6 on the Lanczos Jacobi matrix);
    - `fold_poles`;
    - the loaders `band_systems`, `infrared_windows`, `bound_cells` and `derived_provenance`.
  - The source record `golub_1969`.
  - Tests:
    - `tests/test_infrared.py`: the rules and the fold, on synthetic measures;
    - `tests/test_infrared_windows.py`: §8's invariants of the committed tables;
    - `tests/test_hitran_windows_tool.py`: the tool's numerics on synthetic lines.
  - {doc}`../../data`, *Infrared bands (HITRAN)*: a paragraph on the rules and windows.
- **Deviations from the proposal:**
  - **The optical ranges are split only by absorption (TC-749, J. C. Travers).** At
    σ = 1 cm⁻¹ the lines of weak near-infrared bands give the reference GVD structure
    that no pole model follows: CO2's at 1064 nm is 1600 electronic GVDs, at 1542 nm 39.
    - Enforced there, the tolerance would have ended CO2's, N2O's and CH4's ranges near 0.9 µm,
      with new warnings at fitted data and common laser lines.
    - Each fit's `range_nm` therefore keeps its data's qualification, split where the pure gas absorbs.
    - The tolerance qualifies the windows beyond it.
    - The bound is stored everywhere and reported inside the range.
    - An in-window band inside the range is always folded.

    §4's "split where a near-infrared band absorbs or cannot be represented" holds beyond the
    ranges only.
  - **The analysis starts at 12.5 µm** (800 cm⁻¹), the long end of §9's absorption map.
    No window longer than that is sought.
  - **An aggregated fold.** The weak bands of one window are folded together onto the
    window's two neighbouring nodes. The fold is linear, so this equals folding each band. Above the highest
    window no node is left; the parent's lowest electronic energy takes the upper pole,
    keeping M_0 and M_-2 exactly.
  - **The keep choice** is made among keeping none or the bands above 1e-5, 1e-6,
    1e-7 or 1e-8 of the molecule's M_0, per candidate window. It is the widest span these give,
    not an exhaustive search over bands.
  - **The bound's stored cells** are geometric (ratio 1.001), each holding the largest error of the
    0.25 cm⁻¹ grid cells inside it and their checks, rounded up to four digits.
  - **The rotational system includes CO2's** 836 rotational lines of its asymmetric
    isotopologues. §7 says CO2 has none; their static sum is 4.5 × 10⁻⁸ au (TC-751).
  - **A gap's reason** is the condition failing over most of it. The record gives every
    condition's share (`causes`), since absorption and representation alternate within
    one gap.
- **Decisions:** Q1, Q2, Q3 and Q6 are settled with the recomputed windows (TC-750, TC-748).
  - The rotational system is one rule per gas, left out of the static and electric validators in Phase 2 (TC-751).
  - marchetti_2006 is a reproduction check of marchetti_2007, not an independent validator (TC-753).
- **Register entries:** TC-747 to TC-756.
- **Verification:**
  - The full gate, `uv run poe check`, on the branch before the merge.
  - On synthetic lines the reference and its three derivatives agree with the direct sum over every line to
    1e-8, and the absorption with the direct Gaussian sum to 1e-5.
  - `feasibility.py` reproduces §8's first table for CO's 2.52-4.34 µm window, and its second
    table's 98 % (core 100 %).
  - The committed tables pass §8's invariants: positive strengths, nodes in their spans, each system's
    M_0 and M_-2 kept, the rules' totals equal to the band tables', no pole in a window, and
    the bound tiling every window.
  - Every cell of every window was checked between its grid points (TC-748).
  - The collisional widths were assessed at 1, 3 and 10 atm (TC-755).
  - Not run: `test-fits` and `test-library`, since no fit and no extraction changed.
- **Sources checked:**
  - Golub and Welsch's Eq. 2.6 and the Theorem, against the page image and the text.
  - Marchetti and Simili's 2006 and 2007 apparatus and uncertainties, against their records.
  - The rotational and vibrational static sums, against TC-156, bishop_1982 and the band tables.
  - The literature search for infrared refractivity found no measurement of CO or N2O
    between 2 and 12 µm, and one broadband set for CH4 (Rollefson and Havens, 1940). None was curated, since the papers are not in hand (TC-756).
- **Limitations:**
  - The model is the background at σ = 1 cm⁻¹, not the resolved lines.
  - The windows use the present fits' electronic terms and statistical sigma. Phase 3 recomputes them after the refits.
  - Line mixing, far-wing shapes, missing weak lines and HITRAN's intensity uncertainties
    are stated, not quantified.
  - The HITRAN line lists refer to 296 K (Q8).
  - SF6 absorbs inside its declared range (TC-754); it is deferred.
- **Open questions:**
  - Q18. Should the optical ranges be qualified against a coarser reference (a second σ)
    rather than not enforced? For the roadmap.
  - Q19. Old's 1695 and 1817 nm points and CO2's 1.4-2.3 µm lie in a representation gap
    (TC-752). Are they fitted, compared, or is the gap narrowed by keeping bands there? For
    Phase 3.
  - Q20. The candidate infrared data of TC-756 are to be requested and curated as validators. For Phase 3.
- **Handoff:**
  - Phase 2 adds the configuration keys and the stored fits' `intervals_nm`, reading
    `infrared_windows`, and `representation_bound_au` from `bound_cells`.
  - Phase 3 takes the fixed terms from `band_systems`: every rule and fold, with the rotational system as its own term class.
  - Rerun `tools/hitran/windows.py` after the refits, since the tolerance follows the parent's
    statistical sigma. The tables' provenance names the parent's configuration digest, so a
    freshness test can compare them.

### Phase 2 — Several intervals per model

- **Scope:** Phase 2 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0003-intervals`, merged into `main`. The machinery the refits use; no
  model or fit changes.
- **Delivered:**
  - `[model] band_nodes` and `windows = true` in the fit configurations
    (`optmat.gas.fitting`);
  - the rotational system as a term class of its own;
  - the stored record's `intervals_nm`, `gap_reasons`, `n_rotational`, `band_nodes`
    and `derived_sha256`, written only by a fit that uses the keys;
  - `optmat.gas.infrared.qualified_windows` and `derived_digest`;
  - `FittedModel.gap_reasons`, `n_rotational`, `derived_sha256` and
    `representation_bound_au`, and the gap's reason in the range warning;
  - `tests/test_infrared_intervals.py`; the per-interval stored-fit test; the fitting
    and identity pages.
- **Deviations from the proposal:**
  - `band_nodes` is added beside `infrared_max_wavenumber` rather than replacing it,
    so SF6's and O3's configurations keep theirs, as §6's plan foresaw.
  - The derivation's settings stay in the tables' provenance rather than in a
    separate file.
  - `band_nodes` must reach the tables' node count: a smaller one needs the tables
    rederived.
- **Decisions:** the rotational terms enter refractive observations and harmonic
  differences, and not electric ones or the moments (TC-770, TC-751).
- **Register entries:** TC-770.
- **Verification:**
  - The full gate. Every stored fit's digest, values and warnings are unchanged (the
    golden values).
  - The keys are tested on temporary copies of CO's and CO2's configurations, and the
    bound and the gap warning on CO2's fit given the derived windows.
- **Sources checked:** none new.
- **Limitations:** no stored fit uses the keys until Phase 3's refits.
- **Open questions:** none new.
- **Handoff:** Phase 3 adds `band_nodes = 8` and `windows = true` to CO, CO2, CH4 and
  N2O, refits them in dependency order, and reruns `tools/hitran/windows.py` after.
  Its merge waits for the author, since it changes recommended values.

### Phase 3 — The refitted molecules

- **Scope:** Phase 3 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0003-refits`, merged into `main` on the author's acceptance (Q21).
- **Delivered:**
  - CO, CO2, CH4 and N2O refitted with `band_nodes = 8` and `windows = true`;
  - their robustness records, and the regenerated pages;
  - the molecules' intervals table on {doc}`../../gas_models`, tested against the
    derived tables;
  - the static convention on {doc}`../../conventions`;
  - the gallery example *CO2 from the ultraviolet to its infrared windows*.
- **Deviations from the proposal:**
  - The mixture rule of Q13 is not implemented. It needs the mixture's own
    representation tolerance, which no OEP defines. The full CIPM-2007 air therefore
    takes the intersection of its constituents' intervals, now four, set by its trace
    CH4 and N2O; the reduced air keeps one.
  - No infrared validator was added: the candidates of TC-756 are not in hand (Q20).
- **Decisions made:**
  - Christiansen et al.'s electronic-plus-vibrational static CO is compared as an
    electric observation.
  - A windowed model's golden grid is its first interval.
  - Tests that sampled a model's envelope sample its intervals.
- **Results (TC-797):**
  - The fits' reduced chi-squared: CO 0.367 -> 0.367, CO2 3.38 -> 4.01 (Marchetti's
    values now validators, Q22), CH4 0.696 -> 0.685, N2O 3.014 -> 3.015.
  - Inside the new intervals the polarizability moves by at most 0.05, 0.26, 0.67 and
    1.86 sigma; CH4's above-one-sigma change is a finding.
  - CH4's optical range ends at 1313 nm, and N2O's at 1950 nm.
  - CO's and N2O's static values gain their orientation part.
  - CO2's last window is 8.4-8.9 µm: Marchetti's 9.2-10.7 µm values lie beyond it,
    where the laser bands are kept. They are now validators (Q22).
- **Register entries:** TC-797.
- **Verification:**
  - The full gate, `test-fits`, the comparison with the present fits over the optical
    intervals, and the robustness records.
  - The intervals table recomputed from the derived tables.
- **Sources checked:** none new.
- **Limitations:**
  - The windows were derived with the present fits' statistical sigma (Phase 1). They
    are not rederived after the refits; the fits' sigma changes by 0 to 5 % (CO2 +4-5 %),
    which would move the windows' tolerance by as much.
  - Q18 (a coarser sigma for the optical ranges) and Q19 (Old's points) remain open.
- **Questions decided by the author (2026-10-01):**
  - Q21: the refits are accepted.
  - Q22: Marchetti's values become validators; the laser bands stay kept.
  - Q23: the windows are not rederived.
- **Handoff:** OEP 2 is to be built on these parents once they are accepted, and solprop's
  parity fixtures for the four molecules are to be compared again (OEP 1, §6.11).

### Before Final — the open questions

- **Scope:** the dispositions of Q13 and Q18-Q20, 2026-10-01, by Opus 5.5 (AI) for John
  Travers, on `feature/oep-0003-final`, on the author's acceptance of the
  recommendations. With them every question has a disposition, and the OEP is Final.
- **Q13, a mixture's intervals (TC-809):** implemented. The intersection widens where
  the mixture passes the three checks at its mole fractions, with the mixture's
  tolerance defined as this OEP's own: 0.1 of its statistical sigma in $\alpha$ (the
  constituents' combined as independent) and 1 % of the larger of its $|D_2|$ and its
  electronic $|D_2|$. A new table per molecule, `<gas>_gap_bound.csv`
  (`tools/hitran/gap_bound.py`), holds the stored fit's representation error and the
  absorption in its gaps, against the same smoothed reference.
  - Deviation: Phase 1 stored the bound in the windows only, so the gaps' tables are new
    derived data, made from the refitted parents.
  - Result: the full CIPM-2007 air gains 1499-1516 and 1548-1562 nm. Its trace gases'
    other gaps stay closed, because their band-rule poles fill them; below N2O's
    320 nm nothing is widened, since its ultraviolet absorption is not assessed.
    Crossing them needs the traces in another representation (the roadmap).
- **Q18, a coarser resolution (TC-808):** decided: no change to the intervals.
  `tools/hitran/resolution.py` measures the fits' error inside their optical ranges at
  1, 10 and 100 cm⁻¹. The group-velocity dispersion holds to 1 % (this OEP's tolerance)
  up to about 900 nm at every resolution; at 100 cm⁻¹ to 1.47 µm (CO), 1.18 µm (CO2)
  and 1.14 µm (N2O), and for CH4 not beyond 940 nm, where the fit's own poles bound its
  narrow gaps. {doc}`../../accuracy` states it.
- **Q19, Old's points (TC-807):** compared, not fitted. Old et al.'s 1695 and 1817 nm
  points are a validator of `fit/CO2` with a wavelength window, at the dataset's fitted
  scale (a validator may now take `range_nm`, and one that is also fitted is compared
  at its fitted scale). They lie 68 and 96 ppm below the fit (10 and 14 of their
  sigma), about the gap's representation error, so they test the bands' rules there
  rather than the electronic terms.
- **Q20, the infrared validators (TC-802 to TC-806):** the five candidates are curated.
  Rollefson and Havens's CH4 (1.68-14.8 µm; inside the intervals only at 11.46 and
  11.89 µm), Matsumoto's CO2 at 10.6 µm and Tsankova et al.'s static CO are validators.
  Legay (1955) measured HCl alone, and Mitrofanov et al.'s group delay of air is
  recorded without a dataset (no quantity for a group delay whose zero is not stated).
  Against the fits: CH4 +0.70 % on average (the authors' stated level), CO2 at
  10.6 µm -0.91 % (within its 2 %), CO's static value +110 ppm (0.02 sigma). The fits
  are unchanged; TC-807.
- **Verification:** the full gate, `test-fits`, the library lane for the new
  extractions and the release test.
- **Sources checked:** the five papers, each against its page images (TC-802 to TC-806).
- **Limitations:** SF6 stays deferred (Q7, Q14), as does a complex index for the weak
  bands (Option F); both are on the roadmap. Statescu (1915), CO2's infrared index by
  the prism method, is wanted.

## Part III — History

- 2026-09-30: Drafted, with the feasibility calculation of §2 and §8. Takes up OEP 2's
  former Q3 (a band's second moment) and the CO2 window deferred from OEP 2.
- 2026-09-30: Aligned with OEP 1 and OEP 4. OEP 1 declares `intervals` (one per model)
  and the related identity, temperature and override provisions (its Q30); Q5 resolved
  there. Q9 added: the validity of mixtures and air. Phase 2 now builds on OEP 1's
  Phase 1, and the header requires OEP 1.
- 2026-09-30: GPT-6 (AI) reviewed all of Part I against main at
  `ba368030570ff3f88322f6a7aabfb3c6dc9911df`, the related OEPs, model and fit
  configurations, source records, manual, register and proposal process. Recorded
  unresolved findings as Q10-Q17; the recommendation is to revise before acceptance.
  The Gaussian quadrature construction outside a band's support is supported, but the
  treatment of in-window lines and the definition of qualification require revision.
  Numerical checks read the cached HITRAN line lists, independently reproduced the
  principal convergence trend for CO and CO2, restored the omitted in-window lines,
  estimated a CO2 absorption peak, summed the dynamic rotational contributions, and
  inspected the SF6 cross-section input. They do not reproduce every entry of §8 or
  implement the proposed folding algorithm. Primary sources checked: HITRAN's official
  line-parameter and cross-section definitions (linked in Q11 and Q14), its HITRAN2024
  citation, and Golub & Welsch [9], especially the eigenvector-weight relation in
  Eq. (2.6), checked in the
  [paper](https://csclub.uwaterloo.ca/~pbarfuss/GolubWelsch.pdf). The numerical-method
  reference remains uncurated. The measurement tables were not all freshly checked;
  no fits were regenerated and the full development gate was not run for this review.
  The read-only proposal/register checks returned 47 passed and one expected skip.
  No data, fit or physical convention was changed. The human review remains required
  for `Reviewed`; this record does not accept the proposal or resolve its questions.
- 2026-09-30: John Travers (human) reviewed Part I. With GPT-6's review, the review
  status is `Reviewed`.
- 2026-09-30: Revised for GPT-6's review (Opus 5.5 drafting); every one of Q10-Q17 has
  a disposition.
  - **The target** is the background at a stated resolution: the line sum smoothed by a
    Gaussian, recommended $\sigma = 1$ cm⁻¹, against a reference with every line.
  - **In-window bands** are folded or kept per band. The fold goes to two poles, one on
    each side of the window.
  - **Transparency** is the smoothed total absorption. The band average is only a
    screening.
  - **The rotational band** enters as its own system, and the static convention is
    kept.
  - **The GVD tolerance** is relative to the larger of the total and the electronic GVD.
  - **The representation bound** replaces the representation sigma, with a linear
    combination rule and a specified check between grid points.
  - **SF6** is deferred.
  - **Mixtures** widen their intervals only under the three checks of Q13.
  - **Provenance** carries the raw-input digests, and synthetic tests cover the higher
    moments.
  - §8 gains the feasibility against the complete reference.
    - CO qualifies almost entirely.
    - CO2 qualifies over most of 2.2-2.6 and 2.9-4.2 µm, and over 5.54-6.95 µm.
    - N2O's and CH4's mid-infrared candidates mostly do not.
    - Pure CO2 is not qualified at its laser lines.
  - Q1-Q9 each carry a phase. OEP 1, OEP 2 and OEP 4 were aligned in the same change.
  - The calculations were made from the cached HITRAN2024 line lists, SF6's cached
    cross-sections and the stored fits. No data, fit or physical convention was
    changed.
- 2026-09-30: Accepted by John Travers. Part I is frozen.
- 2026-09-30: Implementation started. Phase 1 (the representation and the windows)
  merged, with its report in Part II (TC-747 to TC-756).
- 2026-10-01: Phase 2 (several intervals per model) merged, with its report in Part II
  (TC-770).
- 2026-10-01: Phase 3 (the refitted molecules) merged, with its report in Part II
  (TC-797).
- 2026-10-01: Q13 and Q18-Q20 settled on the author's acceptance of the
  recommendations (TC-802 to TC-809). Final.

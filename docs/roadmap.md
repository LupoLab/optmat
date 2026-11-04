# Roadmap

```{note}
Planned work is tracked in [GitHub issues](https://github.com/LupoLab/optmat/issues), which can be filtered by type, priority and
area. It is not listed here. This page holds only the status lists of accepted and draft
OEPs below.

The entries this page used to hold moved to issues on 2026-11-04, and the papers and data
we are looking for moved to {doc}`their own page <wanted>`. Where an older page, a source
record or an archived OEP says that something is "on the roadmap", look in those. The
[page as it stood before the move](https://github.com/LupoLab/optmat/blob/24672b1fcc5e7209406b4a37cd50437befa50445/docs/roadmap.md)
shows what it said.
```

optmat began with the gases used in gas-filled fibres and cells, and with the properties
a simulation of nonlinear propagation in them needs first: the refractive index across
the spectrum and the Kerr response. The aim is the same for every property and every
material that follows. Take the data from the primary literature with their provenance,
check them on record, fit physical models with honest uncertainties, and keep the
published alternatives available for comparison.

A major item starts as an {doc}`enhancement proposal <oeps/index>` (OEP). That means a new
medium, a new family of properties or models, a changed default or convention, or a
redesigned interface. Accepted proposals are listed first, then drafts.

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
  its section 5. Tracked in [#35](https://github.com/LupoLab/optmat/issues/35).
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
  target (TC-1056). Start the solprop investigation ([#10](https://github.com/LupoLab/optmat/issues/10)) with OEP 6's continuum coefficients:
  compare pump-ionisation baselines and evaluate perturbative VUV pair-spectral loss
  on the generated fields. OEP 8's lowest bound lines are at different wavelengths.
  This diagnostic need not wait for OEP 8, and it does not provide the continuum's
  missing time-dependent loss current or decide which mechanism dominates.
  optmat supplies material properties and basic data; solprop owns the resulting
  polarisation, electron accounting, plasma coupling and propagation energy balance
  (TC-1057). The internal atomic fixtures test material coefficients and conventions.
  Tracked in [#36](https://github.com/LupoLab/optmat/issues/36).

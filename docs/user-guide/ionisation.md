# Ionisation

Begin with the parameters that define a rate, then read the noble-gas and molecular
models with their validation. The two-photon chapter describes continuum cross
sections and the curated information about bound excitation; these have different
uses and limitations from a strong-field rate.

The {doc}`/quickstart` shows how to evaluate a rate. Check the field convention,
wavelength, polarisation, and stated support before applying it. The guide then
continues with {doc}`solids`.

```{toctree}
:maxdepth: 1

../ionisation
../ionisation_rates
ionisation-theory
ionisation-static
ionisation-laser
ionisation-one-electron
ionisation-measured
ionisation-accuracy
../molecular_ionisation_rates
molecular-ionisation-validation
../two_photon
two-photon-data
two-photon-excitation
```

See the {doc}`ionisation gallery </auto_gallery/ionisation/index>` for comparisons
with calculations and measurements. Reference interfaces include
{mod}`optmat.gas.ionisation_rates`, {mod}`optmat.gas.molecular_ionisation_rates`,
{mod}`optmat.gas.two_photon`, and {mod}`optmat.gas.two_photon_bound`.

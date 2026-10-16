# Conventions

These conventions apply throughout the user guide. The {doc}`quickstart` shows
them in a calculation; {doc}`introduction/curation` explains how a source's own
conventions are preserved in its record. Continue with {doc}`identity` to see how
a model records its provenance, status, and supported range.

## Units

Every interface is SI: metres, seconds, pascals, kelvin, number densities in m$^{-3}$
and polarizabilities in C m$^2$/V. Wavelengths are **vacuum** wavelengths; angular
frequency is $\omega = 2\pi c/\lambda_\mathrm{vac}$. The units of the literature
(atmospheres, torr, ångström, atomic units, wavenumbers) are converted once, where a
source's values are loaded. One atomic unit of polarizability is
$4\pi\varepsilon_0 a_0^3$ (`optmat.constants.AU_POLARIZABILITY`).

## States and densities

A measurement of the refractive index is made at some pressure and temperature and is
usually *reduced* to a reference state: 0 °C and 760 Torr in most of the older
literature, 15 °C and 1 atm in atmospheric science, 0 °C and 1 bar in Börzsönyi et al.
(2008). A reduction is made either to the ideal-gas density at that state or to the
real-gas density. To compare the two kinds, every value is turned into a polarizability
with the number density its author meant (`optmat.gas.eos.number_density`, `ideal=True`
or not). Real-gas densities are $p/(Z k_B T)$ with the compressibility factor $Z$ of the
CoolProp equation of state for that gas. CoolProp comes with the `gas` extra
(`pip install "optmat[gas]"`); without it, every model, air's included, still gives
its polarizability and its index at the ideal-gas density (`ideal=True`), and asking
for the real gas raises an `ImportError` that names the extra (TC-1092).

**The equation of state as an object.** `optmat.gas.eos.equation_of_state(gas)` returns
it with its identity (CoolProp's version, the fluid and the equation's reference key;
{doc}`identity`), its forward law, its speed of sound and its inverse,
`pressure(N, T)`, which solves $p = N k_B T\,Z(p, T)$ to a round trip of 10$^{-12}$ in $N$
(OEP 1, section 4.2). Both keep the CODATA convention: $N$ is $p/(Z k_B T)$ with
CODATA's $k_B$, not the equation's mass density, which carries the gas constant of its
publication (argon's, 5.7 ppm from CODATA's).

**The molar mass.** `molar_mass` is the molar mass the backend's equation uses, kg/mol:
the one its mass density and its speed of sound are formed with, for a propagation
code's mass terms. The number densities do not use it. For the pseudo-pure dry air of
`air_equation_of_state()` it is CoolProp's 28.96546 g/mol, which is not the
28.9586 ± 0.0002 g/mol Lemmon et al. (2000) print, the mean of their composition; the
two differ by 2.4 × 10⁻⁴ (TC-832). A mixture's equation of state gives the mean of its
constituents' molar masses ({doc}`mixtures`).

**The gas phase domain.** The gas equation-of-state interface accepts gases and
supercritical fluids: the forward law and its inverse share
one domain, the gas below its vapour pressure (to 10$^{-5}$ below it, where CoolProp's
pressure-temperature solve still converges) and the supercritical fluid up to the
equation's stated pressure, within its temperatures above the triple point
(`phase_domain(T)`). Where the melting pressure in CoolProp is lower than that, the fluid
ends 10$^{-5}$ below the melting pressure instead: nitrogen up to 284 K, argon up to
254 K, methane up to 256 K and CO2 up to 328 K (TC-811). A condensed or solid state, or
one beyond the equation, raises `StateError` in both directions: SF6 at 294 K
condenses at 21.4 bar, CO2 at 280 K at 41.6 bar. None of the states of the curated
datasets is condensed (TC-731).

```python
from optmat.conditions import State
from optmat.gas.eos import equation_of_state

sf6 = equation_of_state("SF6")
n = sf6.number_density(State(10e5, 294.0))
assert abs(float(sf6.pressure(n, 294.0)) / 10e5 - 1) < 1e-12
low, high = sf6.phase_domain(294.0)
assert 21e5 < high < 22e5
```

Wavelengths given in *standard air* (dry air with 0.03 % CO$_2$ at 15 °C and
101 325 Pa, Leonard's and Cuthbertson's convention) are converted to vacuum with
Edlén's 1966 formula (`optmat.conditions.vacuum_wavelength_from_standard_air`).

## Temperature

A state's temperature enters the index through the density (the equation of state) and
Lorentz-Lorenz. The polarizability itself is taken as independent of temperature. For an
atom that holds at the package's accuracy. A molecule's zero-density polarizability
changes a little with temperature, through its rotational and vibrational averaging:
per K, 1.18(4) ppm for N2 (Egan & Yang 2023), 2.11 ppm for O2 and 4.22 ppm for CO2 (Egan
2026, at 633 nm), and 0.03 % for H2 between 273 and 296 K. A fit that declares
`reference_temperature_k` refers to that temperature, and its datasets measured
elsewhere are brought to it where the coefficient is measured. N2, O2, CO2, N2O, H2 and
D2 refer to 293.15 K; CH4 and CO rest on 297.5 K refractometry and SF6 mostly on 373 K
values, with no curated coefficient (TC-730). For pulse propagation the effect is
usually negligible: 20 K from the reference is at most about 85 ppm (CO2), against the
few 0.1 % to which a gas fill's density is typically known.

A fitted model exposes these as metadata, and they enter its identity
({doc}`identity`): `reference_temperature`, `temperature_coefficient` (the measured
coefficient with its uncertainty, wavelength and dataset, for N2, O2 and CO2) and
`infrared_temperature`, the temperature of its fixed infrared terms (296 K for HITRAN's
intensities, 298.1 K for SF6's cross-section). Evaluation does not follow them (OEP 1,
section 4.9).

```python
from optmat.gas.models.registry import model

co2 = model("fit/CO2")
assert co2.reference_temperature == 293.15
assert co2.infrared_temperature == 296.0
assert co2.temperature_coefficient.coefficient_per_k == 4.22e-6
```

## Polarizability and refractive index

The dispersion of a gas is a property of its atoms or molecules: the dynamic dipole
polarizability $\alpha(\omega)$. At number density $N$ the refractive index follows from
the Lorentz-Lorenz relation

$$
\frac{n^2-1}{n^2+2} = \frac{N\alpha}{3\varepsilon_0} \equiv y, \qquad
\chi = n^2 - 1 = \frac{3y}{1-y}.
$$

Linear density scaling, $\chi = N\alpha/\varepsilon_0$, is its dilute limit and differs
from it by the factor $1/(1-y) \approx 1 + \chi/3$: 2×10$^{-4}$ relative for argon at
1 bar, but about 2 % at 100 bar.

Beyond Lorentz-Lorenz, pairs of atoms or molecules add the second refractivity virial
coefficient $B_R$:

$$
\frac{n^2-1}{n^2+2} = A_R\,\rho + B_R\,\rho^2 + \dots,
$$

with $\rho$ the molar density and $A_R = N_A\alpha/3\varepsilon_0$ the molar
refractivity. The correction $B_R\rho/A_R$ is $-6\times10^{-6}$ for helium and
$+1.9\times10^{-5}$ for argon at 1 atm, and grows in proportion to the density. Every
conversion in the package (a measured refractivity to a polarizability, a model's
polarizability to $n$) includes it, with the $B_R$ chosen for each gas in
`data/gas/refractivity_virial.toml` (Achtermann et al. 1991 at 633 nm where measured,
Achtermann et al. 1993 for neon and krypton, Burns et al. 1986 for xenon). Oxygen takes
Egan et al.'s (2026) published estimate, 0.210 cm³/mol at 633 nm, with 50 %
uncertainty: no measured value is known. Gases without a curated value (deuterium,
nitrous oxide, ozone, carbon monoxide) get no pair term.

A mixture follows the same law with a sum over species and pairs,
$y = \sum_s A_s\rho_s + \sum_{s,t} B_{st}\,\rho_s\rho_t$, inverted once for the whole
mixture. The unlike terms $B_{st}$ are estimates with a declared uncertainty
({doc}`mixtures`).

## The infrared terms and the static limit

A molecule's polarizability carries its infrared bands as fixed poles: Gaussian rules of
HITRAN's lines for CO, CO2, CH4 and N2O (OEP 3; {doc}`fitting`). The polar molecules,
CO and N2O, also carry their pure rotational band as a system of its own. At optical
frequencies it contributes nothing that matters. In the static limit it is the
orientation polarizability $\mu^2/3k_BT$: 0.65 a.u. for CO and 1.45 a.u. for N2O at
296 K, 4.7 and 6.5 % of their static polarizabilities.

The models' static value is therefore the static dielectric one: electronic,
vibrational and orientational. A static value without the orientational part (an ab
initio electronic-plus-vibrational polarizability) is compared as an *electric*
observation, which leaves the rotational system out. Electronic values alone are
compared without the infrared terms at all.

## The magnetic term

At optical frequencies $n^2 = \varepsilon\mu$, so a refractive index measures the electric
polarizability plus the atom's (diamagnetic) magnetic susceptibility $\chi_0$: $-2.12
\times 10^{-5}$ a.u. for helium, $-8.5 \times 10^{-5}$ for neon, $-2.30 \times 10^{-4}$ for
argon, 15-30 ppm of the polarizability. The package's polarizability is the *refractive*
one, the sum; ab initio calculations and dielectric-constant measurements give the
electric one alone, and the fits add $\chi_0$ to them (TC-71).

The refractive polarizability is what the forward-scattering amplitude gives, and beyond
the dipole approximation that amplitude has one more term: its retardation, the
quadrupole and magnetic-dipole parts at order $(a_0 k)^2$. For helium, Pachucki &
Puchalski (2019) give it as $\kappa\,\alpha_\mathrm{nrel}\,(a_0 k)^2$ with $\kappa =
-0.507$: in atomic units $-3.73\times10^{-5}\,\omega^2$, $-0.14$ ppm of $n - 1$ at 633 nm
and $-24$ ppm of the group-velocity dispersion. `fit/He` adds it, like $\chi_0$, to its
refractive data and to its model and not to its electric data or moments (TC-705). For the
other gases it is not known; the same relative size would be about $10^{-5}$ of their
dispersion, below their uncertainties, and it is neglected.

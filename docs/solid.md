# Solids: conventions and models

Glasses, crystals and metals share one subpackage, `optmat.solid` (TC-318). A solid has no gas
state: its optical response is a property of the material, its grade and its
temperature. This chapter defines the conventions and common complex-index
interface. Start with the solid example in {doc}`quickstart`, then follow
{doc}`silica`, {doc}`silicon`, {doc}`crystals` or {doc}`metals` for model choices and evidence.
The {doc}`source catalogue <the-data/sources>` and {doc}`data` format provide the
recorded citations, sample conditions and corrections.

## Conventions

- **The complex index:**
  - n + ik, with k >= 0 for absorption; the intensity falls as exp(-4 pi k z / lambda).
  - The relative permittivity is epsilon = (n + ik)^2, with Im epsilon >= 0.
  - This is the convention of Palik's handbook and of refractiveindex.info.
  - Henke's n = 1 - delta - i beta (X-ray tables) has the opposite sign of the
    imaginary part.
- **Wavelengths** are vacuum wavelengths, in m at every interface.
- **Absolute indices.** Every model returns the absolute index (relative to vacuum).
  Minimum-deviation refractometry in room air gives the index relative to that air,
  about 2.7e-4 lower. A dataset or formula says which with `index_medium` (`vacuum` or
  `air`, with the air's temperature and pressure); a formula relative to air is
  multiplied by the index of that air, from Edlen's dry-air formula at its temperature
  and pressure ({func}`optmat.conditions.air_refractivity`; TC-398). Malitson's and
  Rodney & Spindler's indices are relative to air (TC-399), and so are Dodge's.
  Daimon & Masumura measured CaF2 in nitrogen: a formula says `index_medium =
  "nitrogen"`, with the nitrogen's temperature and pressure, and is multiplied by the
  index of optmat's N2 model at the ideal-gas density (TC-545).
- **Frequencies:** omega = 2 pi c / lambda, in rad/s.
- **Photon energies:** E = hbar omega / e, in eV (CODATA 2018 through scipy).
- **Temperature:** a solid dataset records its sample temperature as
  `parameters.temperature_k` when the source states it. A model states the temperature
  it refers to.
- **Grade:** a glass's grade (Corning 7980, Suprasil 1, Lithosil Q2, ...) is recorded
  per dataset. Grades of fused silica differ in OH and metal impurities and in thermal
  history. This moves the ultraviolet absorption edge and the 2.7 µm OH band, and n at
  the level of 1e-5.
- **Anisotropic crystals:** a dataset of a uniaxial or biaxial crystal names its axis
  (`o`, `e`; `n_p`, `n_m`, `n_g` or `x`, `y`, `z`), and so does a model.
- **Films** differ from bulk material (density, stoichiometry, porosity). Their data
  carry `form = "film"`, and a bulk model uses them only as stated in its register
  entry.

The materials are registered in `data/materials.toml` ({mod}`optmat.materials`); see
{doc}`data`, *Materials*.

## The model interface

Every model of a solid ({class}`optmat.solid.models.base.Model`) gives, at vacuum
wavelengths:

- `permittivity(wavelength)`: epsilon;
- `dpermittivity_domega(wavelength)`: the exact derivative d epsilon / d omega;
- `refractive_index(wavelength)`: n + ik, computed so that a small k keeps its full
  relative precision (a plain complex square root returns k = 0 below about
  1e-16 n);
- `dindex_domega`, `group_index` and `gvd`.

A model warns ({class}`~optmat.exceptions.ModelRangeWarning`) outside the range its
source supports. This matches what a propagation code needs of a bulk medium or a
capillary cladding: epsilon(omega) and its derivative, a name, a citation and a range.

The models of a material, and the recommended one with the reason, come from
{mod}`optmat.solid.models.registry` (`models`, `model`, `default_model`,
`default_reason`; `data/solid/defaults.toml`). An anisotropic crystal's recommendation
is per axis, `default_model("MgF2", "o")`, and asking without the axis is an error
(TC-546).

## Published formulas

A source record's `[[formulas]]` entry for a solid becomes a model when it has one of
these forms ({mod}`optmat.solid.models.formulas`):

| Form | Formula | Terms |
|---|---|---|
| `n2_sellmeier` | n^2 - 1 = A + sum_i B_i lambda^2 / (lambda^2 - C_i), lambda in um | `A` (optional), `B`, `C` (um^2) |
| `n2_sellmeier_wavelength` | the same, with the poles printed as wavelengths, C_i = L_i^2 | `B`, `L` (um) |
| `n2_sellmeier_temperature` | n^2 - 1 = sum_i S_i(T) lambda^2 / (lambda^2 - lambda_i(T)^2), S_i and lambda_i polynomials in T (K) | `S1`, `S2`, ... and `L1`, `L2`, ... (lowest power first), with `reference_temperature_k` and `temperature_range_k` |
| `udm` | Franta et al.'s universal dispersion model (below) | `parameter_file` |

The coefficients are the printed strings, corrected where the record says why, as for
the gases.

A temperature-dependent formula is evaluated at its `reference_temperature_k`;
`at_temperature(T)` gives it at another temperature, and warns outside the stated
range (Leviton & Frey's fused silica, 30-300 K: `leviton_2008/sellmeier`).

## Franta's universal dispersion model

The universal dispersion model (UDM) of Franta, Nečas and co-workers builds the complex
dielectric function of a solid from its excitations (Proc. SPIE 9628, 96281U (2015),
`franta_2015b`):

- the interband transitions, a symmetric term plus excitonic terms;
- the higher-energy excitations of the valence electrons;
- the Urbach tail;
- Gaussian-broadened absorption on localised states and phonons.

Each contribution's imaginary part is a simple model of a joint density of states,
normalised to unit transition strength (the integral of epsilon_i E over E). Its real
part is the Kramers-Kronig transform in closed form. The model is therefore causal and
passive by construction, and the transition strengths satisfy the Thomas-Reiche-Kuhn
sum rule.

{mod}`optmat.solid.udm` implements it from the paper's equations. Its parameters are the
values of the authors' parameter files from their program newAD2, which the source
records hold (TC-320). It is checked against the authors' own tables and program,
against 50-digit arithmetic, and against a numerical Kramers-Kronig transform (TC-321):

- above 1 eV it reproduces the published n and k of all seven files to half a unit of
  their seventh printed digit;
- below 1 eV the published tables carry the authors' program's approximations (at most
  2 units of the seventh digit for silica), and optmat's values are the exact ones.

Two details follow the program rather than the paper, because the program reproduces
the published tables and the parameters were fitted with it:

- the Urbach tail ends at E_h, where the paper's Eq. (33) defines E_m = (E_g + E_h)/2;
- the phonon wavenumbers are converted with the program's constants.

**Why a model and not a table.** A spline of a tabulated n and k is continuous, but its
derivatives are the spline's, not the material's. It gives noisy group-velocity
dispersion and wrong derivatives near absorption features. Franta's table *is* a
sampling of this analytic model; evaluating the model itself removes the problem
exactly, with no smoothing.

- **Differentiability:** the closed forms are analytic between the few band edges where
  a joint density of states has a kink (E_g, E_h, E_xi, ...). There epsilon remains
  continuous, and its higher derivatives have at most logarithmic singularities, at
  isolated energies inside absorption bands.
- **Derivatives:** d epsilon / dE is exact to rounding (complex-step differentiation).
- **Low energies:** far below each band the real part is evaluated from its moment
  series, which avoids the cancellation that costs the closed forms digits.

## Tables of optical constants

A metal has no single model from the ultraviolet to the far infrared. Its constants are
measured piecewise, each range by a different technique: electron energy-loss spectra,
ellipsometry, reflectance with a Kramers-Kronig analysis, far-infrared reflectance. A
**table model** ({mod}`optmat.solid.models.tabulated`, named `table/<material>`) joins
curated datasets, one pair per wavelength range, as `data/solid/tables/<material>.toml`
chooses them, each range with its reason (TC-544, TC-645):

```toml
material = "Ag"

[[segments]]
real = "babar_2015/..."        # refractive_index, refractive_index_decrement or
imaginary = "babar_2015/..."   # permittivity_real; its extinction_coefficient or
range_nm = [70.5, 12400.0]     # permittivity_imaginary
reason = "..."
```

- A segment takes the wavelengths its two datasets share. An x-ray table may print n to
  some energy and 1 - n (`refractive_index_decrement`) above it, with k throughout; a
  permittivity pair is converted to n + ik at the points (TC-645).
- In the soft x-ray range a segment may instead be Henke et al.'s atomic scattering
  factors (`henke = "<element>"`, CXRO's files in `data/external/henke/`) with a curated
  `atom_density` dataset: n + ik = 1 - r_e λ² N (f1 - i f2)/(2π). This independent-atom
  limit fails near the plasma and absorption edges; no recommended table uses it yet
  (TC-645).
- The ranges may not overlap. Nothing is averaged: where two sources cover the same
  wavelengths, the configuration chooses one and says why.
- The model interpolates log n and log k against log omega with a monotone
  piecewise-cubic Hermite interpolant (PCHIP). Unlike a natural cubic spline, it
  neither overshoots between points nor rings at a join between sources. Its first
  derivative is continuous and exact for the interpolant, so epsilon and
  d epsilon / d omega agree.
- Outside the tabulated range the constants are held at the end points, with a
  {class}`~optmat.exceptions.ModelRangeWarning`.

Group-velocity dispersion from a table is the interpolant's, not the material's:
between the points it carries the table's noise and the interpolant's piecewise second
derivative. For a metal this rarely matters (a metal is a reflector or a cladding, not
a propagation medium), but the caveat of the previous section applies.

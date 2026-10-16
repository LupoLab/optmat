# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

# Record Luna.jl's material parameters of the noble gases, H2 and D2 as a test fixture.
#
#     julia --project=<Luna checkout> tools/luna_reference.jl
#
# runs in *Luna's* environment (this package never calls Julia) from the root of this
# repository and writes `tests/data/luna_master.json`. The Luna checkout should be its
# master branch; the commit is read with git, or from the environment variable
# LUNA_COMMIT when the checkout has no git metadata (a `git archive` export, whose commit
# `git get-tar-commit-id` prints). For each gas, at 1 bar and 20 degrees C unless stated:
#
# - `wavelength`, `refractivity` (n - 1) and `beta2` (s^2/m): Luna's `ref_index` and
#   `dispersion(2, ...)` on a grid from 200 nm to 2500 nm;
# - `refractivity_800nm_by_pressure`: n - 1 at 800 nm at 1, 10 and 100 bar, which shows
#   Luna's density scaling;
# - `density`: Luna's number density (m^-3);
# - `gamma3` (m^5/V^2, chi3 = gamma3 * N) and `ionisation_potential` (J);
# - `rate`: the ionisation rate Luna's `prop_capillary` uses with plasma on (PPT with its
#   default Stark-shift and dipole factors for the noble gases, ADK for H2 and D2), driven
#   at 800 nm: at each peak intensity (W/cm^2, the field amplitude F from
#   I = c eps0 F^2 / 2), the rate at the peak field and its average over a cycle of
#   F cos(phi); and `peak_bare`, the rate at the peak field of one electron (occupancy
#   1) without Luna's extra factors (PPT: m = 0 only, no Stark shift, no dipole factor;
#   ADK: no threshold), for comparison with the forms of the literature here;
# - `raman` (H2 and D2): the numerical fields of Luna's `raman_parameters`.
#
# `tools/luna_tables.py` and `tests/test_luna_reference.py` compare them with the models here.

using Luna
import Luna.PhysData: ref_index, dispersion, density, γ3_gas, ionisation_potential,
    raman_parameters, c, ε_0
import Luna.Ionisation: IonRatePPT, IonRateADK

const GASES = (:He, :Ne, :Ar, :Kr, :Xe, :H2, :D2)
const ADK_GASES = (:H2, :D2)  # Luna's Interface.makeplasma! uses ADK for these
const P = 1.0  # bar
const T = 293.15
const OUT = joinpath(@__DIR__, "..", "tests", "data", "luna_master.json")

function commit()
    haskey(ENV, "LUNA_COMMIT") && return ENV["LUNA_COMMIT"]
    readchomp(`git -C $(pkgdir(Luna)) rev-parse HEAD`)
end

# A minimal JSON writer: a Vector of Pairs keeps the keys in order.
json(x::AbstractString) = "\"" * escape_string(x) * "\""
json(x::Symbol) = json(string(x))
json(x::Bool) = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::AbstractFloat) = isfinite(x) ? repr(Float64(x)) : "null"
json(x::AbstractVector{<:Pair}, ind="") = "{\n" * join(
    ["$ind $(json(k)): $(json(v, ind * " "))" for (k, v) in x], ",\n") * "\n$ind}"
json(x::AbstractVector, ind="") = "[" * join(json.(x), ", ") * "]"
json(x, ind) = json(x)

function cycle_average(rate, F; n=2000)
    # the rate is even in the field: average over a quarter cycle, at midpoints
    ϕ = ((1:n) .- 0.5) .* (π / 2 / n)
    sum(rate(F * cos(p)) for p in ϕ) / n
end

function rates(gas)
    λ0 = 800e-9
    adk = gas in ADK_GASES
    rate = adk ? IonRateADK(gas) : IonRatePPT(gas, λ0)
    bare = adk ? IonRateADK(gas; occupancy=1, threshold=false) :
        IonRatePPT(gas, λ0; stark_shift=false, dipole_corr=false, msum=false, occupancy=1)
    intensity = 10 .^ range(12.5, 15.5, length=61)  # W/cm^2
    F = sqrt.(2 .* intensity .* 1e4 ./ (c * ε_0))
    [
        "form" => (adk ? "IonRateADK" : "IonRatePPT"),
        "wavelength" => λ0,
        "intensity_w_cm2" => collect(intensity),
        "peak" => [rate(f) for f in F],
        "cycle_averaged" => [cycle_average(rate, f) for f in F],
        "peak_bare" => [bare(f) for f in F],
    ]
end

# Luna's field names, in ASCII: Δα is delta_alpha, dαdQ dalpha_dq, ρ rho, Ω omega, μ mu.
ascii(name) = replace(string(name), "Δα" => "delta_alpha", "dαdQ" => "dalpha_dq",
    "ρ" => "rho", "Ω" => "omega", "μ" => "mu")

function raman(gas)
    rp = raman_parameters(gas)
    [ascii(k) => v for (k, v) in pairs(rp) if v isa Real && !(v isa Bool)]
end

function entry(gas)
    λ = collect(200e-9:25e-9:2500e-9)
    e = [
        "wavelength" => λ,
        "refractivity" => [real(ref_index(gas, l, P, T)) - 1 for l in λ],
        "beta2" => [dispersion(2, gas, l, P, T) for l in λ],
        "refractivity_800nm_by_pressure" => [
            "$p" => real(ref_index(gas, 800e-9, p, T)) - 1 for p in (1.0, 10.0, 100.0)],
        "density" => density(gas, P, T),
        "gamma3" => γ3_gas(gas),
        "ionisation_potential" => ionisation_potential(gas),
        "rate" => rates(gas),
    ]
    gas in (:H2, :D2) && push!(e, "raman" => raman(gas))
    e
end

function main()
    data = [
        "luna_commit" => commit(),
        "pressure_bar" => P,
        "temperature_k" => T,
        "gases" => ["$g" => entry(g) for g in GASES],
    ]
    write(OUT, json(data) * "\n")
end

main()

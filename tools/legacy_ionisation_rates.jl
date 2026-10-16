# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

# Record the ionisation rates Luna.jl's propagation applies to the noble gases at 800 nm,
# the first of the three steps that write `tests/data/legacy_ionisation_rates.json`:
#
#     julia --project=<Luna checkout> tools/legacy_ionisation_rates.jl
#     (cd <solprop checkout> && uv run python <optmat checkout>/tools/legacy_ionisation_rates.py solprop)
#     uv run python tools/legacy_ionisation_rates.py optmat
#
# The first runs in *Luna's* environment (this package never calls Julia) from the root of
# this repository. It writes the fixture afresh, with the field grid and the pulses the
# other two steps read, so it must run first; the second and third add solprop's and
# optmat's rates on the same fields (`tools/legacy_ionisation_rates.py`). The Luna commit
# is read with git, with whether its tree is clean, or from the environment variable
# LUNA_COMMIT when the checkout has no git metadata.
#
# The rate is the one `Interface.prop_capillary` uses with `plasma = true` (the default)
# for He, Ne, Ar, Kr and Xe: `Interface.makeplasma!` takes PPT for an atom and builds
# `Ionisation.IonRatePPTCached(gas, grid.referenceλ; PPT_options...)` with
# `PPT_options` empty, that is `IonRatePPTAccel(gas, λ0; cache = true)` with its
# defaults (stark_shift, dipole_corr, msum, sum_integral, occupancy 2, the cubic spline
# of log(rate) on 2^16 fields up to twice the barrier-suppression field, above which it
# throws). `Nonlinear.PlasmaCumtrapz` evaluates it at the instantaneous field. Here the
# table is built into a temporary directory, so that nothing is written to Luna's
# cache; it is the same table, and where Luna's own cache already holds a table under
# the same key, it is read and compared.
#
# For each gas:
# - `rate`: the rate (1/s) at the fields of `field_v_per_m`, the amplitudes of the
#   cycle-averaged intensities `intensity_w_cm2` (25 points per two decades from 1e13
#   to 1e15 W/cm^2, with 5e13, 2e14 and 5e14 added; I = c eps0 F^2 / 2), taken
#   as instantaneous values; null where Luna refuses the
#   field, with the reason in `refused`;
# - `yield`: 1 - exp(-int W dt) for Gaussian pulses of optmat's
#   `ionisation_validation.gaussian_pulse` (reimplemented here: `field = -dA/dt`, 96
#   samples per cycle over three FWHM either side of the peak, CEP 0), by the trapezoid
#   rule, at the peak intensities `peak_intensity_w_cm2` (15 from 2e13 to 1e15 W/cm^2,
#   with 5e13, 1e14, 2e14 and 5e14 added); null where the pulse reaches a field Luna
#   refuses.

using Luna
import Luna.PhysData: c, ε_0, ionisation_potential, polarisability,
    polarisability_difference, quantum_numbers
import Luna.Ionisation: IonRatePPT, IonRatePPTAccel, barrier_suppression
import Luna: Utils
import HDF5

const GASES = (:He, :Ne, :Ar, :Kr, :Xe)
const λ0 = 800e-9
# W/cm^2: log grids with the intensities the comparison is quoted at added
const QUOTED = [5e13, 1e14, 2e14, 5e14]
const INTENSITY = sort(unique(vcat(10 .^ range(13, 15, length=25), QUOTED)))
const PEAK_INTENSITY = sort(unique(vcat(10 .^ range(log10(2e13), 15, length=15), QUOTED)))
const FWHM = (6e-15, 30e-15)  # s, of the intensity
const SAMPLES_PER_CYCLE = 96  # optmat's ionisation_validation._SAMPLES_PER_CYCLE
const HALF_WIDTH = 3.0  # FWHM either side of the peak (_PULSE_HALF_WIDTH)
const OUT = joinpath(@__DIR__, "..", "tests", "data", "legacy_ionisation_rates.json")

function provenance()
    dir = pkgdir(Luna)
    haskey(ENV, "LUNA_COMMIT") && return [
        "commit" => ENV["LUNA_COMMIT"], "branch" => nothing, "clean" => nothing]
    [
        "commit" => readchomp(`git -C $dir rev-parse HEAD`),
        "branch" => readchomp(`git -C $dir rev-parse --abbrev-ref HEAD`),
        "clean" => isempty(readchomp(`git -C $dir status --porcelain`)),
    ]
end

# A minimal JSON writer, as in tools/luna_reference.jl: a Vector of Pairs keeps the keys
# in order; `nothing` is null.
json(::Nothing) = "null"
json(x::AbstractString) = "\"" * escape_string(x) * "\""
json(x::Symbol) = json(string(x))
json(x::Bool) = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::AbstractFloat) = isfinite(x) ? repr(Float64(x)) : "null"
json(x::AbstractVector{<:Pair}, ind="") = "{\n" * join(
    ["$ind $(json(k)): $(json(v, ind * " "))" for (k, v) in x], ",\n") * "\n$ind}"
json(x::AbstractVector, ind="") = "[" * join(json.(x), ", ") * "]"
json(x, ind) = json(x)

amplitude(intensity_w_cm2) = sqrt(2 * intensity_w_cm2 * 1e4 / (c * ε_0))

"""The times and the field of optmat's `gaussian_pulse` of peak amplitude F0."""
function gaussian_pulse(F0, fwhm; cep=0.0)
    ω = 2π * c / λ0
    period = λ0 / c
    half = HALF_WIDTH * fwhm
    count = 2 * ceil(Int, half / period * SAMPLES_PER_CYCLE) + 1
    t = collect(range(-half, half, length=count))
    shape = @. exp(-2 * log(2) * t^2 / fwhm^2)
    slope = @. -4 * log(2) * t / fwhm^2 * shape
    phase = @. ω * t + cep
    # -dA/dt of A = -(F0 / ω) e(t) sin(phase)
    t, @. F0 * (shape * cos(phase) + slope / ω * sin(phase))
end

trapezoid(y, t) = sum((y[1:end-1] .+ y[2:end]) ./ 2 .* diff(t))

"""The production rate, with the name of the file Luna's cache keys it under."""
function production_rate(gas)
    tmp = mktempdir()
    rate = IonRatePPTAccel(gas, λ0; cache=true, cachedir=tmp)
    rate, only(readdir(tmp))
end

function cached_comparison(rate, key, fields)
    path = joinpath(Utils.cachedir(), "pptcache", key)
    isfile(path) || return ["exists" => false]
    E, W = HDF5.h5open(path, "r") do file
        (read(file["E"]), read(file["rate"]))
    end
    cached = IonRatePPTAccel(E, W)
    inside = [f for f in fields if f <= min(rate.Emax, cached.Emax)]
    Pair{String,Any}["exists" => true,
        "max_relative_difference" => maximum(abs(cached(f) / rate(f) - 1) for f in inside)]
end

function rates(rate, fields)
    out = Union{Float64,Nothing}[]
    refused = Pair{String,Any}[]
    for (I, f) in zip(INTENSITY, fields)
        try
            push!(out, rate(f))
        catch err
            err isa ErrorException || rethrow()
            push!(out, nothing)
            push!(refused, "$I" => err.msg)
        end
    end
    out, refused
end

function yields(rate)
    entries = Pair{String,Any}[]
    for fwhm in FWHM
        values = Union{Float64,Nothing}[]
        refused = Pair{String,Any}[]
        for I in PEAK_INTENSITY
            t, field = gaussian_pulse(amplitude(I), fwhm)
            try
                push!(values, -expm1(-trapezoid(rate.(field), t)))
            catch err
                err isa ErrorException || rethrow()
                push!(values, nothing)
                push!(refused, "$I" => err.msg)
            end
        end
        push!(entries, "$(round(Int, fwhm * 1e15))fs" => [
            "yield" => values, "refused" => refused])
    end
    entries
end

function entry(gas, fields)
    rate, key = production_rate(gas)
    # the same rate evaluated directly, without the spline
    direct = IonRatePPT(gas, λ0)
    inside = [f for f in fields if f <= rate.Emax]
    W, refused = rates(rate, fields)
    Δα = polarisability_difference(gas)
    α_ion = polarisability(gas, true)
    [
        "ionisation_potential_j" => ionisation_potential(gas),
        "quantum_numbers_n_l_z" => collect(quantum_numbers(gas)),
        "stark_delta_alpha_si" => ismissing(Δα) ? nothing : Δα,
        "dipole_alpha_ion_si" => ismissing(α_ion) ? nothing : α_ion,
        "table_field_max_v_per_m" => rate.Emax,
        "table_field_max_over_barrier_suppression" =>
            rate.Emax / barrier_suppression(ionisation_potential(gas), 1),
        "spline_max_relative_difference_from_direct" =>
            maximum(abs(rate(f) / direct(f) - 1) for f in inside),
        "luna_cache_file" => key,
        "luna_cache_existing" => cached_comparison(rate, key, inside),
        "rate" => W,
        "refused" => refused,
        "pulses" => yields(rate),
    ]
end

function main()
    fields = amplitude.(INTENSITY)
    data = [
        "description" => "Instantaneous ionisation rates (1/s) and Gaussian-pulse " *
            "yields of the noble gases at 800 nm in Luna.jl, solprop and optmat, on the " *
            "same fields: written by tools/legacy_ionisation_rates.jl, then " *
            "tools/legacy_ionisation_rates.py solprop, then ... optmat.",
        "wavelength" => λ0,
        "intensity_w_cm2" => INTENSITY,
        "field_v_per_m" => fields,
        "pulses" => [
            "peak_intensity_w_cm2" => PEAK_INTENSITY,
            "peak_amplitude_v_per_m" => amplitude.(PEAK_INTENSITY),
            "fwhm_s" => collect(FWHM),
            "cep" => 0.0,
            "form" => "optmat.gas.ionisation_validation.gaussian_pulse: field = -dA/dt, " *
                "A = -(F0/omega) e(t) sin(omega t + cep), e = exp(-2 ln2 t^2/fwhm^2); " *
                "$SAMPLES_PER_CYCLE samples per cycle over $HALF_WIDTH FWHM either side",
            "checks" => ["$(round(Int, f * 1e15))fs" => begin
                t, field = gaussian_pulse(1.0, f)
                Pair{String,Any}["n" => length(t), "t_first" => t[1],
                    "sum_field_squared_unit_peak" => sum(abs2, field)]
            end for f in FWHM],
            "yield" => "1 - exp(-int W dt), trapezoid rule on the pulse's own grid",
        ],
        "luna" => [
            provenance()...,
            "julia_version" => string(VERSION),
            "function" => "Luna.Ionisation.IonRatePPTCached(gas, 800e-9) = " *
                "IonRatePPTAccel(gas, 800e-9; cache=true), as built by " *
                "Interface.makeplasma! for plasma=true (PPT_options empty)",
            "options" => "stark_shift=true, dipole_corr=true (both inactive for Kr and " *
                "Xe: no polarisability in Luna), msum=true, occupancy=2, " *
                "sum_integral=true, cycle_average=false, Cnl from the asymptotic " *
                "(Hartree) expression with l* = n* - 1; spline N=2^16 up to twice the " *
                "barrier-suppression field, error above it",
            "gases" => ["$g" => entry(g, fields) for g in GASES],
        ],
    ]
    write(OUT, json(data) * "\n")
end

main()

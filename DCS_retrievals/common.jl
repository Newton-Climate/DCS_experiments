#=

Paths and loaders shared by the analysis scripts (compare_*.jl, rmse_*.jl, calibrate_noise.jl)
for the batch runs of configs/*.toml (SpectralFits v0.2). Run the scripts from this directory:

    julia --project=../../SpectralFit.jl compare_hit24.jl

=#

using SpectralFits, Dates, Statistics, LinearAlgebra, Printf
import JLD2

# ---------------------------------------------------------------- locations

# override with the environment variables when running elsewhere (for example on fluo)
const DATA_DIR    = get(ENV, "DCS_DATA_DIR", joinpath(@__DIR__, "..", "data"))
const SPECTRA_DIR = get(ENV, "DCS_SPECTRA_DIR", joinpath(@__DIR__, "..", "spectra"))
const RUNS_DIR    = joinpath(DATA_DIR, "runs")     # output_root of configs/*.toml
const FIG_DATA    = joinpath(@__DIR__, "..", "fignewtons", "timeseries", "data")   # v0.1-era results

spectra_path(file) = isabspath(file) ? file : joinpath(SPECTRA_DIR, file)

"""Seconds between the DCS machine time and the `t` stored in the v0.1-era results in
fignewtons/timeseries/data (hit08_results.jld2, OCO_results.jld2, obs.jld2). The true
LabVIEW-to-Unix offset is 2082844800 s, so those `t` values are one day late; this constant
reproduces them so new and old records can be matched."""
const OLD_T_OFFSET_S = 2082758400.0

# ---------------------------------------------------------------- batch results

"""
    results_file(experiment; root = RUNS_DIR, input = "20160921") -> String or nothing

Results file for input `<input>.h5` in the most recent run folder of `experiment`
(`find_run`), or `nothing` when there is no such run or file.
"""
function results_file(experiment::AbstractString; root = RUNS_DIR, input = "20160921")
    folder = try find_run(root, experiment) catch e; e isa ArgumentError || rethrow(); return nothing end
    path = joinpath(folder, "results", "$(input)_results.jld2")
    return isfile(path) ? path : nothing
end

has_results(experiment; kwargs...) = results_file(experiment; kwargs...) !== nothing

"""
    load_run(experiment; root, input) -> Vector or nothing

Results of an experiment (`load_results` NamedTuples) without the `Failed` ones; `nothing`,
with a printed note, when the experiment has not been run.
"""
function load_run(experiment; kwargs...)
    path = results_file(experiment; kwargs...)
    path === nothing && (println("note: no results for experiment \"$experiment\", skipped"); return nothing)
    return filter(r -> r.status != Failed, load_results(path).results)
end

"""Time of a result in the convention of the v0.1-era results (see `OLD_T_OFFSET_S`)."""
old_t(r) = Float64(r.metadata["machine_time"]) - OLD_T_OFFSET_S

"""Pressure (hPa), temperature (K) and path length (cm) of a result: the retrieved p and T
when they are in the state, otherwise the conditions it was retrieved with (the measured or,
with a [conditions] table, the tower values)."""
function pT(r)
    c = r.metadata["conditions"]
    p = hasproperty(r.x, :pressure_hPa) ? r.x.pressure_hPa : Float64(c["pressure_hPa"])
    T = hasproperty(r.x, :temperature_K) ? r.x.temperature_K : Float64(c["temperature_K"])
    return (pressure_hPa = p, temperature_K = T, path_length_cm = Float64(c["path_length_cm"]))
end

"""Dry-air column (molecules/cm2) of a result: air column at its p, T and path minus H2O."""
function dry_column(r)
    c = pT(r)
    return homogeneous_air_column(c.pressure_hPa, c.temperature_K, c.path_length_cm) - r.x.H2O
end

"""Index of state label `label` (for example \"CH4\") in a loaded result."""
label_index(r, label) = findfirst(==(label), r.state_labels)

"""Observations not de-weighted by a mask: diagonal precision within 1e-3 of the largest."""
unmasked(r) = (d = diag(r.obs_precision); d .> 1e-3 * maximum(d))

"""Noise sigma of the first full-weight observation of a result."""
noise_sigma(r) = (k = findfirst(unmasked(r)); 1 / sqrt(r.obs_precision[k, k]))

# ---------------------------------------------------------------- v0.1-era records

"""First DCS record time of the campaign in the v0.1 convention (hit08_results.jld2 t[1]);
obs.jld2 stores t_pT and t_pic as days since it."""
campaign_t0() = JLD2.load(joinpath(FIG_DATA, "hit08_results.jld2"))["t"][1]

"Tower p, T (obs.jld2) nearest to each v0.1-convention time in `t`."
function tower_pT(t)
    o  = JLD2.load(joinpath(FIG_DATA, "obs.jld2"))
    tp = campaign_t0() .+ 86400 .* o["t_pT"]
    k  = [argmin(abs.(tp .- ti)) for ti in t]
    return o["p"][k], o["T"][k]
end

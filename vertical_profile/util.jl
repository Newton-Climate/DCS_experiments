#=

Shared helpers for the analysis and table-building scripts of the synthetic vertical-profile
experiments. The experiments themselves are ProfileOSSE configurations in configs/ (see
README.md); their run folders are in runs/. Run the scripts from this directory with the
SpectralFit.jl project:

    julia --project=../../SpectralFit.jl <script>.jl

=#

using SpectralFits, Statistics, LinearAlgebra, Printf
import SpectralFits.Apps.ProfileOSSE

const CONFIG_DIR  = joinpath(@__DIR__, "configs")
const RUNS_DIR    = joinpath(@__DIR__, "runs")
const SPECTRA_DIR = joinpath(@__DIR__, "..", "spectra")

# ---------------------------------------------------------------- configurations

"""
    profile_config(name, variant) -> ProfileOSSEConfig

Variant `variant` of configs/`name`.toml, for example `profile_config("us_standard_ch4", :isobar)`.
Loads the configured cross sections.
"""
profile_config(name, variant) =
    ProfileOSSE.profile_osse_config(joinpath(CONFIG_DIR, "$(name).toml"); variant = String(variant))

# ---------------------------------------------------------------- run folders

"""
    load_run(experiment; root = RUNS_DIR) -> NamedTuple

The most recent run of `experiment` (for example `"us_standard_co2_isobar_tables"`) in `root`:
`folder`, `result` (the retrieval: `x`, `K`, `averaging_kernel`, `dfs_per_element`,
`posterior_covariance`, `chi2_obs_per_dof`, ...), `attributes` (everything the app saved),
`x_true` and `prior_sigma` (as `ComponentVector`s), `vcd` (layer dry-air columns,
molecules/cm²), `z_km`, `pressure_hPa` and `temperature_K` (layers, after the mode
collapse), `mode` and `target`.
"""
function load_run(experiment::AbstractString; root = RUNS_DIR)
    folder = find_run(root, experiment)
    loaded = load_results(joinpath(folder, "results", "$(experiment)_results.jld2"))
    a = loaded.attributes
    return (; experiment, folder, result = only(loaded.results), attributes = a,
              x_true = state_from_labels(a["x_true"], a["x_true_labels"]),
              prior_sigma = state_from_labels(a["prior_sigma"], a["prior_sigma_labels"]),
              vcd = a["dry_air_column_cm2"], z_km = a["height_m"] ./ 1000,
              pressure_hPa = a["pressure_hPa"], temperature_K = a["temperature_K"],
              mode = Symbol(a["mode"]), target = Symbol(a["target"]))
end

# ---------------------------------------------------------------- diagnostics

"""Retrieved and true dry mole fractions of `gas` in run `r` (from `load_run`)."""
retrieved_vmr(r, gas::Symbol) = getproperty(r.result.x, gas) ./ r.vcd
true_vmr(r, gas::Symbol) = getproperty(r.x_true, gas) ./ r.vcd

"""Column averaging kernel of `gas` in run `r`, v0.1 convention `(H' A) ./ H'`, H = vcd / Σ vcd."""
column_ak(r, gas::Symbol) =
    column_averaging_kernel(r.result.averaging_kernel, state_indices(r.result.x, gas), r.vcd)

"""Per-layer degrees of freedom for signal of `gas` in run `r`."""
layer_dfs(r, gas::Symbol) = r.result.dfs_per_element[state_indices(r.result.x, gas)]

rms(v) = sqrt(mean(abs2, v))

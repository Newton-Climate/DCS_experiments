#=

Analysis of the cross-line-list synthetic profile experiments (configs/crosslist_ch4.toml,
configs/crosslist_co2.toml). The synthetic truth spectrum is computed with one line list for
the target gas and retrieved with another:

    CH4: truth HITRAN 2008      → retrieved with HITRAN 2008 (control), 2016, 2020, TCCON 2020
    CO2: truth OCO ABSCO 5.1    → retrieved with ABSCO (control), HITRAN 2016, 2020, TCCON 2020
    CO2: truth HITRAN 2020      → retrieved with HITRAN 2020 (control), 2016, TCCON 2020

Every run has the same noise realization (seed 1), so a retrieval minus its control (the
first retrieval of its truth, same mode) isolates the line-list error from noise and
smoothing. Run the experiments first (any subset; missing runs are skipped):

    julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/crosslist_co2.toml

then this script, which reads the most recent run folder of each variant in runs/:

    julia --project=../../SpectralFit.jl linelist_errors.jl [CH4|CO2]

The averaging kernel, per-layer DOF and posterior covariance are the ones the package stored
(Jacobi-scaled inversions); the hand-written σ-scaled diagnostics of the v0.2 hand-port are
no longer needed. `time_s` is the run time from the run folder (results file written minus
configuration written).

Outputs in runs/linelist_analysis/, per species (stem = species, or species_<modes> when not
every mode has run):
    linelist_summary_<stem>.csv   one row per run
    linelist_profiles_<stem>.csv  one row per run and layer
(the v0.2 hand-port outputs in linelist_results/ have the same columns).

=#

include("util.jl")

const OUT_DIR = joinpath(RUNS_DIR, "linelist_analysis")
const MODES = (:profile, :isotherm, :isobar)

# (truth, retrievals); the first retrieval is the control with the truth's own list
const CROSSLIST_EXPERIMENTS = Dict(
    :CH4 => [("HIT08", ["HIT08", "HIT16", "HIT20", "TCCON"])],
    :CO2 => [("ABSCO51", ["ABSCO51", "HIT16", "HIT20", "TCCON"]),
             ("HIT20",   ["HIT20", "HIT16", "TCCON"])],
)
const SCALE = Dict(:CH4 => 1e9, :CO2 => 1e6)          # ppb, ppm

const BANDS = (("0_2", 0.0, 2.0), ("2_5", 2.0, 5.0), ("5_10", 5.0, 10.0 + 1e-9))   # km, [lo, hi)

"""Experiment name of a crosslist variant (configs/crosslist_<species>.toml)."""
crosslist_experiment(species, truth, retrieval, mode) =
    "crosslist_$(lowercase(String(species)))_truth-$(truth)_retr-$(retrieval)_$(mode)"

# ---------- diagnostics from the stored retrieval ----------

function diagnostics(r, species)
    res = r.result
    t = state_indices(res.x, species)
    A, S = res.averaging_kernel, res.posterior_covariance
    col_true = sum(getproperty(r.x_true, species))
    col_ret = sum(getproperty(res.x, species))
    colAK = vec(sum(A[t, t], dims = 1))                 # Σᵢ A[i, j] = ∂column/∂x_j
    layer_dof = diag(A)[t]
    in_band(lo, hi) = findall(z -> lo <= z < hi, r.z_km)

    # condition numbers of the normal matrix, raw and in prior-σ-scaled coordinates
    N = Symmetric(Matrix(res.K' * (res.obs_precision * res.K)) + Matrix(res.prior_precision))
    d = collect(r.prior_sigma)
    folder_time = mtime(joinpath(r.folder, "results", "$(r.experiment)_results.jld2")) -
                  mtime(joinpath(r.folder, "configuration.toml"))

    return (; col_err_pct = 100 * (col_ret / col_true - 1), col_σ_pct = 100 * sqrt(sum(S[t, t])) / col_true,
              dof_total = sum(layer_dof),
              dof_band = [sum(layer_dof[in_band(lo, hi)]) for (_, lo, hi) in BANDS],
              colAK_band = [mean(colAK[in_band(lo, hi)]) for (_, lo, hi) in BANDS],
              χ² = res.chi2_obs_per_dof,
              vmr_true = SCALE[species] .* true_vmr(r, species),
              vmr_ret = SCALE[species] .* retrieved_vmr(r, species),
              layer_σ = SCALE[species] .* sqrt.(diag(S)[t]) ./ r.vcd, layer_dof, colAK,
              cond_raw = cond(N), cond_scaled = cond(Symmetric(d .* N .* d')), time = folder_time)
end

# ---------- CSV output ----------

const SUMMARY_HEADER = "species,truth,retrieval,mode,col_err_pct,dcol_vs_control_pct,col_sigma_pct," *
    "shape_rms_pct,dshape_rms_vs_control_pct,dshape_maxabs_vs_control_pct,dprofile_rms_vs_control_pct," *
    "dof_total,dof_0_2km,dof_2_5km,dof_5_10km,colAK_0_2km,colAK_2_5km,colAK_5_10km,chi2_per_dof," *
    "cond_raw,cond_scaled,time_s"
const PROFILE_HEADER = "species,truth,retrieval,mode,z_km,p_hPa,T_K,vmr_true,vmr_retrieved,err,err_pct," *
    "dvmr_vs_control,dvmr_vs_control_pct,posterior_sigma,posterior_sigma_pct,dof_layer,column_AK"

g(x) = @sprintf("%.6g", x)

function summary_row(species, truth, key, mode, dg, ctrl)
    pct  = 100 .* (dg.vmr_ret .- dg.vmr_true) ./ dg.vmr_true
    dpct = 100 .* (dg.vmr_ret .- ctrl.vmr_ret) ./ dg.vmr_true
    dcol = dg.col_err_pct - ctrl.col_err_pct
    vals = [dg.col_err_pct, dcol, dg.col_σ_pct, rms(pct .- dg.col_err_pct), rms(dpct .- dcol),
            maximum(abs, dpct .- dcol), rms(dpct), dg.dof_total, dg.dof_band..., dg.colAK_band..., dg.χ²,
            dg.cond_raw, dg.cond_scaled, dg.time]
    return join([String(species), truth, key, String(mode), g.(vals)...], ",")
end

function profile_rows(species, truth, key, mode, r, dg, ctrl)
    map(eachindex(dg.vmr_true)) do j
        e = dg.vmr_ret[j] - dg.vmr_true[j]
        dc = dg.vmr_ret[j] - ctrl.vmr_ret[j]
        join([String(species), truth, key, String(mode),
              g.([r.z_km[j], r.pressure_hPa[j], r.temperature_K[j], dg.vmr_true[j], dg.vmr_ret[j], e,
                  100e / dg.vmr_true[j], dc, 100dc / dg.vmr_true[j], dg.layer_σ[j],
                  100dg.layer_σ[j] / dg.vmr_true[j], dg.layer_dof[j], dg.colAK[j]])...], ",")
    end
end

write_csv(path, header, rows) = open(io -> (println(io, header); foreach(l -> println(io, l), rows)), path, "w")

# ---------- driver ----------

function analyse(species)
    summary, profiles, modes_found = String[], String[], Set{Symbol}()
    for (truth, retrievals) in CROSSLIST_EXPERIMENTS[species], mode in MODES
        ctrl = nothing
        for key in retrievals
            experiment = crosslist_experiment(species, truth, key, mode)
            r = try
                load_run(experiment)
            catch e
                e isa ArgumentError || rethrow()
                @info "no run, skipped" experiment
                continue
            end
            dg = diagnostics(r, species)
            if ctrl === nothing
                key == first(retrievals) || @warn "control $(first(retrievals)) has not run; Δ columns use $key" experiment
                ctrl = dg
            end
            push!(modes_found, mode)
            push!(summary, summary_row(species, truth, key, mode, dg, ctrl))
            append!(profiles, profile_rows(species, truth, key, mode, r, dg, ctrl))
            @printf "%s %-8s truth %-7s retr %-7s col err %+.3f%% (Δctrl %+.3f%%)  DOF %.2f [%s]  χ²/dof %.3f\n" species mode truth key dg.col_err_pct dg.col_err_pct - ctrl.col_err_pct dg.dof_total join(map(x -> @sprintf("%.2f", x), dg.dof_band), ", ") dg.χ²
        end
    end
    isempty(summary) && return @info "no crosslist runs of $species in $RUNS_DIR"
    found = [m for m in MODES if m in modes_found]
    stem = length(found) == length(MODES) ? String(species) : join([String(species); String.(found)], "_")
    mkpath(OUT_DIR)
    write_csv(joinpath(OUT_DIR, "linelist_summary_$(stem).csv"), SUMMARY_HEADER, summary)
    write_csv(joinpath(OUT_DIR, "linelist_profiles_$(stem).csv"), PROFILE_HEADER, profiles)
    println("wrote ", joinpath(OUT_DIR, "linelist_{summary,profiles}_$(stem).csv"))
end

species_list = filter(s -> String(s) in ARGS, (:CH4, :CO2))
foreach(analyse, isempty(species_list) ? (:CH4, :CO2) : species_list)

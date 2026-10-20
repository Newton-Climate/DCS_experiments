#=
Spectral-fit RMSE (measurement minus model, full-weight channels: `r.rmse`) for the HITRAN
2024 comparison experiments on 20160921 (configs/hit24.toml and hit24_co2.toml). OCO ABSCO
residuals are not stored in the campaign results, so its RMSE is recovered exactly from the
stored chi squared: sqrt(chi2 (n-m)/n) sigma (the CO2 window has no masks and the same noise
as co2_hit24). Experiments that have not been run are skipped with a note.
compare_hit24.jl includes this file for spectral_rmse / oco_rmse.
Usage: julia --project=../../SpectralFit.jl rmse_hit24.jl
=#
isdefined(@__MODULE__, :results_file) || include(joinpath(@__DIR__, "common.jl"))

"""Experiment names of the HITRAN 2024 comparison, from its configs."""
hit24_experiments() = [c["experiment"] for f in ("hit24.toml", "hit24_co2.toml")
                       for c in SpectralFits.read_config(joinpath(@__DIR__, "configs", f))]

"Per-record spectral RMSE of an experiment (full-weight channels), or nothing when not run."
function spectral_rmse(n)
    rs = load_run(n)
    rs === nothing && return nothing
    rel = [(k = unmasked(r); 100 * sqrt(mean(((r.y .- r.y_model) ./ r.y)[k] .^ 2))) for r in rs]
    return (rs = rs, rmse = [r.rmse for r in rs], rel = rel, σ = noise_sigma.(rs), t = old_t.(rs),
            time = [r.time for r in rs], χ² = [r.chi2_obs_per_dof for r in rs])
end

"OCO ABSCO spectral RMSE from its stored chi squared, using the noise and size of a HITRAN CO2 run `ref`."
function oco_rmse(ref)
    oco = JLD2.load(joinpath(FIG_DATA, "OCO_results.jld2")); n, m = size(ref.rs[1].K)
    j = [argmin(abs.(oco["t"] .- t)) for t in ref.t]
    return sqrt.(oco["co2_chi2"][j] .* (n - m) / n) .* ref.σ
end

"Per-record ratio a.rmse / base (base: values on the records of `b`), matched by measurement time."
function rmse_ratio(a, b, base = b.rmse)
    at = Dict(zip(b.time, base))
    return [x / at[t] for (x, t) in zip(a.rmse, a.time) if haskey(at, t)]
end

if abspath(PROGRAM_FILE) == @__FILE__
    s(x) = @sprintf("%.4g", x)
    c08, c24, k24 = spectral_rmse("ch4_hit08"), spectral_rmse("ch4_hit24"), spectral_rmse("co2_hit24")
    rm_oco = k24 === nothing ? nothing : oco_rmse(k24)
    if k24 !== nothing
        chk = sqrt.(k24.χ² .* (size(k24.rs[1].K, 1) - size(k24.rs[1].K, 2)) / size(k24.rs[1].K, 1)) .* k24.σ
        println("formula check on HIT24 CO2: median ratio formula/direct = ", s(median(chk ./ k24.rmse)))
    end
    h08 = JLD2.load(joinpath(FIG_DATA, "hit08_results.jld2"))
    if c08 !== nothing
        j8 = [argmin(abs.(h08["t"] .- t)) for t in c08.t]
        println("local vs stored HIT08 chi2: median ratio = ", s(median(c08.χ² ./ h08["ch4_chi2"][j8])))
    end

    println("```csv")
    println("experiment,median_rmse_transmittance_units,mean_rmse,median_rmse_percent_of_signal,median_noise_sigma,median_rmse_over_sigma,median_rmse_ratio_to_voigt_baseline")
    for n in hit24_experiments()
        has_results(n) || (println("# not run: ", n); continue)
        r = spectral_rmse(n)
        # baseline: HITRAN 2008 (CH4) or OCO ABSCO from its stored chi squared (CO2)
        base = startswith(n, "ch4") ? (c08 === nothing ? NaN : median(rmse_ratio(r, c08))) :
                                      (k24 === nothing ? NaN : median(rmse_ratio(r, k24, rm_oco)))
        println(join((n, s(median(r.rmse)), s(mean(r.rmse)), s(median(r.rel)), s(median(r.σ)),
                      s(median(r.rmse ./ r.σ)), s(base)), ","))
    end
    rm_oco === nothing ||
        println(join(("co2_oco (from stored chi2)", s(median(rm_oco)), s(mean(rm_oco)), "NaN", s(median(k24.σ)), s(median(rm_oco ./ k24.σ)), "1"), ","))
    println("```")
    if c08 !== nothing && c24 !== nothing
        q = rmse_ratio(c24, c08)
        @printf("CH4: HIT08 RMSE smaller than HIT24 in %d of %d records; median ratio HIT24/HIT08 = %.4f\n", count(>(1), q), length(q), median(q))
    end
    if k24 !== nothing
        @printf("CO2: OCO RMSE smaller than HIT24 in %d of %d records; median ratio HIT24/OCO = %.4f\n", sum(rm_oco .< k24.rmse), length(rm_oco), median(k24.rmse ./ rm_oco))
    end
    if c24 !== nothing
        p, T = tower_pT(c24.t)
        @printf("tower p mean %.2f hPa, T mean %.2f K over the %d ch4_hit24 records\n", mean(p), mean(T), length(c24.t))
    end
end

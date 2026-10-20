#=

Compare the HITRAN 2024 retrievals on 20160921 (configs/hit24.toml, hit24_co2.toml) against
the existing campaign results in fignewtons/timeseries/data:
  CH4: HITRAN 2008 (hit08_results.jld2) vs HITRAN 2024 (ch4_hit24)
  CO2: OCO ABSCO   (OCO_results.jld2)   vs HITRAN 2024 (co2_hit24)
plus the Picarro in situ record (obs.jld2) as an independent reference, and the local re-run
of HITRAN 2008 (ch4_hit08) as a reproduction check. `variants` adds one summary table of
every comparison experiment (line mixing _lm / _sdvlm, fixed p and T _fixpT, ...).

Conventions follow fignewtons/timeseries/process_results.py:
  dry column = air column(p, T, path) - H2O, using p, T, H2O retrieved in the same window
  (fixed-p/T runs: the tower p and T they were retrieved with, metadata["conditions"])
  records kept when 0.2 < chi squared < 3
Old results store t; new results are matched on machine time - OLD_T_OFFSET_S (common.jl),
to within 60 s. Experiments that have not been run are skipped with a note.

Usage: julia --project=../../SpectralFit.jl compare_hit24.jl [ch4] [co2] [variants]
Per-record tables: ../data/runs/hit24_comparison/

=#

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "rmse_hit24.jl"))     # spectral_rmse, oco_rmse, hit24_experiments

const out_dir = mkpath(joinpath(RUNS_DIR, "hit24_comparison"))

vmr(col, p, T, h2o, scale; path = 195017.0) = scale * col / (homogeneous_air_column(p, T, path) - h2o)
good(χ²) = isfinite(χ²) && 0.2 < χ² < 3.0

"New SpectralFits results -> (t, vmr, sigma, χ², p, T), or nothing when not run."
function load_new(name, gas, scale)
    rs = load_run(name)
    rs === nothing && return nothing
    i = label_index(rs[1], gas)
    return map(rs) do r
        c, dry = pT(r), dry_column(r)
        (t = old_t(r), vmr = scale * r.x[i] / dry, σ = scale * sqrt(r.posterior_covariance[i, i]) / dry,
         χ² = r.chi2_obs_per_dof, p = c.pressure_hPa, T = c.temperature_K)
    end
end

"Old campaign results; `gas` = \"ch4\" uses the CH4-window p, T, H2O."
function load_old(file, gas, scale)
    d = JLD2.load(joinpath(FIG_DATA, file))
    p, T, h = gas == "ch4" ? (d["p_ch4"], d["T_ch4"], d["h2o_ch4"]) : (d["p"], d["T"], d["h2o"])
    return [(t = d["t"][k], vmr = vmr(d[gas][k], p[k], T[k], h[k], scale), σ = NaN,
             χ² = d["$(gas)_chi2"][k], p = p[k], T = T[k]) for k in eachindex(d["t"])]
end

"Picarro mole fraction at the DCS record times (t_pic is days since the first DCS record)."
function picarro(gas, t_dcs)
    o = JLD2.load(joinpath(FIG_DATA, "obs.jld2"))
    tp = campaign_t0() .+ 86400 .* o["t_pic"]
    return map(t_dcs) do t
        k = argmin(abs.(tp .- t))
        abs(tp[k] - t) < 450 ? o[gas][k] : NaN
    end
end

function match_records(a, b)
    pairs = Tuple{Int,Int}[]
    for (i, x) in enumerate(a)
        j = argmin(abs.([y.t for y in b] .- x.t))
        abs(b[j].t - x.t) < 60 && push!(pairs, (i, j))
    end
    return pairs
end

csvrow(io, q, vals...) = println(io, q, ",", join([@sprintf("%.4f", v) for v in vals], ","))

function compare(gas, unit, scale, old, new, label_old, label_new; repro = nothing)
    new === nothing && return println("note: $(uppercase(gas)) comparison skipped")
    prs = match_records(new, old)
    keep = [(i, j) for (i, j) in prs if good(new[i].χ²) && good(old[j].χ²)]
    A = [old[j] for (i, j) in keep]; B = [new[i] for (i, j) in keep]
    g(s, k) = [getfield(v, k) for v in s]
    d = g(B, :vmr) .- g(A, :vmr)
    pic = picarro(gas, g(B, :t))
    hp = .!isnan.(pic)

    println("\n$(uppercase(gas)), 20160921: $label_new minus $label_old. ",
            "Matched records $(length(prs)) of $(length(new)); kept after chi squared filter $(length(keep)); with Picarro $(sum(hp)).")
    println("```csv")
    println("quantity,", label_old, ",", label_new, ",difference_new_minus_old")
    csvrow(stdout, "mean $(uppercase(gas)) ($unit)",           mean(g(A, :vmr)),   mean(g(B, :vmr)),   mean(d))
    csvrow(stdout, "median $(uppercase(gas)) ($unit)",         median(g(A, :vmr)), median(g(B, :vmr)), median(d))
    csvrow(stdout, "std across records ($unit)",               std(g(A, :vmr)),    std(g(B, :vmr)),    std(d))
    csvrow(stdout, "median chi squared per dof",               median(g(A, :χ²)),  median(g(B, :χ²)),  median(g(B, :χ²) .- g(A, :χ²)))
    csvrow(stdout, "mean fitted pressure (hPa)",               mean(g(A, :p)),     mean(g(B, :p)),     mean(g(B, :p) .- g(A, :p)))
    csvrow(stdout, "mean fitted temperature (K)",              mean(g(A, :T)),     mean(g(B, :T)),     mean(g(B, :T) .- g(A, :T)))
    if any(hp)
        csvrow(stdout, "mean minus Picarro ($unit)",           mean(g(A, :vmr)[hp] .- pic[hp]), mean(g(B, :vmr)[hp] .- pic[hp]), NaN)
        csvrow(stdout, "rms minus Picarro ($unit)",            sqrt(mean((g(A, :vmr)[hp] .- pic[hp]).^2)), sqrt(mean((g(B, :vmr)[hp] .- pic[hp]).^2)), NaN)
    end
    println("```")
    @printf("relative difference new minus old: mean %+.3f percent, std %.3f percent; time-series correlation r = %.4f; mean posterior sigma (new) = %.3f %s\n",
            100 * mean(d ./ g(A, :vmr)), 100 * std(d ./ g(A, :vmr)), cor(g(A, :vmr), g(B, :vmr)), mean(g(B, :σ)), unit)

    if repro !== nothing
        rp = match_records(repro, old)
        dr = [repro[i].vmr - old[j].vmr for (i, j) in rp if good(repro[i].χ²) && good(old[j].χ²)]
        @printf("reproduction check (local re-run minus stored %s): mean %+.4f %s, max abs %.4f %s over %d records\n",
                label_old, mean(dr), unit, maximum(abs.(dr)), unit, length(dr))
    end

    csv = joinpath(out_dir, "$(gas)_$(lowercase(label_old))_vs_$(lowercase(label_new)).csv")
    open(csv, "w") do io
        println(io, "old_t,vmr_old,vmr_new,picarro,sigma_new,chi2_old,chi2_new,p_old,p_new,T_old,T_new")
        for k in eachindex(A)
            println(io, join((B[k].t, A[k].vmr, B[k].vmr, pic[k], B[k].σ, A[k].χ², B[k].χ², A[k].p, B[k].p, A[k].T, B[k].T), ","))
        end
    end
    println("per-record table: ", csv)
end

"""
One row per experiment: XGAS bias and RMSE vs Picarro, spectral RMSE, chi squared per dof,
and fitted p, T minus tower. Baselines are the stored campaign results (HITRAN 2008 for
CH4, OCO ABSCO for CO2); their spectral RMSE comes from the local HITRAN 2008 re-run and
from the stored OCO chi squared (see rmse_hit24.jl), NaN when those have not been run.
"""
function variant_table()
    ref = load_new("ch4_hit24", "CH4", 1e9)
    ref === nothing && return println("note: variant table needs ch4_hit24, skipped")
    ref_t = [x.t for x in ref]                         # the 94 records of 20160921
    k24 = spectral_rmse("co2_hit24")
    c08 = spectral_rmse("ch4_hit08")
    rows = [("ch4", "stored HITRAN2008 (fit p T)", () -> load_old("hit08_results.jld2", "ch4", 1e9), c08 === nothing ? nothing : c08.rmse);
            [(n[1:3], n, () -> load_new(n, uppercase(n[1:3]), n[1:3] == "ch4" ? 1e9 : 1e6), n) for n in hit24_experiments()];
            ("co2", "stored OCO_ABSCO (fit p T)", () -> load_old("OCO_results.jld2", "co2", 1e6), k24 === nothing ? nothing : oco_rmse(k24))]
    println("```csv")
    println("gas,variant,records,mean_bias_vs_picarro,rmse_vs_picarro,bias_percent,median_spectral_rmse,median_chi2_per_dof,mean_p_minus_tower_hPa,mean_T_minus_tower_K")
    for (gas, label, loader, spec) in rows
        startswith(label, "stored") || has_results(label) || (println("# not run: ", label); continue)
        recs = loader()
        recs = [recs[argmin(abs.([x.t for x in recs] .- t))] for t in ref_t]   # same records for every row
        recs = filter(x -> good(x.χ²), recs)
        v, t = [x.vmr for x in recs], [x.t for x in recs]
        pic = picarro(gas, t); h = .!isnan.(pic)
        pt, Tt = tower_pT(t)
        srm = spec isa AbstractString ? spectral_rmse(spec).rmse : spec === nothing ? [NaN] : spec
        @printf("%s,%s,%d,%.3f,%.3f,%+.3f,%.4g,%.3f,%+.2f,%+.2f\n", uppercase(gas), label, length(recs),
                mean(v[h] .- pic[h]), sqrt(mean((v[h] .- pic[h]).^2)), 100 * mean((v[h] .- pic[h]) ./ pic[h]),
                median(srm), median([x.χ² for x in recs]), mean([x.p for x in recs] .- pt), mean([x.T for x in recs] .- Tt))
    end
    println("```")
    println("units: CH4 ppb, CO2 ppm; spectral RMSE in measured-intensity units (noise sigma about 0.00109)")
end

gases = isempty(ARGS) ? ["ch4", "co2"] : ARGS
"variants" in gases && variant_table()

"ch4" in gases && compare("ch4", "ppb", 1e9, load_old("hit08_results.jld2", "ch4", 1e9), load_new("ch4_hit24", "CH4", 1e9),
                          "HITRAN2008", "HITRAN2024"; repro = load_new("ch4_hit08", "CH4", 1e9))
"co2" in gases && compare("co2", "ppm", 1e6, load_old("OCO_results.jld2", "co2", 1e6), load_new("co2_hit24", "CO2", 1e6),
                          "OCO_ABSCO", "HITRAN2024")

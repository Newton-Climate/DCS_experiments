#=
Spectral-fit RMSE on absorption lines only, for pairs of the HITRAN 2024 comparison
experiments on 20160921 (configs/hit24.toml, hit24_co2.toml). The baseline (continuum)
channels are excluded.

A channel counts as "on a line" when the target gas's modelled absorption depth,
1 - exp(-column * cross_section(p, T)), exceeds a threshold in either retrieval
of the pair. The same channels are then used for both line lists in that record,
so the comparison is like-for-like. Channels masked in the fit (CH4 masks) are
always excluded. A second mask uses total absorption by all gases,
1 - model / baseline, where baseline is the fitted Legendre polynomial.

The target-gas lookup table of each experiment is the one its config names; the depths are
evaluated on the measured wavenumbers (`r.coordinate`), without the instrument interpolation,
at the retrieved p and T (fixed-p/T runs: the conditions they were retrieved with). Pairs
with an experiment that has not been run are skipped with a note.

Usage: julia --project=../../SpectralFit.jl rmse_lines_hit24.jl
=#

include(joinpath(@__DIR__, "common.jl"))

const thresholds = (0.01, 0.05, 0.10)

"Absorber file of `gas` in the config of `experiment` (configs/hit24.toml or hit24_co2.toml)."
function table_file(experiment, gas)
    for f in ("hit24.toml", "hit24_co2.toml"), c in SpectralFits.read_config(joinpath(@__DIR__, "configs", f))
        c["experiment"] == experiment && return c["absorbers"][String(gas)]["file"]
    end
    error("no config for experiment $experiment")
end

"Target-gas absorption depth and total absorption depth for one record."
function depths(r, gas, lut)
    c, g = pT(r), r.coordinate
    n_air = homogeneous_air_column(c.pressure_hPa, c.temperature_K, c.path_length_cm)
    q = r.x.H2O / (n_air - r.x.H2O)                    # dry H2O mole fraction
    τ = getproperty(r.x, gas) .* cross_section(lut, g, GasConditions(c.pressure_hPa, c.temperature_K, q))
    b = r.x.baseline.window_1
    base = legendre_matrix(length(g), length(b)) * b
    return 1 .- exp.(-τ), 1 .- r.y_model ./ base
end

rmse(r, k) = sqrt(mean((r.y .- r.y_model)[k] .^ 2))

function pair(gas, name_a, name_b, label_a, label_b)
    A, B = load_run(name_a), load_run(name_b)
    (A === nothing || B === nothing) && return println("note: pair $name_a / $name_b skipped")
    # records of both runs, matched by measurement time
    common_t = intersect(Set(r.time for r in A), Set(r.time for r in B))
    A = sort([r for r in A if r.time in common_t]; by = r -> r.time)
    B = sort([r for r in B if r.time in common_t]; by = r -> r.time)
    la = load_cross_section_model(LookupTable, table_file(name_a, gas))
    lb = load_cross_section_model(LookupTable, table_file(name_b, gas))
    keep = unmasked(A[1])
    σ = median(noise_sigma.(A))
    rows = String[]
    allch = [a.rmse for a in A], [b.rmse for b in B]
    push!(rows, @sprintf("%s,all unmasked channels,%d,%.3g,%.3g,%.4f,%d of %d", gas, sum(keep),
                         median(allch[1]), median(allch[2]), median(allch[2] ./ allch[1]), sum(allch[1] .< allch[2]), length(A)))
    for (kind, which) in (("target-gas depth", 1), ("total depth", 2)), thr in thresholds
        ra, rb, nch = Float64[], Float64[], Int[]
        for (a, b) in zip(A, B)
            da, db = depths(a, gas, la)[which], depths(b, gas, lb)[which]
            k = keep .& ((da .> thr) .| (db .> thr))
            push!(ra, rmse(a, k)); push!(rb, rmse(b, k)); push!(nch, sum(k))
        end
        push!(rows, @sprintf("%s,%s > %d percent,%d,%.3g,%.3g,%.4f,%d of %d", gas, kind, round(Int, 100thr), round(Int, median(nch)),
                             median(ra), median(rb), median(rb ./ ra), sum(ra .< rb), length(A)))
    end
    println("\n$gas: $label_a ($name_a) vs $label_b ($name_b), median noise sigma = $(@sprintf("%.3g", σ))")
    println("```csv")
    println("gas,channels,median_n_channels,median_rmse_$(label_a),median_rmse_$(label_b),median_ratio_$(label_b)_over_$(label_a),records_where_$(label_a)_lower")
    foreach(println, rows)
    println("```")
end

pair(:CH4, "ch4_hit08", "ch4_hit24", "HITRAN2008", "HITRAN2024")
pair(:CH4, "ch4_hit24", "ch4_hit24_lm", "HITRAN2024", "HITRAN2024_LM")
pair(:CH4, "ch4_hit24", "ch4_hit24_sdvlm", "HITRAN2024", "HITRAN2024_SDVLM")
pair(:CO2, "co2_oco", "co2_hit24", "OCO_ABSCO", "HITRAN2024")

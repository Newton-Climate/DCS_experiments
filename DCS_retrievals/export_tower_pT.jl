#=

Tower pressure and temperature of the 2016 campaign as a conditions table for the batch app
(`[conditions] file = ...` in configs/hit24.toml, the *_fixpT variants).

Reads fignewtons/timeseries/data/obs.jld2 (t_pT, p, T) and writes
data/tower_pT_20160921.csv with columns time,pressure_hPa,temperature_K, time in ISO 8601 UTC.

obs.jld2 stores t_pT as days since the first DCS record, whose time hit08_results.jld2 stores
as t[1]. Those v0.1-era times are one day late (common.jl, OLD_T_OFFSET_S), so the true Unix
time is t[1] + 86400 t_pT - 86400 s. The batch app matches each 15-minute average on its first
sample time, so this is the same matching run_hit24.jl _fixpT did on machine_time - OLD_T_OFFSET_S.
Rows with non-finite p or T are dropped.

Usage: julia --project=../../SpectralFit.jl export_tower_pT.jl

=#

using JLD2, Dates, Printf

const FIG_DATA = joinpath(@__DIR__, "..", "fignewtons", "timeseries", "data")
const OUT_FILE = joinpath(@__DIR__, "..", "data", "tower_pT_20160921.csv")

obs = JLD2.load(joinpath(FIG_DATA, "obs.jld2"))
t0  = JLD2.load(joinpath(FIG_DATA, "hit08_results.jld2"))["t"][1]
unix = t0 .+ 86400 .* obs["t_pT"] .- 86400          # old convention -> true Unix time
p, T = obs["p"], obs["T"]
keep = findall(k -> isfinite(unix[k]) && isfinite(p[k]) && isfinite(T[k]), eachindex(unix))

open(OUT_FILE, "w") do io
    println(io, "# tower p and T, 2016 DCS campaign; from fignewtons/timeseries/data/obs.jld2 by DCS_retrievals/export_tower_pT.jl")
    println(io, "time,pressure_hPa,temperature_K")
    for k in keep
        t = unix2datetime(round(unix[k]; digits = 3))
        println(io, Dates.format(t, dateformat"yyyy-mm-ddTHH:MM:SS.sss"), ",",
                @sprintf("%.4f", p[k]), ",", @sprintf("%.4f", T[k]))
    end
end
println("wrote $(length(keep)) of $(length(unix)) rows to $(normpath(OUT_FILE)): ",
        unix2datetime(unix[keep[1]]), " to ", unix2datetime(unix[keep[end]]))

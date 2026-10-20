#=

Collect the linelist_summary_<species>.csv files written by linelist_errors.jl
(runs/linelist_analysis/, or the directory given) into one CSV, linelist_summary_all.csv in
the same directory, and print markdown tables for LINELIST_PROFILE_SUMMARY.md. Pure Julia,
no packages.

Run:  julia linelist_tables_md.jl [directory] > runs/linelist_analysis/summary_tables.md
      (linelist_results/ holds the v0.2 hand-port outputs, same columns)

=#

using Printf

const DIR = isempty(ARGS) ? joinpath(@__DIR__, "runs", "linelist_analysis") : ARGS[1]

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    return [Dict(zip(header, split(l, ","))) for l in lines[2:end] if !isempty(l)]
end

files = filter(f -> occursin(r"^linelist_summary_(CH4|CO2)(_(profile|isotherm|isobar))*\.csv$", f), readdir(DIR))
rows = reduce(vcat, [read_csv(joinpath(DIR, f)) for f in sort(files)])
header = split(readline(joinpath(DIR, first(files))), ",")
open(joinpath(DIR, "linelist_summary_all.csv"), "w") do io
    println(io, join(header, ","))
    foreach(r -> println(io, join([r[h] for h in header], ",")), rows)
end

num(r, k) = parse(Float64, r[k])
f2(x) = @sprintf("%.2f", x)
f3(x) = @sprintf("%+.3f", x)

println("| Species | Truth | Retrieval | Mode | Column err % | Δcol vs control % | Shape RMS vs control % | Max shape vs control % | Col σ % | DOF total | DOF 0–2 | DOF 2–5 | DOF 5–10 | χ²/dof |")
println("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    println("| ", join([r["species"], r["truth"], r["retrieval"], r["mode"],
        f3(num(r, "col_err_pct")), f3(num(r, "dcol_vs_control_pct")),
        @sprintf("%.3f", num(r, "dshape_rms_vs_control_pct")), @sprintf("%.3f", num(r, "dshape_maxabs_vs_control_pct")),
        @sprintf("%.3f", num(r, "col_sigma_pct")), f2(num(r, "dof_total")), f2(num(r, "dof_0_2km")),
        f2(num(r, "dof_2_5km")), f2(num(r, "dof_5_10km")), @sprintf("%.3f", num(r, "chi2_per_dof"))], " | "), " |")
end

#=

Compare MERLIN vs HITRAN-2008 retrievals of CH4 over the same window,
sample, and inversion setup — only the CH4 line list differs
(configs/merlin.toml: r6_merlin vs r6_hit08_narrow). Exits with a note
when one has not been run.

Usage: julia --project=../../SpectralFit.jl compare_merlin_vs_hit08.jl
=#

include(joinpath(@__DIR__, "common.jl"))

merlin_all = load_run("r6_merlin")
hit08_all  = load_run("r6_hit08_narrow")
(merlin_all === nothing || hit08_all === nothing) && exit()

# records converged (not Failed) in both runs, matched by measurement time
times  = intersect(Set(r.time for r in merlin_all), Set(r.time for r in hit08_all))
merlin = [r for r in merlin_all if r.time in times]
hit08  = [r for r in hit08_all if r.time in times]
N = length(merlin)

println("State on one result: ", join(merlin[1].state_labels, ", "))
println()

extract(r) = (ch4 = r.x.CH4, h2o = r.x.H2O, χ2 = r.chi2_obs_per_dof)

@printf "%3s | %-22s | %-22s | %8s\n" "i" "MERLIN CH4  [molec/cm²]" "HIT08 CH4  [molec/cm²]" "Δ%"
println(repeat("-", 76))
Δpct = Float64[]
for i in 1:N
    m = extract(merlin[i])
    h = extract(hit08[i])
    dp = 100 * (m.ch4 - h.ch4) / h.ch4
    push!(Δpct, dp)
    if i <= 5 || i > N - 5
        @printf "%3d | %20.4e   | %20.4e   | %+7.3f\n" i m.ch4 h.ch4 dp
    elseif i == 6
        println("   ...")
    end
end
println(repeat("-", 76))
@printf "N = %d records\n" N
@printf "MERLIN − HIT08 CH4 column, mean   = %+7.3f %%\n" mean(Δpct)
@printf "MERLIN − HIT08 CH4 column, median = %+7.3f %%\n" median(Δpct)
@printf "MERLIN − HIT08 CH4 column, std    =  %7.3f %%\n" std(Δpct)

# χ² comparison
χm = [extract(r).χ2 for r in merlin]
χh = [extract(r).χ2 for r in hit08]
@printf "\nχ² median: MERLIN = %.3f    HIT08 = %.3f\n" median(χm) median(χh)
@printf "χ² mean:   MERLIN = %.3f    HIT08 = %.3f\n" mean(χm) mean(χh)

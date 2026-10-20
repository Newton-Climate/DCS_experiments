#=

Three-way comparison over the MERLIN R(6) window (6076–6078 cm⁻¹):
  Voigt-HIT08     — 8 lines, Voigt
  Voigt-MERLIN    — 40 lines, Voigt (profile mismatch for MERLIN)
  HT-MERLIN       — 40 lines, Hartmann–Tran + line mixing (physically correct)

All three use identical poly_degree=5, identical state vector, identical
window and data cut. Only the LUT differs (configs/merlin.toml: experiments
r6_hit08_narrow, r6_merlin, r6_merlin_HT). Exits with a note when one has not been run.

Usage: julia --project=../../SpectralFit.jl compare_all.jl
=#

include(joinpath(@__DIR__, "common.jl"))

runs = (voigt_hit08 = load_run("r6_hit08_narrow"), voigt_merlin = load_run("r6_merlin"),
        ht_merlin = load_run("r6_merlin_HT"))
any(isnothing, runs) && exit()

# records converged (not Failed) in all three runs, matched by measurement time
times = intersect((Set(r.time for r in rs) for rs in runs)...)
voigt_hit08, voigt_merlin, ht_merlin = ([r for r in rs if r.time in times] for rs in runs)
N = length(voigt_hit08)

extract(r) = (ch4 = r.x.CH4, h2o = r.x.H2O, χ2 = r.chi2_obs_per_dof)

# Per-record arrays
ch4_h08 = [extract(voigt_hit08[i]).ch4  for i in 1:N]
ch4_vm  = [extract(voigt_merlin[i]).ch4 for i in 1:N]
ch4_hm  = [extract(ht_merlin[i]).ch4    for i in 1:N]
χ_h08   = [extract(voigt_hit08[i]).χ2   for i in 1:N]
χ_vm    = [extract(voigt_merlin[i]).χ2  for i in 1:N]
χ_hm    = [extract(ht_merlin[i]).χ2     for i in 1:N]

println()
println("N = ", N, " records | poly_degree = 5 for all three")
println()
@printf "%-16s  %-14s  %-14s  %-14s\n" "" "χ² median" "χ² mean" "χ² min–max"
println(repeat("-", 66))
@printf "%-16s  %-14.3f  %-14.3f  %.3f – %.3f\n" "Voigt HIT08"  median(χ_h08) mean(χ_h08) minimum(χ_h08) maximum(χ_h08)
@printf "%-16s  %-14.3f  %-14.3f  %.3f – %.3f\n" "Voigt MERLIN" median(χ_vm)  mean(χ_vm)  minimum(χ_vm)  maximum(χ_vm)
@printf "%-16s  %-14.3f  %-14.3f  %.3f – %.3f\n" "HT MERLIN"    median(χ_hm)  mean(χ_hm)  minimum(χ_hm)  maximum(χ_hm)
println()

@printf "%-24s  %-22s  %-22s  %-22s\n" "" "median" "mean" "std"
println(repeat("-", 94))
@printf "%-24s  %20.4e    %20.4e    %20.4e\n" "CH4 col: Voigt HIT08"  median(ch4_h08) mean(ch4_h08) std(ch4_h08)
@printf "%-24s  %20.4e    %20.4e    %20.4e\n" "CH4 col: Voigt MERLIN" median(ch4_vm)  mean(ch4_vm)  std(ch4_vm)
@printf "%-24s  %20.4e    %20.4e    %20.4e\n" "CH4 col: HT MERLIN"    median(ch4_hm)  mean(ch4_hm)  std(ch4_hm)
println()

pct(a, b) = 100 .* (a .- b) ./ b
Δvm_h08 = pct(ch4_vm, ch4_h08)
Δhm_h08 = pct(ch4_hm, ch4_h08)
Δhm_vm  = pct(ch4_hm, ch4_vm)

@printf "Column bias vs Voigt HIT08 [%%]:\n"
@printf "  Voigt MERLIN − HIT08 : mean %+7.3f  median %+7.3f  std %6.3f\n" mean(Δvm_h08) median(Δvm_h08) std(Δvm_h08)
@printf "  HT    MERLIN − HIT08 : mean %+7.3f  median %+7.3f  std %6.3f\n" mean(Δhm_h08) median(Δhm_h08) std(Δhm_h08)
@printf "  HT    MERLIN − Voigt MERLIN : mean %+7.3f  median %+7.3f\n" mean(Δhm_vm) median(Δhm_vm)

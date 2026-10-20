#=

Column averaging kernels and retrieved profiles of the six synthetic profile experiments
(standard, isobaric, isothermal × CO2 window, CH4 window), read from the run folders in runs/.

    julia --project=../../SpectralFit.jl error_analysis.jl [legacy|us_standard] [suffix]

`legacy` (default) reads the legacy_ch4_<mode> and legacy_co2_<mode> runs
(configs/legacy_ch4.toml, configs/legacy_co2.toml); `us_standard` reads us_standard_ch4_<mode>
and us_standard_co2_<mode>. `suffix` is appended to every experiment name, for example
`_tables` for the us_standard table variants. Each panel uses the layer columns and truth
stored with its own run; H2O is taken from the CO2-window runs (v0.1 used the CH4 standard
and the CO2 isobar/isotherm runs).

Outputs (in runs/): error_analysis_<set><suffix>.csv, one row per gas, mode and layer
(gas,window,mode,z_km,column_ak,vmr_retrieved,vmr_true,dof_layer); with Plots in the
environment stack (for example `julia -e 'using Pkg; Pkg.add("Plots")'`)
also column_averaging_kernel_<set><suffix>.pdf and profile_retrieval_vmr_<set><suffix>.pdf.

=#

include("util.jl")

const SET = length(ARGS) >= 1 ? ARGS[1] : "legacy"
const SUFFIX = length(ARGS) >= 2 ? ARGS[2] : ""
SET in ("legacy", "us_standard") || error("set must be legacy or us_standard, got $SET")

const MODES = (profile = "standard", isobar = "isobaric", isotherm = "isothermal")

runs = Dict((window, mode) => load_run("$(SET)_$(lowercase(String(window)))_$(mode)$(SUFFIX)")
            for window in (:CH4, :CO2), mode in keys(MODES))

# (gas, window run it comes from, plot title, mole-fraction scale)
const PANELS = ((:H2O, :CO2, "H2O", 1e2), (:CO2, :CO2, "CO2", 1e6), (:CH4, :CH4, "CH4", 1e9))

## table of everything that is plotted
rows = String[]
for (gas, window, _, scale) in PANELS, mode in keys(MODES)
    r = runs[(window, mode)]
    cak, ret, tru, dfs = column_ak(r, gas), scale .* retrieved_vmr(r, gas), scale .* true_vmr(r, gas), layer_dfs(r, gas)
    for j in eachindex(r.z_km)
        push!(rows, join([String(gas), String(window), String(mode),
                          (@sprintf("%.6g", v) for v in (r.z_km[j], cak[j], ret[j], tru[j], dfs[j]))...], ","))
    end
    @printf "%-3s %-8s DOF %.2f  column error %+.3f %%  column AK 0–2 km %.3f\n" gas mode sum(dfs) 100 * (sum(getproperty(r.result.x, gas)) / sum(getproperty(r.x_true, gas)) - 1) mean(cak[r.z_km .< 2])
end
csv = joinpath(RUNS_DIR, "error_analysis_$(SET)$(SUFFIX).csv")
open(csv, "w") do io
    println(io, "gas,window,mode,z_km,column_ak,vmr_retrieved,vmr_true,dof_layer")
    foreach(l -> println(io, l), rows)
end
println("wrote ", csv)

## figures (only when Plots is available)
const HAS_PLOTS = Base.find_package("Plots") !== nothing
if HAS_PLOTS
    using Plots
else
    println("Plots is not in the environment stack; figures skipped")
end

if HAS_PLOTS
    colors = (profile = :blue, isobar = :red, isotherm = :green)
    titles = (H2O = "H₂O", CO2 = "CO₂", CH4 = "CH₄")
    function panel_plot(value, xlabel_text)
        plots = map(enumerate(PANELS)) do (k, (gas, window, _, scale))
            plt = plot(title = titles[gas], ylabel = k == 1 ? "height (meters)" : "")
            for mode in keys(MODES)
                r = runs[(window, mode)]
                plot!(plt, value(r, gas, scale), 1000 .* r.z_km, lw = 2, color = colors[mode],
                      label = k == 1 ? MODES[mode] : false)
            end
            k == 2 && xlabel!(plt, xlabel_text)
            plt, window, gas, scale
        end
        return plots
    end
    tag = "$(SET)$(SUFFIX)"

    ps = panel_plot((r, gas, _) -> column_ak(r, gas), "Column averaging kernel")
    for (plt, window, _, _) in ps
        r = runs[(window, :profile)]
        plot!(plt, ones(length(r.z_km)), 1000 .* r.z_km, label = false, color = :black, lw = 2, ls = :dot)
    end
    plot(first.(ps)..., layout = (1, 3), link = :y)
    plot!(fontfamily = "serif-roman", legendfont = font("Computer Modern", 7), xrotation = 45)
    savefig(joinpath(RUNS_DIR, "column_averaging_kernel_$(tag).pdf"))

    ps = panel_plot((r, gas, scale) -> scale .* retrieved_vmr(r, gas), "VMR")
    for (plt, window, gas, scale) in ps
        r = runs[(window, :profile)]
        plot!(plt, scale .* true_vmr(r, gas), 1000 .* r.z_km, label = false, color = :black, lw = 2, ls = :dot)
    end
    plot(first.(ps)..., layout = (1, 3), link = :y)
    plot!(fontfamily = "serif-roman", legendfont = font("Computer Modern", 7), xrotation = 45)
    savefig(joinpath(RUNS_DIR, "profile_retrieval_vmr_$(tag).pdf"))
    println("wrote the figures to ", RUNS_DIR)
end

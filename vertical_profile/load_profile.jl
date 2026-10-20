#=

Load the most recent run of one synthetic profile experiment from runs/, for example

    julia --project=../../SpectralFit.jl -i load_profile.jl us_standard_co2_isobar_tables

`run` (from `load_run` in util.jl) holds the run `folder`; `result`, the retrieval (fields of a
SpectralFits RetrievalResult: x, x_prior, K, posterior_covariance, averaging_kernel,
dfs_per_element, chi2_obs_per_dof, y, y_model, ...; `result.x.CO2` is the CO2 layer
columns); `attributes`, everything the ProfileOSSE app saved (experiment, target, mode, seed,
noise_sigma, snr, height_m, pressure_hPa, temperature_K, dry_air_column_cm2,
profile_pressure_hPa, x_true / prior_sigma with their _labels); and the rebuilt `x_true`,
`prior_sigma`, `vcd`, `z_km`, `pressure_hPa`, `temperature_K`.

=#

include("util.jl")

run = load_run(isempty(ARGS) ? "us_standard_co2_profile_tables" : ARGS[1])
(; result, attributes, x_true, vcd) = run
gas = run.target

println(run.folder)
@printf "%s: %s, %d iterations, χ²/dof %.3f, SNR %.0f\n" run.experiment result.status result.n_iterations result.chi2_obs_per_dof attributes["snr"]
@printf "%s DOF %.2f, column error %+.4f %%\n" gas sum(layer_dfs(run, gas)) 100 * (sum(getproperty(result.x, gas)) / sum(getproperty(x_true, gas)) - 1)
println("z_km,vmr_retrieved,vmr_true,dof_layer,column_ak")
for (z, a, b, d, c) in zip(run.z_km, retrieved_vmr(run, gas), true_vmr(run, gas), layer_dfs(run, gas), column_ak(run, gas))
    @printf "%.3f,%.6g,%.6g,%.4f,%.4f\n" z a b d c
end

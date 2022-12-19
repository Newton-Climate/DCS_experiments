using SpectralFits, JLD2, Plots, StatsBase
include("util.jl")

snr_vec = collect(1500:250:3000)
sampling_vec = collect(200:-50:50)
data = Array{InversionResults}(undef, (4, length(snr_vec)))

@load "CO2_200GHz_results.jld2"
data[1,:] = result

@load "CO2_150GHz_results.jld2"
data[2,:] = result


@load "CO2_100GHz_results.jld2"
data[3,:] = result

@load "CO2_50GHz_results.jld2"
data[4,:] = result

co2_ind = 1:30

degrees = map(x->DOF(x, co2_ind), data)
co2 = map(x->1e6 .* x.x["CO2"] ./ vcd, data)
co2_truth = collect(range(410, stop=400, length=30))
dco2 = map(x-> co2_truth - x, co2)
error = map(x -> rmsd(x, co2_truth), co2)


### plot DOF as functio nof SNR
width = 2
p1 = plot(snr_vec, degrees[1,:], label="200 GHz", lw=width)
plot!(snr_vec, degrees[2,:], label="150 GHz", lw=width)
plot!(snr_vec, degrees[3,:], label="100 GHz", lw=width)
plot!(snr_vec, degrees[4,:], label="50 GHz", lw=width)
plot!(xlabel="SNR", ylabel="Degrees of Freedom", title="Retrieval Information Content")
    plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))
savefig("degreesOfFreedom_snr.pdf")

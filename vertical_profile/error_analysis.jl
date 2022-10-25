using SpectralFits, JLD2, LinearAlgebra, Plots, LaTeXStrings
include("util.jl")
@load "CH4_profile.jld2"
vmr_ch4_idealized = ch4_idealized
vmr_h2o_idealized = 1e2*result.x["H2O"] ./ vcd
ch4_idealized = result
vcd_idealized = vcd_idealized

@load "CH4_isobar.jld2"
vmr_h2o_isobar = 1e2*result.x["H2O"] ./ vcd
vmr_ch4_isobar = ch4_isobar
ch4_isobar = result
vcd_isobar = vcd_isobar

@load "CH4_isotherm.jld2"
vmr_h2o_isotherm = 1e2*result.x["H2O"] ./ vcd
vmr_ch4_isotherm = ch4_isotherm
ch4_isotherm = result
vcd_isotherm = vcd_isotherm
vcd = vcd_idealized

@load "CO2_profile.jld2"
vmr_co2_idealized = co2_idealized
co2_idealized = result

@load "CO2_isotherm.jld2"
vmr_co2_isotherm = co2_isotherm
co2_isotherm = result

@load "CO2_isobar.jld2"
vmr_co2_isobar = co2_isobar
co2_isobar = result;

## define some setup params 
num_layers = 35;
n = num_layers
h2o_ind = 1:num_layers
ch4_ind = num_layers+1:2*num_layers
#co2_ind = 2*num_layers+1:3*num_layers
co2_ind = ch4_ind;

# the true VMRs
co2 = 1e6*collect(range(410e-6, stop=400e-6, length=n))
ch4 = 1e9*collect(range(2000e-9, stop=1900e-9, length=n))
h2o = 1e2*collect(range(0.001, stop=0.005, length=n))



# weighting function
H = vcd ./ sum(vcd)

unity = ones(length(p_idealized))

# Averaging Kernals
# column weighted averaging kernal 
cK_h2o = averaging_kernel(ch4_idealized, vcd, h2o_ind)
cK_co2 = averaging_kernel(co2_idealized, vcd, co2_ind)
cK_ch4 = averaging_kernel(ch4_idealized, vcd_isobar, ch4_ind)


cK_ch4_isobar = averaging_kernel(ch4_isobar, vcd_isobar, ch4_ind)
cK_ch4_isotherm = averaging_kernel(ch4_isotherm, vcd_isotherm, ch4_ind)

cK_co2_isobar = averaging_kernel(co2_isobar, vcd_isobar, co2_ind)
cK_co2_isotherm = averaging_kernel(co2_isotherm, vcd_isotherm, co2_ind)

cK_h2o_isobar = averaging_kernel(co2_isobar, vcd_isobar, h2o_ind)
cK_h2o_isotherm = averaging_kernel(co2_isotherm, vcd_isotherm, h2o_ind)

## plot 
p_h2o_ck = plot(cK_h2o, z,  lw=2, color=:blue, label="standard")
plot!(cK_h2o_isobar, z,  lw=2, color=:red, label="isobaric")
plot!(cK_h2o_isotherm, z,  lw=2, color=:green, label="isothermal")
plot!(unity, z, label=:false, color=:black, lw=2, ls=:dot)
plot!(title=L"H_2O", ylabel="height (meters)")


p_ch4_ck = plot(cK_ch4, z, lw=2, label=:false, color=:blue)
plot!(cK_ch4_isobar, z,  lw=2, color=:red, label=:false)
plot!(cK_ch4_isotherm, z,  lw=2, color=:green, label=:false)
plot!(unity, z, label=:false, color=:black, lw=2, ls=:dot)
plot!(title=L"CH_4")

p_co2_ck = plot(cK_co2, z,  lw=2, label=:false, color=:blue)
plot!(cK_co2_isobar, z,  lw=2, color=:red, label=:false)
plot!(cK_co2_isotherm, z,  lw=2, color=:green, label=:false)
plot!(unity, z, label=:false, color=:black, lw=2, ls=:dot)
title!(L"CO_2", xlabel="Column averaging kernel")


plot(p_h2o_ck, p_co2_ck, p_ch4_ck, layout=(1,3), link=:y)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7), xrotation=45)
savefig("column_averaging_kernel.pdf")


p_h2o_ck = plot(vmr_h2o_idealized, z,  lw=2, color=:blue, label="standard")
plot!(vmr_h2o_isobar, z,  lw=2, color=:red, label="isobaric")
plot!(vmr_h2o_isotherm, z,  lw=2, color=:green, label="isothermal")
plot!(h2o, z, label=:false, color=:black, lw=2, ls=:dot)
plot!(title=L"H_2O", ylabel="height (meters)")


p_ch4_ck = plot(vmr_ch4_idealized, z, lw=2, label=:false, color=:blue)
plot!(vmr_ch4_isobar, z,  lw=2, color=:red, label=:false)
plot!(vmr_ch4_isotherm, z,  lw=2, color=:green, label=:false)
plot!(ch4, z, label=:false, color=:black, lw=2, ls=:dot)
plot!(title=L"CH_4")

p_co2_ck = plot(vmr_co2_idealized, z,  lw=2, label=:false, color=:blue)
plot!(vmr_co2_isobar, z,  lw=2, color=:red, label=:false)
plot!(vmr_co2_isotherm, z,  lw=2, color=:green, label=:false)
plot!(co2, z, label=:false, color=:black, lw=2, ls=:dot)
title!(L"CO_2", xlabel="VMR")

plot(p_h2o_ck, p_co2_ck, p_ch4_ck, layout=(1,3), link=:y)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7), xrotation=45)
savefig("profile_retrieval_vmr.pdf")

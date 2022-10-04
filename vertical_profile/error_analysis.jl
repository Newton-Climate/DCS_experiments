using SpectralFits, JLD2, LinearAlgebra, Plots, LaTeXStrings
include("util.jl")
@load "CH4_profile.jld2"
ch4_idealized = result
vcd_idealized = vcd_idealized

@load "CH4_isobar.jld2"
ch4_isobar = result
vcd_isobar = vcd_isobar

@load "CH4_isotherm.jld2"
ch4_isotherm = result
vcd_isotherm = vcd_isotherm

@load "CO2_profile.jld2"
co2_idealized = result

@load "CO2_isotherm.jld2"
co2_isotherm = result

@load "CO2_isobar.jld2"
co2_isobar = result;

## define some setup params 
num_layers = 35;
n = num_layers
h2o_ind = 1:num_layers
ch4_ind = num_layers+1:2*num_layers
#co2_ind = 2*num_layers+1:3*num_layers
co2_ind = ch4_ind;

## calculate gain matrix
# CH4
G1_idealized = calc_gain_matrix(ch4_idealized)
G1_isotherm = calc_gain_matrix(ch4_isotherm)
G1_isobar = calc_gain_matrix(ch4_isobar)

# CO2
G2_idealized = calc_gain_matrix(co2_idealized)
G2_isotherm = calc_gain_matrix(co2_isotherm)
G2_isobar = calc_gain_matrix(co2_isobar)



## calculate Averaging Kernal
# CH4
A1_idealized  = G1_idealized * ch4_idealized.K
A1_isotherm = G1_isotherm * ch4_isotherm.K
A1_isobar = G1_isobar * ch4_isobar.K

# cO2
A2_idealized  = G2_idealized * co2_idealized.K
A2_isotherm = G2_isotherm * co2_isotherm.K
A2_isobar = G2_isobar * co2_isobar.K

## calculate degrees of freedom
# CH4
CH4_DOF_idealized = tr(A1_idealized[ch4_ind, ch4_ind])
CH4_DOF_isobar = tr(A1_isobar[ch4_ind, ch4_ind])
CH4_DOF_isotherm = tr(A1_isotherm[ch4_ind, ch4_ind])
@show CH4_DOF_idealized
@show CH4_DOF_isotherm
@show CH4_DOF_isobar

# CO2
CO2_DOF_idealized = tr(A2_idealized[co2_ind, co2_ind])
CO2_DOF_isotherm = tr(A2_isotherm[co2_ind, co2_ind])
CO2_DOF_isobar = tr(A2_isobar[co2_ind, co2_ind])
@show CO2_DOF_idealized
@show CO2_DOF_isotherm
@show CO2_DOF_isobar

# H2O Degrees of freedom
H2O_DOF_isobar = tr(A1_isobar[h2o_ind, h2o_ind])
H2O_DOF_isotherm = tr(A1_isotherm[h2o_ind, h2o_ind])
H2O_DOF_idealized = tr(A1_idealized[h2o_ind, h2o_ind])
@show H2O_DOF_idealized
@show H2O_DOF_isotherm
@show H2O_DOF_isobar


## get vcd
z = collect(range(0, stop=10000, length=n))
T1(T₀, z) = T₀-6.5e-3*z
p1(z) = 1e3*exp(-z/8.5e3)


p_idealized = p1.(z)
T_idealized = T1.(300, z)
vcd = SpectralFits.make_vcd_profile(p_idealized, T_idealized)

# weighting function
H = vcd ./ sum(vcd)

unity = ones(length(p_idealized))

# Averaging Kernals
# column weighted averaging kernal 
cK_h2o = (H'*A1_idealized[h2o_ind, h2o_ind])' ./ H
cK_co2 = (H'*A2_idealized[co2_ind, co2_ind])' ./ H
cK_ch4 = (H'*A1_idealized[ch4_ind, ch4_ind])' ./ H


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



using SpectralFits
using JLD2, Plots, LinearAlgebra, LaTeXStrings, ColorSchemes

function calc_vcd(p::Real, T::Real, δz::Float64, VMR_H₂O::Real)    
    ρₙ = p*(1-VMR_H₂O) / (r*T)*Nₐ/1.0e4
    vcd = δz*ρₙ
    return vcd
end # function calc_vcd

"""Calculate the vertical column density given pressure, temperature, and layer thickness"""
function calc_vcd(p::Real, T::Real, δz::Float64)
    VMR_H₂O = 0
    ρₙ = p*(1-VMR_H₂O) / (r*T)*Nₐ/1.0e4
    vcd = δz*ρₙ
    return vcd
end #function calc_vcd


@load "ch4_hit20_error.JLD2"
pathlength = 195017.0 # round trip path length in cm DCSA
#p, T = collect(p), collect(T)
#p = collect(400:100:1000)
#T = collect(200:10:300)

# get the CH4 and H2O column amounts
np, nT = size(out)
ch4 = [out[i,j].x["CH4"] for i=1:np,j=1:nT]
ch4_col = copy(ch4)
h2o = [out[i,j].x["H2O"] for i=1:np,j=1:nT]

# get the p and T
p = [out[i,j].x["pressure"] for i=1:np,j=1:nT]
T = [out[i,j].x["temperature"] for i=1:np,j=1:nT]

## true p and T. 
p_true = 400.0:50.0:1000.0
T_true = 200.0:10.0:300.0



### calculate the δp and δT
p_true_grid = ones(size(out))
T_true_grid = ones(size(out))
p_true, T_true = collect(p_true), collect(T_true)

for i=1:length(T_true)
    p_true_grid[:,i] = p_true
end

for i=1:length(p_true)
    T_true_grid[i,:] = T_true
end

dp = p .- p_true_grid
dT = T .- T_true_grid

## calculate the true vcd 
vcd_true = calc_vcd.(p_true_grid, T_true_grid, pathlength)

## calculate the retrieved vcd from retrieved p and T
vcd = calc_vcd.(p, T, pathlength)
vcd_dry = vcd - h2o


## calculate the CH4 VMR
ch4 = 1e9*(ch4 ./ vcd_dry)
d_ch4 = 2000 .- ch4
ch4_true = vcd_true * 2000.0e-9

## relative error in CH4 in percent 
d_ch4_col = 100.0*(ch4_col .- ch4_true) ./ ch4_true

## calcualte H2O errors 
h2o_true = vcd_true*0.01
d_h2o = 100*(h2o .- h2o_true) ./ h2o_true
    
## original instrument noise 
σ = 0.005610022028250306 / sqrt(10)
Se = Diagonal(1/σ^2 * ones(length(out[1,1].grid)))
χ² = zeros((np, nT))


error = 2000.0 .- ch4 

# get the extrema for color limits
clims = extrema([d_ch4_col; dp; dT; d_h2o])
cp = reverse(cgrad(:thermal))
cp2 = :thermal
p1 = heatmap(T_true, p_true, d_ch4_col,
                          xlabel="degrees K", ylabel="HPa", title=L"CH_4 Column Error", c=cp, rever=true)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))


p2 = heatmap(T_true, p_true, dp,
                          xlabel="degrees K", ylabel="HPa", title="Pressure Error", c=cp, rever=true)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7), )

p3 = heatmap(T_true, p_true, dT,
                          xlabel="degrees K", ylabel="HPa", title="Temperature Error", c=cp2)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7), )

p4 = heatmap(T_true, p_true, d_h2o,
                          xlabel="degrees K", ylabel="HPa", title=L"H_2O Error", c=cp, colormap_title="percent error", rever=true)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))

plot(p1, p2, p3, p4, layout=(2,2), link=:all)
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))
savefig("column_retrieval_error_hit20_ch4.pdf") 
# p_all = scatter!([0], [0], zcolor=[NaN], , label="", c=:viridis, colorbar_title="cbar", background_color_subplot=:transparent, markerstrokecolor=:transparent, framestyle=:none, inset=bbox(0.1, 0, 0.6, 0.9, :center, :right), subplot=6)

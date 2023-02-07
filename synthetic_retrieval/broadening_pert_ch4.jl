using Distributed
addprocs(5)
@everywhere begin
using SpectralFits, DistributedData
using vSmartMOM
using Dates, StatsBase, Revise
using Plots, JLD2, LinearAlgebra
end

on_fluo = true 
# define the reetrieval parameters
@everywhere inversion_setup = Dict{String,Any}(
    "poly_degree" => 2,
    "fit_pressure" => true,
    "fit_temperature" => true,
    "use_OCO" => false,
    "use_TCCON" => false,
    "architecture" => CPU(),
"averaging_window" => Minute(15),
    "verbose_mode" => true,
"fit_column" => true,
"linear" => false)

# Just defining the spectral windows for each species
@everywhere ν_grid = 6050:0.005:6120
# Read the DCS DAta 
if on_fluo
    datadir = "/net/fluo/data1/data/NIST/DCS_A/"
else
    datadir = "../data/"
end

data = read_DCS_data(joinpath(datadir, "20160921.h5"))
measurement =  get_measurement(1, data, ν_grid[1], ν_grid[end]) # get 1 measurement 


# Get the HiTran parameters
spec_dir = on_fluo ? "/net/fluo/data1/data/NIST/spectra/" : "../spectra/"


    CH₄ = get_molecule_info("CH4", joinpath(spec_dir, "hit08_12CH4.par"), 6, 1, ν_grid, architecture=CPU())
    CH₄_pert = get_molecule_info("CH4", joinpath(spec_dir, "hit08_12CH4.par"), 6, 1, ν_grid, architecture=CPU())
    H₂O = get_molecule_info("H2O", joinpath(spec_dir, "tccon_2020.par"), 1, 1, ν_grid, architecture=CPU())
spec_true = setup_molecules([H₂O, CH₄])

CH₄_pert.model.hitran.n_air .*= 1.05
spec1 = setup_molecules([H₂O, CH₄_pert])



p_guess, T_guess = 500.0, 250.0
pathlength = 195017.0
measurement.vcd = SpectralFits.calc_vcd(p_guess, T_guess, pathlength) # calc v=the vcd given specified p and TDCSA

# the true p and T for computing psudo-measurements 
p = 400.0:20.0:1000.0
T = 200.0:3.0:300.0
n = length(measurement.grid)
inversion_setup["obs_covariance"] = 1.0*I(n)

# create the array that stores p and T for the maping 
params = Array{Tuple{Float64, Float64}}(undef, (length(p), length(T)))
for i in eachindex(p)
    for j in eachindex(T)
        params[i,j] = (p[i], T[j])
    end
end

# Define a function to do the synthetic retrieval
@everywhere function retrieve(p, T, inversion_setup, measurement, spectra, spec_pert)

            println("T=",T, ", p=", p)
pathlength = 195017.0 # round trip path length in meters DCS
            vcd = SpectralFits.calc_vcd(p, T, pathlength, VMR_H₂O=0.005)
#	    fetch(get_from(1, :spec_true))
#	    fetch(get_from(1, :spec1))



    # true state 
            x_true = StateVector("H2O" => 0.005 * vcd,
                                          "CH4" => 2000e-9 * vcd,
                                          "pressure" => p,
                                          "temperature" => T,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])




            f = generate_forward_model(x_true, measurement, spectra, inversion_setup);
            τ = f(x_true)
            σ = 0.005610022028250306 / sqrt(10)
            ϵ = randn(length(τ)) * σ
            measurement.intensity = τ #.+ ϵ

    # initial guess 
           xₐ = StateVector("H2O" => 0.005 * measurement.vcd,
                                                         "CH4" => 2000e-9 * measurement.vcd,
                  "pressure" => measurement.pressure,
                  "temperature" => measurement.temperature,
                  "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])
 

            f1 = generate_forward_model(xₐ, measurement, spec_pert, inversion_setup);
out = nonlinear_inversion(f1, xₐ, measurement, spec1, inversion_setup)
            println("done")


    return out 
end


out = pmap(x->retrieve(x[1], x[2], inversion_setup, measurement, spec_true, spec1),params)

println("done with all data. saving")
pathlength = 195017.0 # round trip path length in cm DCSA
(np, nT) = size(out)

# load the column amounts 
ch4_col = [out[i,j].x["CH4"] for i=1:np,j=1:nT]    
h2o_col = [out[i,j].x["H2O"] for i=1:np,j=1:nT]

# pressure and temperature 
p_true = p
T_true = T
p = [out[i,j].x["pressure"] for i=1:np,j=1:nT]    
T = [out[i,j].x["temperature"] for i=1:np,j=1:nT]
    
# calculate vcd
vcd = SpectralFits.calc_vcd.(p, T, pathlength)

# VMR of gases 
ch4 = 1e9*ch4_col ./ (vcd - h2o_col)    
h2o = 1e2*h2o_col ./ (vcd - h2o_col)
println("Saving data")    
@save "ch4_broadening_results.jld2" p p_true T T_true ch4_col ch4 h2o_col h2o vcd 
CH₄

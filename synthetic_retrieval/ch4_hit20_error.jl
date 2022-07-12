using Distributed
addprocs(15)
@everywhere begin
using SpectralFits
using vSmartMOM
using Dates, StatsBase, Revise
using Plots, JLD2
end

# define the reetrieval parameters
@everywhere inversion_setup = Dict{String,Any}(
    "poly_degree" => 40,
    "fit_pressure" => true,
    "fit_temperature" => true,
    "use_OCO" => false,
    "use_TCCON" => false,
    "architecture" => CPU(),
"averaging_window" => Minute(15),
    "verbose_mode" => true,
"fit_column" => true)

# Just defining the spectral windows for each species
@everywhere ν_grid = 6050:0.005:6120
# Read the DCS DAta 
@everywhere data = read_DCS_data("/net/fluo/data1/data/NIST/DCS_A/20160921.h5")
@everywhere measurement =  get_measurement(1, data, ν_grid[1], ν_grid[end]) # get 1 measurement 


 # Get the HiTran parameters
@everywhere begin
    datadir = "/net/fluo/data1/data/NIST/spectra/"
    CH₄_08 = get_molecule_info("CH4", joinpath(datadir, "hit08_12CH4.par"), 6, 1, ν_grid)
    CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit20_12CH4.par"), 6, 1, ν_grid)
    H₂O = get_molecule_info("H2O", joinpath(datadir, "tccon_2020.par"), 1, 1, ν_grid)
    molecules_08 = [H₂O, CH₄_08]
molecules = [H₂O, CH₄]
end


p_guess, T_guess = 1000.0, 300.0
pathlength = 195017.0
measurement.vcd = SpectralFits.calc_vcd(p_guess, T_guess, pathlength) # calc v=the vcd given specified p and TDCSA

# the true p and T for computing psudo-measurements 
p = 400.0:20.0:1000.0
T = 200.0:3.0:300.0


# create the array that stores p and T for the maping 
params = Array{Tuple{Float64, Float64}}(undef, (length(p), length(T)))
for i in eachindex(p)
    for j in eachindex(T)
        params[i,j] = (p[i], T[j])
    end
end

# Define a function to do the synthetic retrieval
@everywhere function retrieve(p, T, inversion_setup)

            println("T=",T, ", p=", p)
pathlength = 195017.0 # round trip path length in meters DCS
            vcd = SpectralFits.calc_vcd(p, T, pathlength)
    spec_true = setup_molecules(molecules_08)

    # true state 
            x_true = StateVector("H2O" => 0.01 * vcd,
                                          "CH4" => 2000e-9 * vcd,
                                          "pressure" => p,
                                          "temperature" => T,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])




            f = generate_forward_model(x_true, measurement, spec_true, inversion_setup);
            τ = f(x_true)
            σ = 0.005610022028250306 / sqrt(10)
            ϵ = randn(length(τ)) * σ
            measurement.intensity = τ #.+ ϵ

    # initial guess 
           xₐ = StateVector("H2O" => 0.01 * vcd,
                                                         "CH4" => 2000e-9 * vcd,
                  "pressure" => measurement.pressure,
                  "temperature" => measurement.temperature,
                  "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])
 
            spec1 = setup_molecules(molecules)
            f1 = generate_forward_model(xₐ, measurement, spec1, inversion_setup);
out = nonlinear_inversion(f1, xₐ, measurement, spec1, inversion_setup)
            println("done")


    return out 
end


out = pmap(x->retrieve(x[1], x[2], inversion_setup),params) 
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
@save "ch4_hit20_results.jld2" p p_true T T_true ch4_col ch4 h2o_col h2o vcd 

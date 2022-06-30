using Distributed
addprocs(15)
@everywhere begin
using SpectralFits
using RadiativeTransfer
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
@everywhere ν_CH4 = (6050, 6120)
                        ν_CO2 = (6205, 6255);
ν_HDO = (6310,6380);

# Read the DCS DAta 
#data = read_DCS_data("../../data/DCSA/DCS_A_1/20160926.h5")
@everywhere data = read_DCS_data("/net/fluo/data1/data/NIST/DCS_A/20160921.h5")
@everywhere data = take_time_average(data, δt=Minute(20))
@everywhere measurement =  get_measurement(1, data, ν_CH4[1], ν_CH4[2]) # get 1 measurement 

@everywhere ν_min, ν_max = ν_CH4[1]-1, ν_CH4[2]+1
@everywhere ν_range = 6050:6120

 # Get the HiTran parameter
@everywhere begin
CH₄_20 = get_molecule_info("CH4", "../12CH4_20.par", 6, 1, ν_range)
CH₄_08 = get_molecule_info("CH4", "../12CH4_S_hit08.data", 6, 1, ν_range)
H₂O = get_molecule_info("H2O", "../H2O_20.par", 1, 1, ν_range)
CO₂ = get_molecule_info("CO2", "../12CO2_20.par", 2,1,ν_range)
end

# Calculate the cross-sections and store in dictionary
molecules_08 = [H₂O, CH₄_08]


@everywhere molecules = [H₂O, CH₄_20]







p_guess, T_guess = 1000.0, 300.0
pathlength = 195017.0
measurement.vcd = SpectralFits.calc_vcd(p_guess, T_guess, 195017.0) # calc v=the vcd given specified p and TDCSA





# generate data with hit08
## generate data under new pressure grid 

# the true p and T for computing psudo-measurements 
p = 400.0:50.0:1000.0
T = 200.0:10.0:300.0


params = Array{Tuple{Float64, Float64}}(undef, (length(p), length(T)))
for i in eachindex(p)
    for j in eachindex(T)
        params[i,j] = (p[i], T[j])
    end
end


#out = Array{InversionResults}(undef, (length(p), length(T)))
@everywhere function retrieve(p, T, inversion_setup)
pathlength = 195017.0 # round trip path length in meters DCSApathelngth = 


            println("T=",T, ", p=", p)

            vcd = SpectralFits.calc_vcd(p, T, pathlength)
            spec_true = construct_spectra(molecules, ν_grid=ν_min:0.001:ν_max, p=p, T=T)
            x_true = OrderedDict{Any,Any}("H2O" => 0.01 * vcd,
                                          "CH4" => 2000e-9 * vcd,
                                          "pressure" => p,
                                          "temperature" => T,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])




            f = generate_forward_model(x_true, measurement, spec_true, inversion_setup);
            τ = f(x_true)
            σ = 0.005610022028250306 / sqrt(10)
            ϵ = randn(length(τ)) * σ
            measurement.intensity = τ #i.+ ϵ

           xₐ = OrderedDict{Any,Any}("H2O" => 0.01 * vcd,
                                                         "CH4" => 2000e-9 * vcd,
                  "pressure" => measurement.pressure,
                  "temperature" => measurement.temperature,
                  "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])
 
            spec1 = construct_spectra(molecules, ν_grid=ν_min:0.001:ν_max, p=p, T=T)
            f1 = generate_forward_model(xₐ, measurement, spec1, inversion_setup);
out = nonlinear_inversion(f1, xₐ, measurement, spec1, inversion_setup)
            println("done")


    return out 
end
#@everywhere @sync retrieve 


out = pmap(x->retrieve(x[1], x[2], inversion_setup),params) 
println("done with all data. saving")
@save "ch4_hit20_error.JLD2" out p T


using SpectralFits, vSmartMOM, Dates, DiffResults, ForwardDiff
data = read_DCS_data("../20160921.h5")
data = take_time_average!(data, δt=Minute(15))
measurement = get_measurement(10, data, 6050, 6120)



n = 11 # number of layers
inversion_setup = Dict{String,Any}(
    "poly_degree" => 100,
    "fit_pressure" => false,
    "fit_temperature" => false,
    "use_OCO" => false,
"use_TCCON" => false,
"verbose_mode" => true,
"architecture" => CPU(),
"fit_column" => true)

# Just defining the spectral windows for each species
ν_CH4 = (6050, 6120)
ν_range = ν_CH4[1]:0.002:ν_CH4[2]
ν_min , ν_max = ν_CH4[1]-1, ν_CH4[end]+1

datadir = "../../retrieval/julia/data"

CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit16_12CH4.par"), 6, 1, ν_range, architecture=inversion_setup["architecture"])
H₂O = get_molecule_info("H2O", joinpath(datadir, "hit16_H2O.par"), 1, 1, ν_range, architecture=inversion_setup["architecture"])
CO₂ = get_molecule_info("CO2", joinpath(datadir, "hit20_12CO2.par"), 2,1,ν_range, architecture=inversion_setup["architecture"])

# Calculate the cross-sections and store in dictionary
molecules = [H₂O, CH₄, CO₂]
spec = setup_molecules(molecules)
#spec = construct_spectra(molecules, ν_grid=ν_min:0.003:ν_max, p=1e3, T=290.0)




a = ones(n)
xₐ = OrderedDict{String, Vector{Float64}}("H2O" => 0.01*a,
                                                         "CH4" => 2000e-9*a,
                                                         "CO2" => 400e-6 * a,
#                  "pressure" => measurement.pressure*a,
#                  "temperature" => measurement.temperature*a,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

σ = OrderedDict{String, Vector{Float64}}("H2O" => 0.001*a,
                                                         "CH4" => 20e-9*a,
                                                         "CO2" => 4e-6 * a,
#                  "pressure" => measurement.pressure*a,
#                  "temperature" => measurement.temperature*a,
                  "shape_parameters" => ones(inversion_setup["poly_degree"]))

inversion_setup["σ"] = σ

# just testing the fit itself
z = collect(range(0, stop=10000, length=n))
T(T₀, z) = T₀+6.5e-3*z
p(z) = 1e3*exp(-z/8.5e3)


measurement.pressure = p.(z)
measurement.temperature = T.(300, z)

f = generate_profile_model(xₐ, measurement, spec, inversion_setup);
@time out = f(xₐ)
measurement.intensity = out .+ sqrt(measurement.σ²)
result = profile_inversion(f, xₐ, measurement, spec, inversion_setup)
@save "information_content.jld2" result 

using SpectralFits, vSmartMOM, Dates, DiffResults, ForwardDiff, JLD2, Revise, Statistics


## define spectral window
ν_CH4 = (6050, 6120)
ν_grid = ν_CH4[1]:0.002:ν_CH4[2]

# directory to line-lists or lookup tables 
datadir = "/net/fluo/data1/data/NIST/DCS_A/"

## get a measurement struct 
data = read_DCS_data(joinpath(datadir, "20160921.h5"))
data = take_time_average!(data, δt=Minute(15))
measurement = get_measurement(10, data, ν_grid[1], ν_grid[end])

## some params to customize 
inversion_setup = Dict{String,Any}(
    "poly_degree" => 100,
    "fit_pressure" => false,
    "fit_temperature" => false,
"verbose_mode" => true,
"architecture" => CPU(),
"fit_column" => true)

## get molecular data
datadir = "/net/fluo/data1/data/NIST/spectra/"

CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit20_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "hit20_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "hit20_12CO2.jld2"), ν_grid)

# store molecules and models in dictionary
molecules = [H₂O, CH₄, CO₂]
spec = setup_molecules(molecules)

## customize the atmospheric profile
# species concentrations profile
snr = 2000.0 # Instrument SNR
n = 44 # number of layers
co2 = collect(range(410e-6, stop=400e-6, length=n))
ch4 = collect(range(2000e-9, stop=1800e-9, length=n))
h2o = collect(range(0.001, stop=0.01, length=n))

# define p(z) and T(z) functions
z = collect(range(0, stop=10000, length=n)) # height 
T(T₀, z) = T₀ .- 6.5e-3 .* z
p(z) = 1e3*exp.(-z/8.5e3)

# save custom p and T
measurement.pressure = p.(z)
measurement.temperature = 280.0*ones(n)
vcd = SpectralFits.make_vcd_profile(measurement.pressure, measurement.temperature, vmr_H₂O=h2o)

# true state 
x_true = OrderedDict{String, Vector{Float64}}("H2O" => h2o .* vcd,
                                                         "CH4" => ch4 .* vcd,
                                                         "CO2" => co2 .* vcd,
                                              "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

#              # a priori state vector
xₐ = OrderedDict{String, Vector{Float64}}("H2O" => 0.05 * vcd,
                                                         "CH4" => 1900e-9 * vcd,
                                                         "CO2" => 405e-6 * vcd,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

## define a priori 
a = ones(n)
σ = OrderedDict{String, Vector{Float64}}("H2O" => 0.1*vcd,
                                                         "CH4" => 600e-9*vcd,
                                                         "CO2" => 100e-6 * vcd,
                  "shape_parameters" => ones(inversion_setup["poly_degree"]))

## save in setup dictionary 
inversion_setup["σ"] = σ

### define a customized Sₐ⁻¹
inversion_setup["Sₐ⁻¹"] = SpectralFits.make_prior_error(σ)


# generate synthetic data 
f = generate_profile_model(x_true, measurement, spec, inversion_setup);
@time out = f(x_true)

# save synthetic measurement and noise 
signal = mean(out)
σ = signal/snr
measurement.σ² = σ^2
measurement.intensity = out .+ randn(length(measurement.intensity)) * σ

## retrieval
result = profile_inversion(f, xₐ, measurement, spec, inversion_setup)

## save some variables for later processing 
p0 = measurement.pressure; T0 = measurement.temperature 
ch4_isotherm = 1e9*result.x["CH4"] ./ (vcd - result.x["H2O"])
println(ch4_isotherm)
vcd_isotherm = vcd

@save "CH4_isotherm.jld2" result p0 T0 xₐ vcd_isotherm x_true ch4_isotherm n


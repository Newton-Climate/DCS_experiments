using SpectralFits, vSmartMOM, Dates, DiffResults, ForwardDiff, JLD2, Revise, StatsBase

on_fluo = true; ; # ar e we on Fluo server?

## define spectral window
ν_grid = 6180:0.002:6260.0

# directory to line-lists or lookup tables
if on_fluo
    datadir = "/net/fluo/data1/data/NIST/DCS_A/"
    else
    datadir = "../../retrieval/julia/data"
end

## get a measurement struct 
data = read_DCS_data(joinpath(datadir, "20160921.h5"))
data = take_time_average!(data, δt=Minute(15))
measurement = get_measurement(10, data, ν_grid[1], ν_grid[end])

## some params to customize 
inversion_setup = Dict{String,Any}(
    "poly_degree" => 5,
    "fit_pressure" => false,
    "fit_temperature" => false,
"verbose_mode" => true,
    "architecture" => CPU(),
    "γ" => 50.0,
    "linear" => true,
"fit_column" => true)


## get molecular data
if on_fluo
    datadir = "/net/fluo/data1/data/NIST/spectra/"
end

datadir = "../../retrieval/julia/data"
H₂O = get_molecule_info("H2O", joinpath(datadir, "hit20_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "hit20_12CO2.jld2"), ν_grid)

# store molecules and models in dictionary
molecules = [H₂O, CO₂]
spec = setup_molecules(molecules)

## customize the atmospheric profile
# species concentrations profile
snr = 2000.0 # Instrument SNR
n = 35 # number of layers
co2 = collect(range(410e-6, stop=400e-6, length=n))
ch4 = collect(range(2000e-9, stop=1800e-9, length=n))
h2o = collect(range(0.005, stop=0.001, length=n))

# define p(z) and T(z) functions
z = collect(range(0, stop=10000, length=n)) # height
δz = 1e2*mean(diff(z)) # layer thickness in cm
T(T₀, z) = T₀ .- 6.5e-3 .* z
p(z) = 1e3*exp.(-z/8.5e3)

# save custom p and T
vcd = SpectralFits.calc_vcd.(p.(z), T.(280, z), δz, h2o)
measurement.vcd = vcd;
H = vcd ./ sum(vcd)
measurement.pressure = p.(z)
measurement.temperature = (H' * T.(280, z)) * ones(n);
#vcd = SpectralFits.make_vcd_profile(measurement.pressure, measurement.temperature, vmr_H₂O=h2o)



# true state 
x_true = OrderedDict{String, Vector{Float64}}("H2O" => h2o .* vcd,
                                                         "CO2" => co2 .* vcd,
                                              "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

#              # a priori state vector
xₐ = OrderedDict{String, Vector{Float64}}("H2O" => 0.005 * vcd,
                                                         "CO2" => 405e-6 * vcd,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

## define a priori 
a = ones(n)
σ = OrderedDict{String, Vector{Float64}}("H2O" => 0.01 * vcd,
                                                         "CO2" => 100e-6 * vcd,
                  "shape_parameters" => ones(inversion_setup["poly_degree"]))

## save in setup dictionary 
inversion_setup["σ"] = σ

### define a customized Sₐ⁻¹
# Number of gases:
nGases = 2
# Correlation length scale (in pressure) for the n Gases:
pcorr = [0.01, 50.0]
# Call custom make_prior_error function
include("custom_Sa.jl")
inversion_setup["Sₐ⁻¹"] = make_prior_error(σ, nGases, measurement.pressure, pcorr)



# generate synthetic data 
f = generate_profile_model(x_true, measurement, spec, inversion_setup);
@time out = f(x_true)

# save synthetic measurement and noise 
signal = mean(out)
σ = signal/snr
measurement.σ² = σ^2
measurement.intensity = out .+ (σ .* randn(length(out))) ./ out

## retrieval
result = SpectralFits.adaptive_inversion(f, xₐ, measurement, spec, inversion_setup)

## save some variables for later processing 
p0 = measurement.pressure; T0 = measurement.temperature
co2_isotherm = 1e6*result.x["CO2"] ./ vcd
println(co2_isotherm)


@save "CO2_isotherm.jld2" result p0 T0 xₐ vcd x_true co2_isotherm n

include("util.jl")
num_layers = n;
h2o_ind = 1:num_layers
co2_ind = num_layers+1:2*num_layers

degrees = DOF(result, co2_ind)
degrees = DOF(result, co2_ind)
residual = rmsd(result.model, result.measurement)
error = rmsd(1e6co2, co2_isotherm)
@show degrees 
@show error
@show residual 

true_degrees = DOF_at_prior(f, x_true, measurement.σ², result, co2_ind)
@show true_degrees
prior_degrees = DOF_at_prior(f, xₐ, measurement.σ², result, co2_ind)
@show prior_degrees

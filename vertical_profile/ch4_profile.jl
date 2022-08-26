using SpectralFits, vSmartMOM, Dates, JLD2, Revise, Statistics, Plots
include("util.jl")

## define spectral window
ν_CH4 = (6050, 6120)
ν_grid = ν_CH4[1]:0.002:ν_CH4[2]

# directory to line-lists or lookup tables 
#datadir = "../../retrieval/julia/data"
datadir = "/net/fluo/data1/data/NIST/DCS_A/"

# get a measurement struct 
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
    "γ" => 5.0,
"fit_column" => true)

## get molecular data
datadir = "/net/fluo/data1/data/NIST/spectra/"
CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit08_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "tccon_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "hit20_12CO2.jld2"), ν_grid)

# store molecules and models in dictionary
molecules = [H₂O, CH₄, CO₂]
spec = setup_molecules(molecules)

## customize the atmospheric profile
# species concentrations profile
snr = 2000.0 # Instrument SNR
n = 30 # number of layers
co2 = collect(range(410e-6, stop=400e-6, length=n))
ch4 = collect(range(2000e-9, stop=1900e-9, length=n))
h2o = collect(range(0.001, stop=0.01, length=n))

# define p(z) and T(z) functions
z = collect(range(0, stop=10000, length=n)) # height 
T(T₀, z) = T₀ .- 6.5e-3 .* z
p(z) = 1e3*exp.(-z/8.5e3)
δz = 1e2*mean(diff(z)) # layer thickness in cm

# save custom p and T
measurement.pressure = p.(z)
measurement.temperature = T.(300, z)

# calculate dry vcd 
#vcd = SpectralFits.make_vcd_profile(measurement.pressure, measurement.temperature, vmr_H₂O=h2o)
#vcd = ones(n)
vcd = SpectralFits.calc_vcd.(measurement.pressure, measurement.temperature, 2*δz, h2o)
#vcd = ones(n)
# true state
x_true = StateVector("H2O" => h2o .* vcd,
                     "CH4" => ch4 .* vcd,
                     "CO2" => co2 .* vcd,
                                              "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

#              # a priori state vector
xₐ = StateVector("H2O" => 0.01 * vcd,
                 "CH4" => 1950e-9 * vcd,
                 "CO2" => 405e-6 * vcd,
                                          "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

## define a priori 
a = ones(n)
σ = OrderedDict{String, Vector{Float64}}("H2O" => 0.07*vcd,
                                         "CH4" => 400e-9*vcd,
                                         "CO2" => 100e-6 * vcd,
                  "shape_parameters" => ones(inversion_setup["poly_degree"]))

## save in setup dictionary 
inversion_setup["σ"] = σ

### define a customized Sₐ⁻¹

# Number of gases:
nGases = 3
# Correlation length scale (in pressure) for the n Gases:
pcorr = [0.01, 50.0, 80.0]
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
measurement.intensity = out .+ randn(length(measurement.intensity)) * σ

## retrieval
result = SpectralFits.adaptive_inversion(f, xₐ, measurement, spec, inversion_setup)

## save some variables for later processing 
p0, T0 = measurement.pressure, measurement.temperature
ch4_idealized = 1e9*result.x["CH4"] ./ vcd
println(ch4_idealized)
vcd_idealized = vcd;
@save "CH4_profile.jld2" result p0 T0 vcd xₐ x_true ch4_idealized

num_layers = n;
h2o_ind = 1:num_layers
ch4_ind = num_layers+1:2*num_layers
co2_ind = 2*num_layers+1:3*num_layers

degrees = DOF(result, ch4_ind)
@show degrees 
@show result.χ²

p1 = plot(ch4_idealized, measurement.pressure, yflip=true,lw=2, label="retrieved", color=:red)
plot!(1e9*ch4, measurement.pressure, yflip=true,lw=2, label="truth", color=:blue)
title!("VMR")

cK_ch4 = averaging_kernal(result, measurement.vcd, ch4_ind)
p2 = plot(cK_ch4, measurement.pressure, yflip=true,lw=2, label=:false, color=:red)
title!("Averaging Kernal")

plot(p1, p2, layout=(1,2))

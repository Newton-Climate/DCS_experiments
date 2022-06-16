using Distributed
using vSmartMOM, OrderedCollections, Dates, Statistics
addprocs(10)
@everywhere using SpectralFits
@sync @everywhere fit_spectra 

# get data directories 
data_path = "/net/fluo/data1/data/NIST/DCS_A/"
output_path = joinpath(data_path, "TCCON_results/")
datadir = "/net/fluo/data1/data/NIST/spectra/"

# define the retrieval parameters 
inversion_setup = Dict{String,Any}(
    "poly_degree" => 100,
    "fit_column" => true,
    "fit_pressure" => true,
    "fit_temperature" => true,
    "architecture" => CPU(),
    "averaging_window" => Minute(15),
    "verbose_mode" => true)

# Just defining the spectral windows for each species
ν_CH4 = (6050, 6120)
ν_CO2 = (6180, 6260);

# Read the DCS DAta 
data = read_DCS_data(joinpath(data_path, "20160925.h5"))
data = take_time_average!(data, δt=Minute(15))
measurement =  get_measurement(1, data, ν_CO2[1], ν_CO2[2]) # get 1 measurement 

ν_min, ν_max = ν_CO2[1]-1, ν_CO2[2]+1
ν_range = 6000:6400
ν_grid = 6049.0:0.003:6263.0

# Get the HiTran parameters

CH₄ = get_molecule_info("CH4", joinpath(datadir, "tccon_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "tccon_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "tccon_12CO2.jld2"), ν_grid)
HDO = get_molecule_info("HDO", joinpath(datadir, "tccon_HDO.jld2"), ν_grid)
molecules = [H₂O, CH₄, CO₂, HDO]

# define the initial guess 
xₐ = OrderedDict{String,Union{Float64, Vector{Float64}}}("H2O" => 0.01 * measurement.vcd,
    "CH4" => 1900e-9 * measurement.vcd,
                  "CO2" => 400e-6 * measurement.vcd,
                  "HDO" => 1e-4 * measurement.vcd,
                  "pressure" => mean(data.pressure),
                  "temperature" => mean(data.temperature),
                  "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

									 spectral_windows = [ν_CH4, ν_CO2]
process_all_files(xₐ, molecules, inversion_setup, spectral_windows, data_path=data_path, out_path=output_path)

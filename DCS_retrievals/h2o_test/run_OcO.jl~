using Distributed
using vSmartMOM, OrderedCollections, Dates, Statistics
addprocs(15)
using SpectralFits


# get data directories 
data_path = "/net/fluo/data1/data/NIST/DCS_A/"
#data_path = "/Users/newtn/projects/FreqComb/retrieval/julia/data"
output_path = joinpath(data_path, "h2o_tests", "TCCON_results")
#output_path = pwd()
datadir = data_path


# define the retrieval parameters 
inversion_setup = RetrievalSetup(
    "poly_degree" => 100,
    "fit_column" => true,
    "fit_pressure" => true,
    "fit_temperature" => true,
    "architecture" => CPU(),
    "averaging_window" => Minute(15),
    "verbose_mode" => true)

## Just defining the spectral windows for each species
ν_CH4 = (6050, 6120)
ν_CO2 = (6180, 6260);
spectral_windows = [ν_CO2]
ν_grid = 6179.0:0.003:6263.0

## Read the DCS DAta 
data = read_DCS_data(joinpath(data_path, "20160921.h5"))
data = take_time_average!(data, δt=Minute(15))
measurement =  get_measurement(1, data, ν_CO2[1], ν_CO2[2]) # get 1 measurement 

## Get the HiTran parameters
datadir = "/net/fluo/data1/data/NIST/spectra/"
CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit08_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "tccon_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "OCO_12CO2_weakband.jld2"), ν_grid)
HDO = get_molecule_info("HDO", joinpath(datadir, "tccon_HDO.jld2"), ν_grid)
molecules = [H₂O, CH₄, CO₂, HDO]


# define the initial guess 
xₐ = StateVector("H2O" => 0.01 * measurement.vcd,
    "CH4" => 1900e-9 * measurement.vcd,
                  "CO2" => 400e-6 * measurement.vcd,
                  "HDO" => 1e-4 * measurement.vcd,
                  "pressure" => mean(data.pressure),
                  "temperature" => mean(data.temperature),
                  "shape_parameters" => [maximum(measurement.intensity); zeros(inversion_setup["poly_degree"]-1)])

									 
process_all_files(xₐ, molecules, inversion_setup, spectral_windows, data_path=data_path, out_path=output_path)

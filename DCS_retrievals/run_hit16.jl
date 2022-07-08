using Distributed
using vSmartMOM, OrderedCollections, Dates, Statistics
#addprocs(10)
@everywhere using SpectralFits
@sync @everywhere fit_spectra 

## get data directories 
data_path = "/net/fluo/data1/data/NIST/DCS_A/"
output_path = joinpath(data_path, "hit16_results/")

## define lines to filter out
masked_range_ch4 = [6109.211 6109.411; 6116.137 6116.464; 6119.62 6120; 6088.826 6089.246; 6087.358 6087.652; 6076.742 6078.55; 6063.577 6064.077; 6058.666 6059.153; 6109.885 6110.071 ];

## define the retrieval parameters 
inversion_setup = RetrievalSetup(
    "poly_degree" => 100,
    "fit_column" => true,
    "fit_pressure" => true,
    "fit_temperature" => true,
    "architecture" => CPU(),
    "averaging_window" => Minute(15),
    "verbose_mode" => true,
     "masked_windows" => masked_range_ch4)

## Just defining the spectral windows for each species
ν_CH4 = (6050, 6120)
ν_CO2 = (6180, 6260);
spectral_windows = [ν_CH4, ν_CO2]
ν_grid = 6049.0:0.003:6263.0 # grid to pass into retrieval over bowth CO2 and CH4 ranges 



## Read the DCS DAta 
data = read_DCS_data(joinpath(data_path, "20160925.h5"))
data = take_time_average!(data, δt=setup["averaging_window"])
measurement =  get_measurement(1, data, ν_CO2[1], ν_CO2[2]) # get 1 measurement 


## Get the HiTran parameters
datadir = "/net/fluo/data1/data/NIST/spectra/"
CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit16_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "tccon_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(datadir, "hit16_12CO2.jld2"), ν_grid)
HDO = get_molecule_info("HDO", joinpath(datadir, "tccon_HDO.jld2"), ν_grid)
molecules = [H₂O, CH₄, CO₂, HDO]

## define the initial guess 
xₐ = StateVector("H2O" => 0.01 * measurement.vcd,
    "CH4" => 1900e-9 * measurement.vcd,
                  "CO2" => 400e-6 * measurement.vcd,
                  "HDO" => 1e-4 * measurement.vcd,
                  "pressure" => mean(data.pressure),
                  "temperature" => mean(data.temperature),
                  "shape_parameters" => prior_shape_params(data, inversion_setup))

## run the retrievals over all files 
process_all_files(xₐ, molecules, inversion_setup, spectral_windows, data_path=data_path, out_path=output_path)

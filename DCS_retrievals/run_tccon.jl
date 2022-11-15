using Distributed
using vSmartMOM, OrderedCollections, Dates, Statistics
addprocs(10)
@everywhere using SpectralFits
@sync @everywhere fit_spectra 

## get data directories 
on_fluo = true
if on_fluo
    data_dir = "/net/fluo/data1/data/NIST/DCS_A/"
    spectra_dir = "/net/fluo/data1/data/NIST/spectra/"
else
    data_dir = "../data/"
    spectra_dir = "../spectra/"
end

output_path = joinpath(data_dir, "TCCON_results/")

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
    "masked_windows" => masked_range_ch4,
"linear" => false)

## Just defining the spectral windows for each species
ν_CH4 = (6050, 6115)
ν_CO2 = (6180, 6260);
spectral_windows = [ν_CH4, ν_CO2]
ν_grid = 6049.0:0.003:6263.0 # grid to pass into retrieval over bowth CO2 and CH4 ranges 


## Get the HiTran parameters
#
CH₄ = get_molecule_info("CH4", joinpath(spectra_dir, "hit08_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(spectra_dir, "tccon_H2O.jld2"), ν_grid)
CO₂ = get_molecule_info("CO2", joinpath(spectra_dir, "hit20_12CO2.jld2"), ν_grid)
HDO = get_molecule_info("HDO", joinpath(spectra_dir, "tccon_HDO.jld2"), ν_grid)
molecules = [H₂O, CH₄, CO₂, HDO]

# extract a dataset to get the maximum of the intensity for prior on shape params
data = read_DCS_data(joinpath(data_dir, "20160921.h5"))
vcd = SpectralFits.calc_vcd(850.0, 285.0,data.pathlength) 

## define the initial guess 
xₐ = StateVector("H2O" => 0.01 * vcd,
                  "CH4" => 1900e-9 * vcd,
                  "CO2" => 400e-6 * vcd,
                  "HDO" => 1e-4 * vcd,
                  "pressure" => 850.0,
                  "temperature" => 285.0,
                  "shape_parameters" => prior_shape_params(data, inversion_setup))

## run the retrievals over all files 
process_all_files(xₐ, molecules, inversion_setup, spectral_windows, data_path=data_dir, out_path=output_path)

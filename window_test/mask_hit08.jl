using vSmartMOM, Distributed, Revise
addprocs(10)
@everywhere using SpectralFits
using Dates, Statistics
using JLD2
@sync @everywhere fit_spectra

"""
Runs the script to mask out CH₄ lines 
Using Hitran 08 line-list
"""


## The windows to mask out 
masked_range_ch4 = [6109.211 6109.411; 6116.137 6116.464; 6119.62 6120; 6088.826 6089.246; 6087.358 6087.652; 6076.742 6078.55; 6063.577 6064.077; 6058.666 6059.153; 6109.885 6110.071 ];

## define the retrieval parameters
setup = RetrievalSetup(
"poly_degree" => 100,
"fit_pressure" => true,
"fit_temperature" => true,
"masked_windows" => masked_range_ch4,
"architecture" => CPU(),
"averaging_window" => Minute(15),
"verbose_mode" => true,
"fit_column" => true)



## defining the spectral windows for the retrieval and spectral calculation 
ν_CH4 = (6050, 6120)
ν_grid = ν_CH4[1]:0.005:ν_CH4[end]

## Read the DCS DAta
# datafile locations 
#datadir = "/Users/newtn/projects/FreqComb/retrieval/julia/data"
datadir = "/net/fluo/data1/data/NIST/DCS_A/"
data = read_DCS_data(joinpath(datadir, "20160921.h5"))
data = take_time_average!(data, δt=setup["averaging_window"])
measurement =  get_measurement(1, data, ν_CH4[1], ν_CH4[end]) # get 1 measurement

## Get the HiTran parameters
# where are the spectra stored?
# can use either .par files or jld2 lookup tables
#datadir = "/Users/newtn/projects/FreqComb/retrieval/julia/data"
datadir = "/net/fluo/data1/data/NIST/spectra"
CH₄ = get_molecule_info("CH4", joinpath(datadir, "hit08_12CH4.jld2"), ν_grid)
H₂O = get_molecule_info("H2O", joinpath(datadir, "TCCON_H2O.jld2"), ν_grid)
molecules = [H₂O, CH₄]

## define the initial guess
vcd = measurement.vcd
xₐ = StateVector("H2O" => 0.01 * vcd,
    "CH4" => 2000e-9 * vcd,
                  "pressure" => mean(data.pressure),
                  "temperature" => mean(data.temperature),
                 "shape_parameters" => prior_shape_params(data, setup))

## run the retrieval
# for only one run
#f = generate_forward_model(xₐ, measurement, molecules, setup);
#results = nonlinear_inversion(f, xₐ, measurement, molecules, _setup)

# for the entire dataset
spectral_windows = [ν_CH4]
results = run_inversion(xₐ, data, molecules, setup, spectral_windows)
@save "masked_hit08_lines.jld2" results

using JLD2, SpectralFits, Statistics 

function average_spectra(filename, measurement_range; outfilename="averaged_spectra.jld2")
    ds = jldopen(filename, "r")
    results = ds["results"]
    idx = measurement_range[end]
    
   CH4_model, CO2_model = mean([results[i,1].model for i=1:idx]), mean([results[i,2].model for i=1:idx])
CH4_obs, CO2_obs = mean([results[i,1].measurement for i=1:idx]), mean([results[i,2].measurement for i=1:idx]) 
    
    
    println("saving fits to file ", outfilename)
    @save outfilename CO2_model CO2_obs CH4_model CH4_obs
    return CO2_model, CO2_obs, CH4_model, CH4_obs
end


files = ["masked_hit08.jld2", "masked_hit16.jld2", "masked_hit20.jld2", "masked_TCCON.jld2"]
for file in files
    average_spectra(joinpath("../", dir[i], filename), 1:10, outfilename=dir[i]*".jld2")
end


function State(params::Pair{String, UT}...) where UT
    x = OrderedDict{String, UT}(params[i] for i=1:length(params))
    return x

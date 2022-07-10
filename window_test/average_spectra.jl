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

files = ["masked_hit08_lines.jld2", "masked_hit16_lines.jld2", "masked_hit20_lines.jld2", "masked_TCCON_lines.jld2"]
for file in files
    average_spectra(file, 1:90, outfilename=file*"_spectra.jld2")
end

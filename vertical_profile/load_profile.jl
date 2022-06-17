using JLD2, SpectralFits

@load "vertical_profile.JLD2"

# variables are:

# p T vcd

# ch4_results co2_results
# those contain the inversion results
# try fieldnames(InversionResults) for the fields
# x_true xₐ are the state vectors 

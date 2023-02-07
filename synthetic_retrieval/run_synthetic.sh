#!/bin/bash

julia ch4_hit16_error.jl;
julia ch4_hit20_error.jl;
julia ch4_TCCON_error.jl;
julia co2_hit16_error.jl;
julia co2_hit20_error.jl;
julia co2_TCCON_error.jl;
julia broadening_pert_ch4.jl;
echo "done with all runs"

#!/bin/bash
# Run the synthetic vertical-profile experiments (ProfileOSSE configs in configs/), then the
# analysis. Run from this directory. Pass config names to run only those, e.g.
#   ./run_profiles.sh us_standard_co2 legacy_co2
# Single variants: julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/<name>.toml <variant> ...
set -e
cd "$(dirname "$0")"
pkg=../../SpectralFit.jl
julia=(julia --project="$pkg")
configs=("$@")
[ ${#configs[@]} -eq 0 ] && configs=(us_standard_ch4 us_standard_co2 legacy_ch4 legacy_co2 ch4_profile2 crosslist_ch4 crosslist_co2)
for c in "${configs[@]}"; do
    "${julia[@]}" "$pkg/bin/run_config.jl" "configs/$c.toml"
done
"${julia[@]}" error_analysis.jl legacy
"${julia[@]}" error_analysis.jl us_standard
"${julia[@]}" linelist_errors.jl

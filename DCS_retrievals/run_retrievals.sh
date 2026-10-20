#!/bin/bash
# Campaign line-list retrievals in the background, one SpectralFits batch run per config
# (each starts its own workers, n_workers in the config). Logs go to <config>.out here.
#
# Usage: ./run_retrievals.sh [config ...]     (default: hit08 hit16 hit20 tccon)
# e.g.   ./run_retrievals.sh hit24 merlin

cd "$(dirname "$0")"
sf=../../SpectralFit.jl
configs=${@:-hit08 hit16 hit20 tccon}
for c in $configs; do
    nohup julia --project=$sf $sf/bin/run_config.jl configs/${c%.toml}.toml > ${c%.toml}.out 2>&1 &
done

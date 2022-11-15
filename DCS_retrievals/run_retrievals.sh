#!/bin/bash

nohup julia run_hit08.jl > hit08.out &
nohup julia run_hit16.jl > hit16.out &
nohup julia run_hit20.jl > hit20.out &
nohup julia run_tccon.jl > tccon.out &

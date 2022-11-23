#!/usr/bin/env julia

#
# This script is used to start analytical continuation simulations.
# It will launch only 8 processes.
#
# Usage:
#
#     $ acrun.jl ac.toml
#

using Distributed
addprocs(8)

@everywhere push!(LOAD_PATH, ENV["ACFLOW_HOME"])

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

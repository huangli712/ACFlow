#!/usr/bin/env julia

#
# This script is used to start analytical continuation simulations.
# It will launch only 1 process.
#
# Usage:
#
#     $ acrun.jl ac.toml
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

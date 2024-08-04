#!/usr/bin/env julia
push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src")
#
# This script is used to start analytic continuation simulations.
# It will launch only 1 process.
#
# Usage:
#
#     $ acrun.jl ac.toml
#

using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

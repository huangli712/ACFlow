#!/usr/bin/env julia

#
# This script is used to start analytic continuation simulations.
# It will launch only 1 process.
#
# Usage:
#
#     $ acrun.jl ac.toml
#

#push!(LOAD_PATH, "/home/lihuang/yun/devel/acflow/src")
#push!(LOAD_PATH, "/home/lihuang/Downloads/ACFlow-1.6.2/src")
push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src")
using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

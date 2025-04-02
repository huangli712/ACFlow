#!/usr/bin/env julia

#
# This julia script is used to perform analytic continuation calculations
# sequentially. It suits for the MaxEnt, BarRat, and NevanAC solvers. It
# will launch only 1 process.
#
# Usage:
#
#     $ acrun.jl ac.toml
#

haskey(ENV,"ACFLOW_PATH") && pushfirst!(LOAD_PATH, ENV["ACFLOW_PATH"])

using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

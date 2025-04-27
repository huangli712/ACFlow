#!/usr/bin/env julia

#
# This julia script is used to perform analytic continuation calculations
# parallely. It suits the stochastic analytic continuation methods, such
# as the StochAC, StochSK, StochOM, and StochPX solvers. By default, it
# will launch only 8 processes. Please fix line (18) to change the number
# of used processors.
#
# Usage:
#
#     $ acprun.jl ac.toml
#

haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Distributed
addprocs(8)

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

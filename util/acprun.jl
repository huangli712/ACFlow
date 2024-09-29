#!/usr/bin/env julia

#
# This julia script is used to perform analytic continuation calculations
# parallely. It suits stochastic analytic continuation methods, such as
# the StochAC, StochSK, StochOM, and StochPX solver. By default it will
# launch only 8 processes. Please fix line (16) to change the number of
# used processors.
#
# Usage:
#
#     $ acprun.jl ac.toml
#

using Distributed
addprocs(8)

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

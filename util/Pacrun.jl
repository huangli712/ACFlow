#!/usr/bin/env julia

#
# This script is used to start analytic continuation simulations.
# By default it will launch only 8 processes. Please fix line (14)
# to change the number of used processors.
#
# Usage:
#
#     $ Pacrun.jl ac.toml
#

using Distributed
addprocs(8)

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

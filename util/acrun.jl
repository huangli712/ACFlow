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

using ACFlow

welcome()
overview()
read_param()
solve(read_data())
goodbye()

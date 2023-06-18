#!/usr/bin/env julia

using ACFlow

welcome()

# For MaxEnt solver

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.inp",
    "mtype"  => "gauss",
    "mesh"   => "lorentz",
    "ngrid"  => 100,
    "wmax"   => 12.0,
    "wmin"   => -12.0,
    "beta"   => 5.0,
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
    "nalph"  => 15,
    "alpha"  => 1e11,
)
#
setup_param(B, S)

# Call the solver
mesh, Aout, Gout = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.mem.data", force = true)
cp("Gout.data", "Gout.mem.data", force = true)
cp("repr.data", "repr.mem.data", force = true)

# For StochOM solver

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "solver" => "StochOM"
)
#
# For [StochOM] block
# See types.jl/_PStochOM for default setup
S = Dict{String,Any}(
    "ntry"   => 100000
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aout, Gout = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.som.data", force = true)
cp("Gout.data", "Gout.som.data", force = true)
cp("repr.data", "repr.som.data", force = true)

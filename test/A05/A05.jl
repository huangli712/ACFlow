#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

welcome()



# Setup parameters
C = Dict{String,Any}(
    "finput" => "sigma.inp",
    "mtype"  => "gauss",
    "mesh"   => "tangent",
    "ngrid"  => 300,
    "wmax"   => 30.0,
    "wmin"   => -30.0,
    "beta"   => 38.0,
)
#
S = Dict{String,Any}(
    "nalph"  => 13,
    "alpha"  => 1e12,
    "blur"   => 0.3,
)
#
setup_param(C, S)

# Call the solver
Aout, Gout = solve(read_data())

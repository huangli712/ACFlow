#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using ACFlow

welcome()

# For MaxEnt solver

# Setup parameters
C = Dict{String,Any}(
    "finput" => "sigma.data",
    "mtype"  => "gauss",
    "mesh"   => "tangent",
    "ngrid"  => 100,
    "wmax"   => 20.0,
    "wmin"   => -20.0,
    "beta"   => 10.0,
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
mesh, Aout, Gout = solve(read_data())
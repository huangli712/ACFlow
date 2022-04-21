#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

C = Dict{String,Any}(
    "finput" => "sigma.inp",
    "mtype"  => "gauss",
    "mesh"   => "tangent"
    "ngrid"  => 300,
    "wmax"   => 30.0,
    "wmin"   => -30.0,
    "beta"   => 38.0,
)

S = Dict{String,Any}(
    "nalph"  => 13,
    "alpha"  => 1e12,
    "blur"   => 0.3,
)

welcome()
setup_param(C, S)
Aout, Gout = solve(read_data())
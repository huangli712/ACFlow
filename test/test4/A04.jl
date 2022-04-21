#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

C = Dict{String,Any}(
    "finput" => "green.inp",
    "ngrid"  => 100,
    "wmax"   => 12.0,
    "wmin"   => -12.0,
    "beta"   => 5.0,
)

S = Dict{String,Any}(
    "nalph"  => 15,
    "alpha"  => 1e11,
)

welcome()
setup_param(C, S)
Aout, Gout = solve(read_data())
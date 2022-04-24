#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

C = Dict{String,Any}(
    "finput" => "green.11.data",
    "mtype"  => "gauss",
    "ngrid"  => 20,
    "nmesh"  => 400,
    "wmax"   => 4.0,
    "wmin"   => -4.0,
    "beta"   => 40.0,
)

S = Dict{String,Any}(
    "nalph"  => 28,
    "alpha"  => 1e18,
)

welcome()
setup_param(C, S)
Aout11, Gout11 = solve(read_data())

cp("Aout.data", "Aout.11.data", force = true)
cp("Gout.data", "Gout.11.data", force = true)
cp("repr.data", "repr.11.data", force = true)

C = Dict{String,Any}(
    "finput" => "green.22.data",
)

S = Dict{String,Any}(
)

setup_param(C, S, false)
Aout22, Gout22 = solve(read_data())

cp("Aout.data", "Aout.22.data", force = true)
cp("Gout.data", "Gout.22.data", force = true)
cp("repr.data", "repr.22.data", force = true)
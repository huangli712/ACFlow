#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

C = Dict{String,Any}(
    "finput" => "sigma.inp",
    "ngrid"  => 100,
    "wmax"   => 20.0,
    "wmin"   => -20.0,
    "beta"   => 5.0,
)

S = Dict{String,Any}(
    "nalph"  => 15,
    "alpha"  => 1e12,
)

welcome()
setup_param(C, S)
Aout, Gout = solve(read_data())

cp("Aout.data", "Aout.data.mem", force = true)
cp("Gout.data", "Gout.data.mem", force = true)
cp("repr.data", "repr.data.mem", force = true)

C = Dict{String,Any}(
    "solver" => "StochOM"
)

S = Dict{String,Any}(
)

setup_param(C, S, false)
Aout, Gout = solve(read_data())
cp("Aout.data", "Aout.data.som", force = true)
cp("Gout.data", "Gout.data.som", force = true)
cp("repr.data", "repr.data.som", force = true)

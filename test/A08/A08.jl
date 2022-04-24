#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

#"finput"  => "green.data",
#"solver"  => "MaxEnt",
#"ktype"   => "fermi",
#"mtype"   => "flat",
#"grid"    => "ffreq",
#"mesh"    => "linear",
#"ngrid"   => 10,
#"nmesh"   => 501,
#"wmax"    => 5.0,
#"wmin"    => -5.0,
#"beta"    => 10.0,
#"offdiag" => false,

C = Dict{String,Any}(
    "finput" => "green11.data",
    "mtype"  => "gauss",
    "ngrid"  => 20,
    "nmesh"  => 400,
    "wmax"   => 4.0,
    "wmin"   => -4.0,
    "beta"   => 40.0,
)

#"method"  => "chi2kink",
#"nalph"   => 12,
#"alpha"   => 1e9,
#"ratio"   => 10.0,
#"blur"    => -1.0,

S = Dict{String,Any}(
    "nalph"  => 28,
    "alpha"  => 1e18,
)

welcome()
setup_param(C, S)
Aout, Gout = solve(read_data())

cp("Aout.data", "Aout.11.data", force = true)
cp("Gout.data", "Gout.11.data", force = true)
cp("repr.data", "repr.11.data", force = true)

error()

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
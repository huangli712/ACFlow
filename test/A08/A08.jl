#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Printf
using ACFlow

welcome()

# For MaxEnt solver

# For diagonal elements: green.11.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.11.data",
    "mtype"  => "gauss",
    "ngrid"  => 20,
    "nmesh"  => 400,
    "wmax"   => 4.0,
    "wmin"   => -4.0,
    "beta"   => 40.0,
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
    "nalph"  => 28,
    "alpha"  => 1e18,
)
#
setup_param(B, S)

# Call the solver
mesh, Aout11, Gout11 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.11.data", force = true)
cp("Gout.data", "Gout.11.data", force = true)
cp("repr.data", "repr.11.data", force = true)

# For diagonal elements: green.22.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.22.data",
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aout22, Gout22 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.22.data", force = true)
cp("Gout.data", "Gout.22.data", force = true)
cp("repr.data", "repr.22.data", force = true)

# For diagonal elements: green.diff.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "gdiff.data",
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aoutdiff, Goutdiff = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.diff.data", force = true)
cp("Gout.data", "Gout.diff.data", force = true)
cp("repr.data", "repr.diff.data", force = true)

# For non-diagonal elements: green.12.data

# Generate model function at first
model_offdiag = sqrt.(Aout11 .* Aout22)
#
open("model.inp", "w") do fout
    for i in eachindex(mesh)
        @printf(fout, "%20.16f %20.16f\n", mesh[i], model_offdiag[i])
    end
end

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.12.data",
    "mtype"  => "file",
    "offdiag"=> true,
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
    "nalph"  => 30,
    "alpha"  => 1e15,
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aout12, Gout12 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.12.data", force = true)
cp("Gout.data", "Gout.12.data", force = true)
cp("repr.data", "repr.12.data", force = true)

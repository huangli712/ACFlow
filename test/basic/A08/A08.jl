#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Printf
using ACFlow

welcome()

# For MaxEnt solver

# For diagonal part: giw.11.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.11.data",
    "mtype"  => "gauss",
    "ngrid"  => 20,
    "nmesh"  => 500,
    "wmax"   => 5.0,
    "wmin"   => -5.0,
    "beta"   => 40.0,
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
    "stype"  => "sj",
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

# For diagonal part: giw.22.data

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

# For auxiliary functions: giw.aux12.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.aux12.data",
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aaux12, Gaux12 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.aux12.data", force = true)
cp("Gout.data", "Gout.aux12.data", force = true)
cp("repr.data", "repr.aux12.data", force = true)

# For auxiliary functions: giw.aux21.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.aux21.data",
)
#
# For [MaxEnt] block
# See types.jl/_PMaxEnt for default setup
S = Dict{String,Any}(
)
#
setup_param(B, S, false)

# Call the solver
mesh, Aaux21, Gaux21 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.aux21.data", force = true)
cp("Gout.data", "Gout.aux21.data", force = true)
cp("repr.data", "repr.aux21.data", force = true)

# Generate model function at first
model_offdiag = sqrt.(Aout11 .* Aout22)
#
open("model.inp", "w") do fout
    for i in eachindex(mesh)
        @printf(fout, "%20.16f %20.16f\n", mesh[i], model_offdiag[i])
    end
end

# For non-diagonal part: giw.12.data

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

# For non-diagonal part: giw.21.data

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "finput" => "giw.21.data",
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
mesh, Aout21, Gout21 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.21.data", force = true)
cp("Gout.data", "Gout.21.data", force = true)
cp("repr.data", "repr.21.data", force = true)

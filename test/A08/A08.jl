#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using Printf
using ACFlow

welcome()

# For MaxEnt solver

# For diagonal elements: green.11.data

# Setup parameters
C = Dict{String,Any}(
    "finput" => "green.11.data",
    "mtype"  => "gauss",
    "ngrid"  => 20,
    "nmesh"  => 400,
    "wmax"   => 4.0,
    "wmin"   => -4.0,
    "beta"   => 40.0,
)
#
S = Dict{String,Any}(
    "nalph"  => 28,
    "alpha"  => 1e18,
)
#
setup_param(C, S)

# Call the solver
mesh, Aout11, Gout11 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.11.data", force = true)
cp("Gout.data", "Gout.11.data", force = true)
cp("repr.data", "repr.11.data", force = true)

# For diagonal elements: green.22.data

# Setup parameters
C = Dict{String,Any}(
    "finput" => "green.22.data",
)
#
S = Dict{String,Any}(
)
#
setup_param(C, S, false)

# Call the solver
mesh, Aout22, Gout22 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.22.data", force = true)
cp("Gout.data", "Gout.22.data", force = true)
cp("repr.data", "repr.22.data", force = true)

# For non-diagonal elements: green.12.data

# Generate model function at first
wmin = -4.0 # Left boundary
wmax = +4.0 # Right boundary
nmesh = 400 # Number of real-frequency points
#
model_offdiag = sqrt.(Aout11 .* Aout22)
#
w_real = collect(LinRange(wmin, wmax, nmesh))
#
open("model.data", "w") do fout
    for i in eachindex(w_real)
        @printf(fout, "%20.16f %20.16f\n", w_real[i], model_offdiag[i])
    end
end

# Setup parameters
C = Dict{String,Any}(
    "finput" => "green.12.data",
    "mtype"  => "file",
    "offdiag"=> true,
)
#
S = Dict{String,Any}(
    "nalph"  => 30,
    "alpha"  => 1e15,
)
#
setup_param(C, S, false)

# Call the solver
mesh, Aout12, Gout12 = solve(read_data())

# Backup calculated results
cp("Aout.data", "Aout.12.data", force = true)
cp("Gout.data", "Gout.12.data", force = true)
cp("repr.data", "repr.12.data", force = true)

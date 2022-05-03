#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using ACFlow
using Printf

welcome()

# Deal with self-energy function
#
# Read self-energy function
dlm = readdlm("sigma.data")
#
# Get grid
grid = dlm[:,1]
#
# Get self-energy function
Σinp = dlm[:,2] + im * dlm[:,3]
Σerr = dlm[:,4] + im * dlm[:,5]
#
# Subtract hartree term
Σ∞  = 1.0
@. Σinp = Σinp - Σ∞

# For MaxEnt solver

# Setup parameters
C = Dict{String,Any}(
    "mtype"  => "gauss",
    "mesh"   => "tangent",
    "ngrid"  => 100,
    "wmax"   => 10.0,
    "wmin"   => -10.0,
    "beta"   => 10.0,
)
#
S = Dict{String,Any}(
    "nalph"  => 15,
    "alpha"  => 1e12,
    "blur"   => -1.0,
)
#
setup_param(C, S)

# Call the solver
mesh, Aout, Gout = solve(grid, Σinp, Σerr)

# Calculate final self-energy function on real axis
#
# Construct final self-energy function
Σout = @. Gout + Σ∞
#
# Write self-energy function
open("sigma.real", "w") do fout
    for i in eachindex(mesh)
        z = Σout[i]
        @printf(fout, "%16.12f %16.12f %16.12f\n", mesh[i], real(z), imag(z))
    end
end

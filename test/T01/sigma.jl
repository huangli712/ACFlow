#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using Printf
using ACFlow

welcome()

# Deal with self-energy function
#
# Read self-energy function
dlm = readdlm("siw.inp")
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
B = Dict{String,Any}(
    "mtype"  => "gauss",
    "mesh"   => "tangent",
    "ngrid"  => 20,
    "nmesh"  => 801,
    "wmax"   => 8.0,
    "wmin"   => -8.0,
    "beta"   => 10.0,
)
#
S = Dict{String,Any}(
    "nalph"  => 15,
    "alpha"  => 1e12,
    "blur"   => -1.0,
)
#
setup_param(B, S)

# Call the solver
mesh, Aout, Σout = solve(grid, Σinp, Σerr)

# Calculate final self-energy function on real axis
#
# Construct final self-energy function
@. Σout = Σout + Σ∞
#
# Write self-energy function
open("sigma.data", "w") do fout
    for i in eachindex(mesh)
        z = Σout[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", mesh[i], real(z), imag(z))
    end
end

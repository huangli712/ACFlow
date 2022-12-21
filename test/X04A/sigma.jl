#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using Printf
using ACFlow

welcome()

# Deal with self-energy function
#
# Read self-energy function
dlm = readdlm("siw.data")
#
# Get grid
grid = dlm[:,1]
#
# Get self-energy function
Σinp = dlm[:,2] + im * dlm[:,3] # Value part
Σerr = dlm[:,4] + im * dlm[:,5] # Error part
#
# Subtract hartree term
Σ∞  = 1.0
@. Σinp = Σinp - Σ∞

# For StochPX solver

# Setup parameters
#
# For [BASE] block
# See types.jl/_PBASE for default setup
B = Dict{String,Any}(
    "solver" => "StochPX", # Choose StochPX solver
    "mtype"  => "gauss",   # Default model function
    "mesh"   => "tangent", # Mesh for spectral function
    "ngrid"  => 20,        # Number of grid points for input data
    "nmesh"  => 801,       # Number of mesh points for output data
    "wmax"   => 8.0,       # Right boundary of mesh
    "wmin"   => -8.0,      # Left boundary of mesh
    "beta"   => 10.0,      # Inverse temperature
)
#
# For [StochPX] block
# See types.jl/_PStochPX for default setup
S = Dict{String,Any}(
    "method" => "mean",    # How to evaluate the final spectral density
    "nfine"  => 100000,    # Number of points of a very fine linear mesh
    "npole"  => 2000,      # Number of poles
    "ntry"   => 2000,      # Number of attempts (tries) to seek the solution
    "nstep"  => 100,       # Number of Monte Carlo steps per attempt / try
    "theta"  => 1e+8,      # Artificial inverse temperature
    "eta"    => 1e-3,      # Tiny distance from the real axis
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
# Write self-energy function to sigma.data
open("sigma.data", "w") do fout
    for i in eachindex(mesh)
        z = Σout[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", mesh[i], real(z), imag(z))
    end
end

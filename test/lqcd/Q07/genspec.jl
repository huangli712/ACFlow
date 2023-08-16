#!/usr/bin/env julia

#
# Well, this script will generate a pseudo-spectrum, which is then used
# to create a nonlinear mesh for the StochPX solver.
#
#

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +9.0  # Left boundary
wmax = +16.  # Right boundary
nmesh = 2001 # Number of real-frequency points
𝑀₁   = 9.60  # Parameters for Gaussian mixture model
𝑀₂   = 11.5
Γ₁   = 0.01
Γ₂   = 0.20
𝐴₁   = 1.00
𝐴₂   = 0.40

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    image[i] =            𝐴₁ * exp(-(rmesh[i] - 𝑀₁) ^ 2.0 / Γ₁)
    image[i] = image[i] + 𝐴₂ * exp(-(rmesh[i] - 𝑀₂) ^ 2.0 / Γ₂)
end

# Write spectral function
open("Aout.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

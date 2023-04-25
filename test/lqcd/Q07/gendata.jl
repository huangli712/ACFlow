#!/usr/bin/env julia

#
# This example is taken from Phys. Rev. D 87, 114019 (2013).
# HTL Wilson line spectral function
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +10.  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
ntau = 501   # Number of imaginary time points
beta = 10.0  # Inverse temperature
𝑀₁   = 1.00  # Parameters for Gaussian mixture model
𝑀₂   = 3.00
𝑀₃   = 7.00
Γ₁   = 0.01
Γ₂   = 0.20
Γ₃   = 2.00
𝐴₁   = 1.00
𝐴₂   = 0.20
𝐴₃   = 0.10

#
# For true spectrum
#

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    image[i] =            𝐴₁ * exp(-(rmesh[i] - 𝑀₁) ^ 2.0 / Γ₁)
    image[i] = image[i] + 𝐴₂ * exp(-(rmesh[i] - 𝑀₂) ^ 2.0 / Γ₂)
    image[i] = image[i] + 𝐴₃ * exp(-(rmesh[i] - 𝑀₃) ^ 2.0 / Γ₃)
end
#
rmesh[1] = 1e-8 # To avoid NaN
image = image ./ rmesh
rmesh[1] = 0.0

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

#
# For Matsubara frequency data
#

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 0.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise = randn(rng, F64, niw) * noise_ampl

# Kernel function
kernel = -2.0 * reshape(rmesh .^ 2.0, (1,nmesh)) ./
         (reshape(iw .^ 2.0, (niw,1)) .+ reshape(rmesh .^ 2.0, (1,nmesh)))
kernel[1,1] = -2.0

# Build green's function
KA = kernel .* reshape(image, (1,nmesh))
chiw = zeros(F64, niw)
for i in eachindex(chiw)
    chiw[i] = trapz(rmesh, KA[i,:]) + noise[i]
end

# Build error
err = ones(F64, niw) * noise_ampl

# Write green's function
open("chiw.data", "w") do fout
    for i in eachindex(chiw)
        z = chiw[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", iw[i], z, err[i])
    end
end

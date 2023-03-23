#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +10.0 # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
ntau = 501   # Number of imaginary time points
beta = 50.0  # Inverse temperature
𝑀    = 2.00  # Parameters for Breit-Wigner model
Γ    = 0.50
𝐴    = 1.00

#
# For true spectrum
#

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    B₁ = (𝑀 ^ 2.0 + Γ ^ 2.0 - rmesh[i] ^ 2.0) ^ 2.0
    B₂ = 4.0 * (Γ ^ 2.0) * (rmesh[i] ^ 2.0)
    image[i] = 4.0 * 𝐴 * Γ * rmesh[i] / (B₁ + B₂)
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

#
# For imaginary time data
#

# Imaginary time mesh
tmesh = collect(LinRange(0, beta, ntau))

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise = randn(rng, F64, ntau) * noise_ampl

# Build green's function
chit = zeros(F64, ntau)
for i = 1:ntau
    tw = exp.(-tmesh[i] * rmesh)
    bw = exp.(-beta * rmesh)
    btw = exp.(-(beta - tmesh[i]) * rmesh)
    K = rmesh .* (tw .+ btw) ./ (1.0 .- bw)
    K[1] = 2.0 / beta
    global KA = K .* image
    chit[i] = trapz(rmesh, KA) + noise[i]
end

# Build error
err = ones(F64, ntau) * noise_ampl

# Write green's function
open("chit.data", "w") do fout
    for i in eachindex(chit)
        @printf(fout, "%16.12f %16.12f %16.12f\n", tmesh[i], chit[i], err[i])
    end
end

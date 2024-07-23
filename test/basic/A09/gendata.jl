#!/usr/bin/env julia

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
ntau = 1001  # Number of imaginary time points
beta = 5.00  # Inverse temperature
ϵ₁   = 2.00  # Parameters for gaussian peaks
ϵ₂   = -2.0
Γ₁   = 0.50
Γ₂   = 0.50

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image  = exp(-(rmesh - ϵ₁) ^ 2.0 / (2.0 * Γ₁ ^ 2.0)) / (sqrt(2.0 * π) * Γ₁)
@. image += exp(-(rmesh - ϵ₂) ^ 2.0 / (2.0 * Γ₂ ^ 2.0)) / (sqrt(2.0 * π) * Γ₂)
#
image = image ./ trapz(rmesh, image)

# Imaginary time mesh
tmesh = collect(LinRange(0, beta, ntau))

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise = randn(rng, F64, ntau) * noise_ampl

# Build green's function
gtau = zeros(F64, ntau)
for i = 1:ntau
    tw = exp.(-tmesh[i] * rmesh)
    bw = exp.(-beta * rmesh)
    gtau[i] = trapz(rmesh, image .* tw ./ (1.0 .+ bw)) + noise[i]
end

# Build error
err = ones(F64, ntau) * noise_ampl * 10.0

# Write green's function
open("gtau.data", "w") do fout
    for i in eachindex(gtau)
        @printf(fout, "%20.16f %20.16f %20.16f\n", tmesh[i], gtau[i], err[i])
    end
end

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

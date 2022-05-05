#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +8.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
ntau = 1000  # Number of imaginary time points
beta = 5.00  # Inverse temperature
ϵ₁   = 1.00  # Parameters for gaussian peaks
ϵ₂   = -1.0
ϵ₃   = 6.00
ϵ₄   = -6.0
A₁   = 1.00
A₂   = 1.00
A₃   = 0.20
A₄   = 0.20
Γ₁   = 1.00
Γ₂   = 1.00
Γ₃   = 0.50
Γ₄   = 0.50

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image  = A₁ * exp(-(rmesh - ϵ₁) ^ 2.0 / (2.0 * Γ₁ ^ 2.0)) / (sqrt(2.0 * π) * Γ₁)
@. image -= A₂ * exp(-(rmesh - ϵ₂) ^ 2.0 / (2.0 * Γ₂ ^ 2.0)) / (sqrt(2.0 * π) * Γ₂)
@. image += A₃ * exp(-(rmesh - ϵ₃) ^ 2.0 / (2.0 * Γ₃ ^ 2.0)) / (sqrt(2.0 * π) * Γ₃)
@. image -= A₄ * exp(-(rmesh - ϵ₄) ^ 2.0 / (2.0 * Γ₄ ^ 2.0)) / (sqrt(2.0 * π) * Γ₄)
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
    btw = exp.(-(beta - tmesh[i]) * rmesh)
    K = 0.5 * rmesh .* (tw .+ btw) ./ (1.0 .- bw)
    K[1] = 1.0 / beta
    KA = K .* image
    gtau[i] = trapz(rmesh, KA)
end

# Build error
err = ones(F64, ntau) * noise_ampl

# Write green's function
open("chit.data", "w") do fout
    for i in eachindex(gtau)
        @printf(fout, "%20.16f %20.16f %20.16f\n", tmesh[i], gtau[i], err[i])
    end
end

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f %20.16f\n",
            rmesh[i], image[i], rmesh[i] * image[i])
    end
end

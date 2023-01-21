#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 20.0  # Inverse temperature
ϵ₁   = 2.00  # Parameters for δ-like peaks
ϵ₂   = -2.0
ϵ₃   = 1.00
ϵ₄   = -1.0
A₁   = 0.25
A₂   =-0.25
A₃   = 0.50
A₄   =-0.15
η    = 1e-2

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 0.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise_abs = randn(rng, F64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Build green's function
giw = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = (
        A₁ / (iωₙ[i] * im - ϵ₁) +
        A₂ / (iωₙ[i] * im - ϵ₂) +
        A₃ / (iωₙ[i] * im - ϵ₃) +
        A₄ / (iωₙ[i] * im - ϵ₄) + noise[i]
    )
end
#
gre = zeros(C64, nmesh)
for i in eachindex(gre)
    gre[i] = (
        A₁ / (ω[i] + η * im - ϵ₁) +
        A₂ / (ω[i] + η * im - ϵ₂) +
        A₃ / (ω[i] + η * im - ϵ₃) +
        A₄ / (ω[i] + η * im - ϵ₄)
    )
end

# Build error
err = ones(F64, niw) * noise_ampl

# Write green's function (Matsubara frequency axis)
open("giw.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iωₙ[i], real(z), imag(z), err[i])
    end
end

# Write green's function (real frequency axis)
open("gre.data", "w") do fout
    for i in eachindex(gre)
        z = gre[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", ω[i], real(z), imag(z))
    end
end

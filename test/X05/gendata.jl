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
beta = 10.0  # Inverse temperature
ϵ₁   = 0.70  # Parameters for δ-like peaks
ϵ₂   = 1.20
A₁   = 1.10
A₂   = 1.335663
B₁   = A₁ / (2.0 * ϵ₁)
B₂   = A₂ / (2.0 * ϵ₂)
η    = 0.05

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
giw2 = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = (
        A₁ / ((iωₙ[i] * im) ^ 2.0 - ϵ₁ ^ 2.0) +
        A₂ / ((iωₙ[i] * im) ^ 2.0 - ϵ₂ ^ 2.0)
    )
    giw2[i] = (
        B₁ / (iωₙ[i] * im - ϵ₁) + (-B₁) / (iωₙ[i] * im + ϵ₁) +
        B₂ / (iωₙ[i] * im - ϵ₂) + (-B₂) / (iωₙ[i] * im + ϵ₂)
    )
end
#
gre = zeros(C64, nmesh)
gre2 = zeros(C64, nmesh)
for i in eachindex(gre)
    gre[i] = (
        A₁ / ((ω[i] + η * im) ^ 2.0 - ϵ₁ ^ 2.0) +
        A₂ / ((ω[i] + η * im) ^ 2.0 - ϵ₂ ^ 2.0)
    )
    gre2[i] = (
        B₁ / (ω[i] + η * im - ϵ₁) + (-B₁) / (ω[i] + η * im + ϵ₁) +
        B₂ / (ω[i] + η * im - ϵ₂) + (-B₂) / (ω[i] + η * im + ϵ₂)
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

open("giw2.data", "w") do fout
    for i in eachindex(giw)
        z = giw2[i]
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

open("gre2.data", "w") do fout
    for i in eachindex(gre)
        z = gre2[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", ω[i], real(z), imag(z))
    end
end

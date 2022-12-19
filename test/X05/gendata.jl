#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -10.0 # Left boundary
wmax = +10.0 # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 10.0  # Inverse temperature
ϵ₀   = 9.00  # Parameters for δ-like peaks
ϵ₁   = 3.00  
ϵ₂   = 0.50
ϵ₃   = 0.00
ϵ₄   =-0.50
ϵ₅   =-3.00
ϵ₆   =-6.00
ϵ₇   =-9.00
A₀   = 0.125
A₁   = 0.125
A₂   = 0.125
A₃   = 0.125
A₄   = 0.125
A₅   = 0.125
A₆   = 0.125
A₇   = 0.125
η    = 1e-2

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

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
        A₀ / (iωₙ[i] * im - ϵ₀) +
        A₁ / (iωₙ[i] * im - ϵ₁) +
        A₂ / (iωₙ[i] * im - ϵ₂) +
        A₃ / (iωₙ[i] * im - ϵ₃) +
        A₄ / (iωₙ[i] * im - ϵ₄) +
        A₅ / (iωₙ[i] * im - ϵ₅) +
        A₆ / (iωₙ[i] * im - ϵ₆) +
        A₇ / (iωₙ[i] * im - ϵ₇) + noise[i]
    )
end
#
gre = zeros(C64, nmesh)
for i in eachindex(gre)
    gre[i] = (
        A₀ / (ω[i] + η * im - ϵ₀) +
        A₁ / (ω[i] + η * im - ϵ₁) +
        A₂ / (ω[i] + η * im - ϵ₂) +
        A₃ / (ω[i] + η * im - ϵ₃) +
        A₄ / (ω[i] + η * im - ϵ₄) +
        A₅ / (ω[i] + η * im - ϵ₅) +
        A₆ / (ω[i] + η * im - ϵ₆) +
        A₇ / (ω[i] + η * im - ϵ₇)
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

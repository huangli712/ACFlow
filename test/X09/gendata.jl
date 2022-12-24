#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -4.0 # Left boundary
wmax = +4.0 # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 20.0  # Inverse temperature
ϵ₀   = 36.0 / 10  # Parameters for δ-like peaks
ϵ₁   = 12.0 / 10
ϵ₂   = 2.00 / 10
ϵ₃   = 0.00 / 10
ϵ₄   =-2.00 / 10
ϵ₅   =-12.0 / 10
ϵ₆   =-24.0 / 10
ϵ₇   =-36.0 / 10
ϵ₈   =-8.00 / 10
ϵ₉   = 10.0 / 10
A₀   = 0.10
A₁   = 0.10
A₂   = 0.10
A₃   = 0.10
A₄   = 0.10
A₅   = 0.10
A₆   = 0.10
A₇   = 0.10
A₈   = 0.10
A₉   = 0.10
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
        A₇ / (iωₙ[i] * im - ϵ₇) +
        A₈ / (iωₙ[i] * im - ϵ₈) +
        A₉ / (iωₙ[i] * im - ϵ₉) + noise[i]
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
        A₇ / (ω[i] + η * im - ϵ₇) +
        A₈ / (ω[i] + η * im - ϵ₈) +
        A₉ / (ω[i] + η * im - ϵ₉)
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

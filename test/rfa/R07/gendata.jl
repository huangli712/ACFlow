#!/usr/bin/env julia

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -1.5  # Left boundary
wmax = +1.5  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 100   # Number of Matsubara frequencies
beta = 50.0  # Inverse temperature
η    = 1e-2  # Parameters for δ-like peaks
ϵ    = [-1.0, -0.3, -0.1, 0.1, 0.3, 1.0]
A    = [-1.0, -2.0, -3.0, 3.0, 2.0, 1.0]

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 0.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-8
noise_abs = randn(rng, F64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Build green's function
giw = zeros(C64, niw)
for i in eachindex(giw)
    for j in eachindex(ϵ)
        giw[i] = giw[i] + A[j] / (iωₙ[i] * im - ϵ[j])
    end
    giw[i] = giw[i] + noise[i]
end
#
gre = zeros(C64, nmesh)
for i in eachindex(gre)
    for j in eachindex(ϵ)
        gre[i] = gre[i] + A[j] / (ω[i] + η * im - ϵ[j])
    end
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

#!/usr/bin/env julia

#
# This example is taken from Phys. Rev. D 106, L051502 (2022)
# Resonance peak + Continuum spectral function
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), one(x))

function ξ(ω, 𝑀, Δ)
    return 1.0 / ( 1.0 + exp( (𝑀 ^ 2.0 - ω ^ 2.0) / (ω * Δ) ) )
end

function ρᵣ(ω, 𝐶, 𝑀, Γ)
    return 𝐶 * ω ^ 2.0 / ( (ω ^ 2.0 / 𝑀 / Γ - 𝑀 / Γ) ^ 2.0 + 1.0 )
end

function ρₜ(ω, 𝐶, 𝑀, β)
    𝑃₁ = 𝐶 * 3.0 * ω ^ 2.0 / 8.0 / π
    𝑃₂ = heaviside(ω ^ 2.0 - 4.0 * 𝑀 ^ 2.0)
    𝑃₃ = tanh(ω * β / 4.0)
    𝑃₄ = sqrt(abs(1.0 - (2.0 * 𝑀 / ω) ^ 2.0))
    𝑃₅ = 2.0 + (2.0 * 𝑀 / ω) ^ 2.0
    return 𝑃₁ * 𝑃₂ * 𝑃₃ * 𝑃₄ * 𝑃₅
end

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +4.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 20    # Number of Matsubara frequencies
ntau = 501   # Number of imaginary time points
beta = 50.0  # Inverse temperature
𝑀ᵣ   = 0.10  # Parameters for model
𝑀ₜ   = 0.05
𝐶ᵣ   = 2.00
𝐶ₜ   = 2.10
Γ    = 0.03
Δ    = 1.00

#
# For true spectrum
#

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
rmesh[1] = 1e-8 # To avoid NaN
for i in eachindex(rmesh)
    ρ₁ = ξ(rmesh[i], 𝑀ᵣ, Γ) * ρᵣ(rmesh[i], 𝐶ᵣ, 𝑀ᵣ, Γ) * (1.0 - ξ(rmesh[i], 𝑀ᵣ + Γ, Γ))
    ρ₂ = ξ(rmesh[i], 𝑀ᵣ + Γ, Γ) * ρₜ(rmesh[i], 𝐶ₜ, 𝑀ₜ, beta)
    image[i] = ρ₁ + ρ₂
end
#
image = image ./ rmesh ./ rmesh ./ rmesh
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

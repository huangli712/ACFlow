#!/usr/bin/env julia

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -10.0 # Left boundary
wmax = +10.0 # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 200   # Number of Matsubara frequencies
beta = 100.0 # Inverse temperature
ϵ₁   = 0.00  # Parameters for gaussian peaks
ϵ₂   = -2.5
ϵ₃   = 2.50
A₁   = 0.50
A₂   = 0.30
A₃   = 0.20
Γ₁   = 0.50
Γ₂   = 0.80
Γ₃   = 0.80

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image  = A₁ / π * Γ₁ / ((rmesh - ϵ₁) ^ 2.0 + Γ₁ ^ 2.0)
@. image += A₂ / π * Γ₂ / ((rmesh - ϵ₂) ^ 2.0 + Γ₂ ^ 2.0)
@. image += A₃ * exp(-0.5 * ((rmesh - ϵ₃) / Γ₃) ^ 2.0) / (sqrt(2.0 * π) *  Γ₃)
#
image = image ./ trapz(rmesh, image)

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-8
noise_abs = randn(rng, F64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(rmesh, (1,nmesh)))

# Build green's function
KA = kernel .* reshape(image, (1,nmesh))
giw = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = trapz(rmesh, KA[i,:]) + noise[i]
end

# Build error
err = ones(F64, niw) * noise_ampl

# Write green's function
open("giw.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err[i])
    end
end

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

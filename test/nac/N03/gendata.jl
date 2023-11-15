#!/usr/bin/env julia

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 1000 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
beta = 100.0  # Inverse temperature
ϵ₁   = 2.50  # Parameters for gaussian peaks
ϵ₂   = -2.5
A₁   = 0.50
A₂   = 0.50
Γ₁   = 0.50
Γ₂   = 0.50

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image  = A₁ * exp(-(rmesh - ϵ₁) ^ 2.0 / (2.0 * Γ₁ ^ 2.0))
@. image += A₂ * exp(-(rmesh - ϵ₂) ^ 2.0 / (2.0 * Γ₂ ^ 2.0))
#
image = image ./ trapz(rmesh, image)

gaussian(x, mu, sigma) = exp(-0.5*((x-mu)/sigma)^2)/(sqrt(2*π)*sigma)
rho(omega) = 0.45*gaussian(omega, -2.5, 0.7) + 0.1*gaussian(omega, 0.0, 0.1) + 0.45*gaussian(omega, 2.5, 0.7)
omegas = LinRange(-10, 10, 1000)

rmesh = omegas
image = rho.(omegas)

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 0.0e-4
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
        @printf(fout, "%30.24f %30.24f %30.24f %30.24f\n", iw[i], real(z), imag(z), err[i])
    end
end

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

open("image1.data", "w") do fout
    image = rho.(omegas)
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", omegas[i], image[i])
    end
end

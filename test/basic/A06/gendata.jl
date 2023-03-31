#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -10.0 # Left boundary
wmax = +10.0 # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
beta = 20.0  # Inverse temperature
ϵ    = 0.50  # Parameters for gaussian peaks
Γ    = 1.00

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image = exp(-(rmesh - ϵ) ^ 2.0 / (2.0 * Γ ^ 2.0))
#
image = image ./ trapz(rmesh, image)

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise = randn(rng, F64, niw) + im * randn(rng, F64, niw)
noise = noise_ampl * noise / sqrt(2.0)

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

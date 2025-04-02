#!/usr/bin/env julia

haskey(ENV,"ACFLOW_PATH") && pushfirst!(LOAD_PATH, ENV["ACFLOW_PATH"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +8.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 100   # Number of Matsubara frequencies
beta = 50.0  # Inverse temperature
W₁   = 0.30  # Parameters for gaussian peaks
W₂   = 0.20
Γ₁   = 0.30
Γ₂   = 1.20
Γ₃   = 4.00
ϵ    = 3.00

#
# For true spectrum
#

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    A = W₁ / (1.0 + ((rmesh[i] - 0) / Γ₁) ^ 2.0) +
        W₂ / (1.0 + ((rmesh[i] - ϵ) / Γ₂) ^ 2.0) +
        W₂ / (1.0 + ((rmesh[i] + ϵ) / Γ₂) ^ 2.0)
    image[i] = A / (1.0 + (rmesh[i] / Γ₃) ^ 6.0)
end
#
image = image ./ trapz(rmesh, image)

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

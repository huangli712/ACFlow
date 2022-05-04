#!/usr/bin/env julia

using Random
using Printf

# Numerical integration by the composite trapezoidal rule
function trapz(x, y, linear::Bool = false)
    # For linear mesh
    if linear
        h = x[2] - x[1]
        value = y[1] + y[end] + 2.0 * sum(y[2:end-1])
        value = h * value / 2.0
    # For non-equidistant mesh
    else
        len = length(x)
        value = 0.0
        for i = 1:len-1
            value = value + (y[i] + y[i+1]) * (x[i+1] - x[i])
        end
        value = value / 2.0
    end

    return value
end

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 10.0  # Inverse temperature
ϵ₁   = 0.50  # Parameters for gaussian peaks
ϵ₂   = -2.5
A₁   = 1.00
A₂   = 0.30
Γ₁   = 0.20
Γ₂   = 0.80

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
@. image =  A₁ * exp(-(rmesh - ϵ₁) ^ 2.0 / (2.0 * Γ₁ ^ 2.0))
@. image += A₂ * exp(-(rmesh - ϵ₂) ^ 2.0 / (2.0 * Γ₂ ^ 2.0))
#
image = image ./ trapz(rmesh, image)

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise_abs = randn(rng, Float64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(rmesh, (1,nmesh)))

# Build green's function
KA = kernel .* reshape(image, (1,nmesh))
giw = zeros(ComplexF64, niw)
for i in eachindex(giw)
    giw[i] = trapz(rmesh, KA[i,:]) + noise[i]
end

# Build error
err = ones(Float64, niw) * noise_ampl

# Write green's function
open("green.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err[i])
    end
end

# Write spectral function
open("exact.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

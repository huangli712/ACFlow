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
wmin = +0.0  # Left boundary
wmax = +8.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
ntau = 1000  # Number of imaginary time points
beta = 20.0  # Inverse temperature
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
w_real = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
spec_real = similar(w_real)
#
for i in eachindex(w_real)
    A = W₁ / (1.0 + (w_real[i] / Γ₁) ^ 2.0) +
        W₂ / (1.0 + ((w_real[i] - ϵ) / Γ₂) ^ 2.0) +
        W₂ / (1.0 + ((w_real[i] + ϵ) / Γ₂) ^ 2.0)
    spec_real[i] = A / (1.0 + (w_real[i] / Γ₃) ^ 6.0)
end
#
spec_real = spec_real ./ trapz(w_real, spec_real)

# Write spectral function
open("exact.data", "w") do fout
    for i in eachindex(spec_real)
        @printf(fout, "%20.16f %20.16f\n", w_real[i], spec_real[i])
    end
end

#
# For Matsubara frequency data
#

# Matsubara frequency mesh
iw = 2.0 * π / beta * collect(0:niw-1)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_amplitude = 1.0e-4
noise = randn(rng, Float64, niw) * noise_amplitude

# Kernel function
kernel = reshape(w_real .^ 2.0, (1,nmesh)) ./
         (reshape(iw .^ 2.0, (niw,1)) .+ reshape(w_real .^ 2.0, (1,nmesh)))
kernel[1,1] = 1.0

# Build green's function
KA = kernel .* reshape(spec_real, (1,nmesh))
gf_mats = zeros(Float64, niw)
for i in eachindex(gf_mats)
    gf_mats[i] = trapz(w_real, KA[i,:]) + noise[i]
end
norm = gf_mats[1]
gf_mats = gf_mats / norm

# Build error
err = ones(Float64, niw) * noise_amplitude / norm

# Write green's function
open("chiw.data", "w") do fout
    for i in eachindex(gf_mats)
        z = gf_mats[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", iw[i], z, err[i])
    end
end

#
# For imaginary time data
#

# Imaginary time mesh
t_mesh = collect(LinRange(0, beta, ntau))

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_amplitude = 1.0e-4
noise = randn(rng, Float64, ntau) * noise_amplitude

# Build green's function
gtau = zeros(Float64, ntau)
for i = 1:ntau
    tw = exp.(-t_mesh[i] * w_real)
    bw = exp.(-beta * w_real)
    btw = exp.(-(beta - t_mesh[i]) * w_real)
    K = 0.5 * w_real .* (tw .+ btw) ./ (1.0 .- bw)
    K[1] = 1.0 / beta
    global KA = K .* spec_real
    gtau[i] = trapz(w_real, KA)
end

# Build error
err = ones(Float64, ntau) * noise_amplitude

# Write green's function
open("chit.data", "w") do fout
    for i in eachindex(gtau)
        @printf(fout, "%16.12f %16.12f %16.12f\n", t_mesh[i], gtau[i], err[i])
    end
end

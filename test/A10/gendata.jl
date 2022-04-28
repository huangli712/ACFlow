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
nmesh = 801  # Number of real-frequency points
ntau = 501   # Number of imaginary time points
beta = 5.0   # Inverse temperature

# Real frequency mesh
w_real = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
spec_real = similar(w_real)
@. spec_real  = exp(-0.5 * (w_real - 1.0) ^ 2.0 / (1.0 ^ 2.0)) / (sqrt(2.0 * π) * 0.5)
@. spec_real -= exp(-0.5 * (w_real + 1.0) ^ 2.0 / (1.0 ^ 2.0)) / (sqrt(2.0 * π) * 0.5)
@. spec_real += exp(-0.5 * (w_real - 6.0) ^ 2.0 / (0.5 ^ 2.0)) / (sqrt(2.0 * π) * 0.5) * 0.2
@. spec_real -= exp(-0.5 * (w_real + 6.0) ^ 2.0 / (0.5 ^ 2.0)) / (sqrt(2.0 * π) * 0.5) * 0.2
spec_real = spec_real ./ trapz(w_real, spec_real)

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
    KA = K .* spec_real
    gtau[i] = trapz(w_real, KA)
end

# Build error
err = ones(Float64, ntau) * noise_amplitude

# Write green's function
open("green.data", "w") do fout
    for i in eachindex(gtau)
        @printf(fout, "%16.12f %16.12f %16.12f\n", t_mesh[i], gtau[i], err[i])
    end
end

# Write spectral function
open("exact.data", "w") do fout
    for i in eachindex(spec_real)
        @printf(fout, "%16.12f %16.12f %16.12f\n",
            w_real[i], spec_real[i], w_real[i] * spec_real[i])
    end
end

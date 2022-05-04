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
ntau = 1000  # Number of imaginary time points
beta = 10.0  # Inverse temperature
Δ    = 0.5   # 2Δ is the size of the gap
W    = 6.0   # Bandwidth of the spectrum

#
# For true spectrum
#

# Real frequency mesh
w_real = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
spec_real = similar(w_real)
for i in eachindex(w_real)
    spec_real[i] = 0.0
    if Δ < abs(w_real[i]) < W/2
        spec_real[i] = abs(w_real[i]) / sqrt(w_real[i] ^ 2.0 - Δ ^ 2.0) / W
    end
end
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
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_amplitude = 1.0e-4
noise_abs = randn(rng, Float64, niw) * noise_amplitude
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(w_real, (1,nmesh)))

# Build green's function
KA = kernel .* reshape(spec_real, (1,nmesh))
gf_mats = zeros(ComplexF64, niw)
for i in eachindex(gf_mats)
    gf_mats[i] = trapz(w_real, KA[i,:]) + noise[i]
end

# Build error
err = ones(Float64, niw) * noise_amplitude

# Write green's function
open("giw.data", "w") do fout
    for i in eachindex(gf_mats)
        z = gf_mats[i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err[i])
    end
end

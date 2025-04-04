#!/usr/bin/env julia

haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
ntau = 1000  # Number of imaginary time points
beta = 10.0  # Inverse temperature
Δ    = 0.50  # 2Δ is the size of the gap
W    = 6.00  # Bandwidth of the spectrum

#
# For true spectrum
#

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    image[i] = 0.0
    if Δ < abs(rmesh[i]) < W/2
        image[i] = abs(rmesh[i]) / sqrt(rmesh[i] ^ 2.0 - Δ ^ 2.0) / W
    end
end
#
image = image ./ trapz(rmesh, image)
_rmesh = deepcopy(rmesh)
_rmesh[1] = _rmesh[1] + 1e-10
@. image  = image / _rmesh

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
giw = zeros(F64, niw)
for i in eachindex(giw)
    giw[i] = trapz(rmesh, KA[i,:]) + noise[i]
end

# Build error
err = ones(F64, niw) * noise_ampl

# Write green's function
open("giw.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", iw[i], z, err[i])
    end
end

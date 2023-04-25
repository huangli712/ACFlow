#!/usr/bin/env julia

#
# This example is taken from Phys. Rev. D 87, 114019 (2013).
# HTL Wilson line spectral function
#
# Note: Please change Araw1.inp to Araw2.inp in L28 to use another input.
# Araw1.inp -> r = 0.066fm (Fig. 5)
# Araw2.inp -> r = 0.466fm (Fig. 5)
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using DelimitedFiles
using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +8.0  # Right boundary
niw  = 50    # Number of Matsubara frequencies
beta = 1.60  # Inverse temperature

#
# For true spectrum
#

# Read raw spectral function
dlm = readdlm("Araw1.inp", F64)
nmesh, _ = size(dlm)

# Real frequency mesh
rmesh = dlm[:,1]

# Spectral function
image = dlm[:,2]
image = image ./ rmesh

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

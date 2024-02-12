#!/usr/bin/env julia

#
# This example is taken from Phys. Rev. D 97, 094503 (2018).
# 
# Stochastic reconstructions of spectral functions: Application to lattice QCD
#
# H.-T. Ding, et al.
#
# Araw.inp -> Fig. 16
#
# Tc = 269.64 MeV, T = 0.75Tc 
#

using DelimitedFiles
using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +10.  # Right boundary
niw  = 50    # Number of Matsubara frequencies
beta = 5.00  # Inverse temperature

#
# For true spectrum
#

# Read raw spectral function
dlm = readdlm("Araw.inp", F64)
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

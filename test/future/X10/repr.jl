#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 20    # Number of Matsubara frequencies
beta = 40.0  # Inverse temperature
ϵ₁   = 0.50  # Parameters for gaussian peaks
ϵ₂   = -0.5
A₁   = -0.10
A₂   = +0.10
η₁   = 1e-2
η₂   = 1e-2

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Initial green's function (in Matsubara axis)
giw = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = (
        A₁ / (iωₙ[i] * im - ϵ₁) + A₂ / (iωₙ[i] * im - ϵ₂)
    )
end

# Initial green's function (in real axis)
gre = zeros(C64, nmesh)
for i in eachindex(gre)
    gre[i] = (
        A₁ / (ω[i] + η₁ * im - ϵ₁) + A₂ / (ω[i] + η₂ * im - ϵ₂)
    )
end

# Write green's function
open("giw.data", "w") do fout
    for i = 1:niw
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", iωₙ[i], real(z), imag(z))
    end
end
#
open("gre.data", "w") do fout
    for i = 1:nmesh
        z = gre[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", ω[i], real(z), imag(z))
    end
end

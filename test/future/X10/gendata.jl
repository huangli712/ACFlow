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
ϵ₁   = 1.00  # Parameters for gaussian peaks
ϵ₂   = 2.00
A₁   = 1.00
A₂   = 1.00
η₁   = 1e-2
η₂   = 1e-2

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Initial green's function (in Matsubara axis)
giw1 = zeros(C64, niw)
for i in eachindex(giw1)
    giw1[i] = (
        A₁ / (iωₙ[i] * im - ϵ₁) + noise[i]
    )
end
#
giw2 = zeros(C64, niw)
for i in eachindex(giw2)
    giw2[i] = (
        A₂ / (iωₙ[i] * im - ϵ₂) + noise[i]
    )
end
#
𝔾iw = zeros(C64, (2,2,niw))
𝔾iw[1,1,:] .= giw1
𝔾iw[2,2,:] .= giw2

# Initial green's function (in real axis)
gre1 = zeros(C64, nmesh)
for i in eachindex(gre1)
    gre1[i] = (
        A₁ / (ω[i] + η * im - ϵ₁)
    )
end
#
gre2 = zeros(C64, nmesh)
for i in eachindex(gre2)
    gre2[i] = (
        A₂ / (ω[i] + η * im - ϵ₂)
    )
end
#
𝔾re = zeros(C64, (2,2,nmesh))
𝔾re[1,1,:] .= gre1
𝔾re[2,2,:] .= gre2

# Rotate green's function to generate non-diagonal element
#
# Set rotation angle
θ = 0.1
#
# Build rotation matrix
ℝ = [cos(θ) sin(θ); -sin(θ) cos(θ)]
#
# Get final green's function (in Matsubara axis)
𝒢iw = zeros(C64, (2,2,niw))
for w = 1:niw
    𝒢iw[:,:,w] = ℝ * 𝔾iw[:,:,w] * ℝ'
end
#
# Get final green's function (in real axis)
𝒢re = zeros(C64, (2,2,nmesh))
for w = 1:nmesh
    𝒢re[:,:,w] = ℝ * 𝔾re[:,:,w] * ℝ'
end

# Build error
err = 1e-5

# Write green's function
# For diagonal part
open("giw.11.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
# For non-diagonal part
open("giw.12.data", "w") do fout
    for i = 1:niw
        z = giw[1,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
# For non-diagonal part
open("giw.21.data", "w") do fout
    for i = 1:niw
        z = giw[2,1,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
# For diagonal part
open("giw.22.data", "w") do fout
    for i = 1:niw
        z = giw[2,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
# Be careful, we use two different auxiliary green's functions here
# For auxiliary green's function
open("giw.aux12.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i] + giw[2,2,i] + 2 * giw[1,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
# For auxiliary green's function
open("giw.aux21.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i] + giw[2,2,i] - 2 * giw[2,1,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end

# Write spectral function
# For diagonal part
open("image.11.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], 𝒜[1,1,i])
    end
end
#
# For non-diagonal part
open("image.12.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], 𝒜[1,2,i])
    end
end
#
# For non-diagonal part
open("image.21.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], 𝒜[2,1,i])
    end
end
#
# For diagonal part
open("image.22.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], 𝒜[2,2,i])
    end
end

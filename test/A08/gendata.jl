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
ϵ₃   = -1.0
ϵ₄   = -2.1
A₁   = 0.50
A₂   = 0.50
A₃   = 0.50
A₄   = 0.50
Γ₁   = 0.20
Γ₂   = 0.70
Γ₃   = 0.25
Γ₄   = 0.60

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Initial spectral function
image1 = similar(rmesh)
@. image1  = A₁ * exp(-(rmesh - ϵ₁) ^ 2.0 / (2.0 * Γ₁ ^ 2.0)) / (Γ₁ * sqrt(2.0 * π))
@. image1 += A₂ * exp(-(rmesh - ϵ₂) ^ 2.0 / (2.0 * Γ₂ ^ 2.0)) / (Γ₂ * sqrt(2.0 * π))
#
image2 = similar(rmesh)
@. image2  = A₃ * exp(-(rmesh - ϵ₃) ^ 2.0 / (2.0 * Γ₃ ^ 2.0)) / (Γ₃ * sqrt(2.0 * π))
@. image2 += A₄ * exp(-(rmesh - ϵ₄) ^ 2.0 / (2.0 * Γ₄ ^ 2.0)) / (Γ₄ * sqrt(2.0 * π))
#
𝐴 = zeros(F64, (2,2,nmesh))
𝐴[1,1,:] .= image1
𝐴[2,2,:] .= image2

# Rotate spectral function to generate non-diagonal element
#
# Set rotation angle
θ = 0.1
#
# Build rotation matrix
𝑅 = [cos(θ) sin(θ); -sin(θ) cos(θ)]
#
# Get final spectral function
true_image = zeros(F64, (2,2,nmesh))
for w = 1:nmesh
    true_image[:,:,w] = 𝑅 * 𝐴[:,:,w] * 𝑅'
end

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(rmesh, (1,nmesh)))

# Build green's function
KA = reshape(kernel, (1,1,niw,nmesh)) .* reshape(true_image, (2,2,1,nmesh))
giw = zeros(C64, (2,2,niw))
for i = 1:2
    for j = 1:2
        for w = 1:niw
            giw[i,j,w] = trapz(rmesh, KA[i,j,w,:])
        end
    end
end

# Build error
err = 1e-5

# Write green's function
open("giw.11.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
open("giw.12.data", "w") do fout
    for i = 1:niw
        z = giw[1,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
open("giw.22.data", "w") do fout
    for i = 1:niw
        z = giw[2,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end
#
open("gdiff.data", "w") do fout
    for i = 1:niw
        z = giw[1,1,i] + giw[2,2,i] - 2 * giw[1,2,i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iw[i], real(z), imag(z), err)
    end
end

# Write spectral function
open("image.11.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], true_image[1,1,i])
    end
end
#
open("image.12.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], true_image[1,2,i])
    end
end
#
open("image.22.data", "w") do fout
    for i in eachindex(rmesh)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], true_image[2,2,i])
    end
end

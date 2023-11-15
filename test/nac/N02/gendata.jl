#!/usr/bin/env julia

using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 1001 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
beta = 100.0 # Inverse temperature
ϵ₁   = -1.9  # Parameters for gaussian peaks
ϵ₂   = 1.90
A₁   = 0.50
A₂   = 0.50
Γ₁   = 0.50
Γ₂   = 0.50

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
gaussian(x, μ, σ) = exp(-0.5 * ( ( x - μ ) / σ ) ^ 2.0) / (sqrt(2*π) * σ)
rho(ω) = A₁ * gaussian(ω, ϵ₁, Γ₁) + A₂ * gaussian(ω, ϵ₂, Γ₂)
image = rho.(rmesh)

# Matsubara frequency mesh
iw = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Kernel function
kernel = 1.0 ./ (im * reshape(iw, (niw,1)) .- reshape(rmesh, (1,nmesh)))

# Build green's function
KA = kernel .* reshape(image, (1,nmesh))
giw = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = trapz(rmesh, KA[i,:])
end

# Build error
err = ones(F64, niw) * 1.0e-4

# Write green's function
open("giw.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%30.24f %30.24f %30.24f %30.24f\n", iw[i], real(z), imag(z), err[i])
    end
end

# Write spectral function
open("image.data", "w") do fout
    for i in eachindex(image)
        @printf(fout, "%20.16f %20.16f\n", rmesh[i], image[i])
    end
end

#!/usr/bin/env julia

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -10.0  # Left boundary
wmax = +10.0  # Right boundary
nmesh = 1000 # Number of real-frequency points
niw  = 50    # Number of Matsubara frequencies
beta = 100.0  # Inverse temperature
ϵ₁   = 2.50  # Parameters for gaussian peaks
ϵ₂   = -2.5
A₁   = 0.50
A₂   = 0.50
Γ₁   = 0.50
Γ₂   = 0.50

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Spectral function
gaussian(x, mu, sigma) = exp(-0.5*((x-mu)/sigma)^2)/(sqrt(2*π)*sigma)
rho(omega) = 0.8*gaussian(omega, -1.0, 1.0) + 0.2*gaussian(omega, 3, 0.7)
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
err = ones(F64, niw) * 0.0e-4

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

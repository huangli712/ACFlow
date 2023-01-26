#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +4.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 50.0  # Inverse temperature
nkx  = 50
nky  = 50
t    = 0.25
mune = -0.5
eta  = 0.05

function fermi(e)
    return 1.0 / ( exp((e - mune) * beta) + 1.0 )
end

function calc_ek(ikx, iky)
    kx = ikx * π / nkx
    ky = iky * π / nky
    ek = -2.0 * t * ( cos(kx) + cos(ky) )
    return ek
end

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))
ek = zeros(F64, 2 * nkx + 1, 2 * nky + 1)
chiw = zeros(C64, nmesh)
for ipx = -nkx:nkx
    for ipy = -nky:nky
        ek[ipx+nkx+1,ipy+nky+1] = calc_ek(ipx,ipy)
    end
end

for m = 1:nmesh
    iqx = 40
    iqy = 0
    r = 0.0
    k = 0
    w = rmesh[m] + im * eta

    for ipx = -nkx:nkx
        for ipy = -nky:nky
            ipqx = ipx + iqx
            if ipqx > +nkx
                ipqx = ipqx - 2 * nkx
            end
            if ipqx < -nkx
                ipqx = ipqx + 2 * nkx
            end

            ipqy = ipy + iqy
            if ipqy > +nky
                ipqy = ipqy - 2 * nky
            end
            if ipqy < -nky
                ipqy = ipqy + 2 * nky
            end

            ep = ek[ipx+nkx+1, ipy+nky+1]
            epq = ek[ipqx+nkx+1, ipqy+nky+1]
            
            k = k + 1
            r = r + ( fermi(ep) - fermi(epq) ) / (ep - epq + w)
        end
    end
    chiw[m] = r / k
end

open("chiw.data", "w") do fout
    for i in eachindex(chiw)
        z = chiw[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", rmesh[i], real(z), imag(z))
    end
end

error()

#
# For true spectrum
#


# Spectral function
image = similar(rmesh)
#
for i in eachindex(rmesh)
    A = W₁ / (1.0 + ((rmesh[i] - 0) / Γ₁) ^ 2.0) +
        W₂ / (1.0 + ((rmesh[i] - ϵ) / Γ₂) ^ 2.0) +
        W₂ / (1.0 + ((rmesh[i] + ϵ) / Γ₂) ^ 2.0)
    image[i] = A / (1.0 + (rmesh[i] / Γ₃) ^ 6.0)
end
#
image = image ./ trapz(rmesh, image)

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

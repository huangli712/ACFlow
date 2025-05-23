#!/usr/bin/env julia

haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

#
# This script is used to generate Lindhard function for square lattice.
#
# See Phys. Rev. B 94, 245140 (2016)
#

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +2.5  # Right boundary
nmesh = 501  # Number of real-frequency points
niw  = 20    # Number of Matsubara frequencies
beta = 50.0  # Inverse temperature
nkx  = 50    # Number of k-points along the x axis
nky  = 50    # Number of k-points along the y axis
t    = 0.25  # Hopping parameter
mune = -0.5  # Chemical potential
eta  = 0.05  # η

# Fermi distribution
function fermi(e)
    return 1.0 / ( exp((e - mune) * beta) + 1.0 )
end

# Define high-symmetry path: Γ - X - M - Γ
function calc_kpath()
    KGX = [(i,0) for i = 0:nkx]    # Γ - X
    KXM = [(nkx,i) for i = 0:nky]  # X - M
    KMG = [(i,i) for i = nkx:-1:0] # M - Γ
    KPATH = union(KGX, KXM, KMG)
    push!(KPATH, (0,0)) # Add the final Γ points
    return KPATH
end

# Band dispersion for square lattice
function calc_ek(ikx, iky)
    kx = ikx * π / nkx
    ky = iky * π / nky
    ek = -2.0 * t * ( cos(kx) + cos(ky) )
    return ek
end

# K-summation
function calc_ksum(iqx, iqy, w, ek)
    k = 0.0 # Counter for k-points
    r = zeros(C64, length(w))

    # Go through the brillouin zone
    for ipx = -nkx:nkx
        for ipy = -nky:nky
            # Increase the counter
            k = k + 1

            # Periodical boundary condition
            # Along x axis
            ipqx = ipx + iqx
            if ipqx > +nkx
                ipqx = ipqx - 2 * nkx
            end
            if ipqx < -nkx
                ipqx = ipqx + 2 * nkx
            end
            #
            # Along y axis
            ipqy = ipy + iqy
            if ipqy > +nky
                ipqy = ipqy - 2 * nky
            end
            if ipqy < -nky
                ipqy = ipqy + 2 * nky
            end

            # Evaluate the Lindhard function
            ep = ek[ipx+nkx+1,ipy+nky+1]
            epq = ek[ipqx+nkx+1,ipqy+nky+1]
            f₁ = fermi(ep) - fermi(epq)
            f₂ = ep - epq
            #
            @. r = r + f₁ / (f₂ + w)
        end
    end
    return r / k
end

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 0.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise_abs = randn(rng, F64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Evaluate KPATH
KPATH = calc_kpath()

# Precompute band dispersion
ek = zeros(F64, 2 * nkx + 1, 2 * nky + 1)
for ipx = -nkx:nkx
    for ipy = -nky:nky
        ek[ipx+nkx+1,ipy+nky+1] = calc_ek(ipx,ipy)
    end
end

# Announce the Lindhard function
chir = zeros(C64, length(KPATH), nmesh)
chiw = zeros(C64, length(KPATH), niw)

# Try to calculate the Lindhard function
for ik in eachindex(KPATH)
    iqx, iqy = KPATH[ik]
    println("ik = ", ik)

    # For real frequency data
    w = @. rmesh + im * eta
    chir[ik,:] = calc_ksum(iqx, iqy, w, ek)

    # For imaginary frequency data
    w = iωₙ * im .+ 1e-10 # Add a small real part to avoid NaN at iωₙ = 0
    chiw[ik,:] = calc_ksum(iqx, iqy, w, ek) + noise
end

# Write the Lindhard function at real frequency
open("chir.data", "w") do fout
    for k in eachindex(KPATH)
        for i in eachindex(rmesh)
            z = chir[k,i]
            @printf(fout, "%5i %20.16f %20.16f %20.16f\n", k, rmesh[i], real(z), imag(z))
        end
    end
end

# Write the Lindhard function at imaginary frequency
open("chiw.data", "w") do fout
    for k in eachindex(KPATH)
        for i in eachindex(iωₙ)
            z = chiw[k,i]
            @printf(fout, "%5i %20.16f %20.16f %20.16f\n", k, iωₙ[i], real(z), imag(z))
        end
    end
end

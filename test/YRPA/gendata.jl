#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = +0.0  # Left boundary
wmax = +4.0  # Right boundary
nmesh = 1001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 50.0  # Inverse temperature
nkx  = 10    # Number of k-points along the x axis
nky  = 10    # Number of k-points along the y axis
t    = 0.25  # Hopping parameter
mune = -0.5  # Chemical potential
eta  = 0.05  # η

# Fermi distribution
function fermi(e)
    return 1.0 / ( exp((e - mune) * beta) + 1.0 )
end

# Band dispersion for square lattice
function calc_ek(ikx, iky)
    kx = ikx * π / nkx
    ky = iky * π / nky
    ek = -2.0 * t * ( cos(kx) + cos(ky) )
    return ek
end

# Define high-symmetry path: Γ - X - M - Γ
KGX = [(i,0) for i = 0:nkx]    # Γ - X
KXM = [(nkx,i) for i = 0:nky]  # X - M
KMG = [(i,i) for i = nkx:-1:0] # M - Γ
#
KPATH = union(KGX, KXM, KMG)
#
push!(KPATH, (0,0)) # Add the final Γ points

# Real frequency mesh
rmesh = collect(LinRange(wmin, wmax, nmesh))

# Lindhard function
chiw = zeros(C64, length(KPATH), nmesh)

# Precompute band dispersion
ek = zeros(F64, 2 * nkx + 1, 2 * nky + 1)
for ipx = -nkx:nkx
    for ipy = -nky:nky
        ek[ipx+nkx+1,ipy+nky+1] = calc_ek(ipx,ipy)
    end
end

for m = 1:nmesh
    w = rmesh[m] + im * eta
    println("w -> ", w)

    for ik in eachindex(KPATH)
        iqx, iqy = KPATH[ik]
        r = 0.0
        k = 0

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

        chiw[ik,m] = r / k
    end
end

open("chiw.data", "w") do fout
    for i in eachindex(rmesh)
        for k in eachindex(KPATH)
            z = chiw[k,i]
            @printf(fout, "%5i %20.16f %20.16f %20.16f\n", k, rmesh[i], real(z), imag(z))
        end
    end
end

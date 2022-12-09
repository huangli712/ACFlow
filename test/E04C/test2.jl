#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using Random
using Printf
using ACFlow

# Setup parameters
wmin = -5.0  # Left boundary
wmax = +5.0  # Right boundary
nmesh = 2001 # Number of real-frequency points
niw  = 10    # Number of Matsubara frequencies
beta = 10.0  # Inverse temperature
ϵ₁   = 2.00  # Parameters for gaussian peaks
ϵ₂   = -2.0
ϵ₃   = 0.00
A₁   = 0.50
A₂   = 0.50
A₃   = 0.00
η    = 1e-2

# Real frequency mesh
ω = collect(LinRange(wmin, wmax, nmesh))

# Matsubara frequency mesh
iωₙ = π / beta * (2.0 * collect(0:niw-1) .+ 1.0)

# Noise
seed = rand(1:100000000)
rng = MersenneTwister(seed)
noise_ampl = 1.0e-4
noise_abs = randn(rng, F64, niw) * noise_ampl
noise_phase = rand(rng, niw) * 2.0 * π
noise = noise_abs .* exp.(noise_phase * im)

# Build green's function
giw = zeros(C64, niw)
for i in eachindex(giw)
    giw[i] = A₁ / (iωₙ[i] * im - ϵ₁) + A₂ / (iωₙ[i] * im - ϵ₂) + A₃ / (iωₙ[i] * im - ϵ₃) + noise[i]
end
#
gre = zeros(C64, nmesh)
for i in eachindex(gre)
    gre[i] = A₁ / (ω[i] + η * im - ϵ₁) + A₂ / (ω[i] + η * im - ϵ₂) + A₃ / (ω[i] + η * im - ϵ₃)
end

# Build error
err = ones(F64, niw) * noise_ampl

# Write green's function (Matsubara frequency axis)
open("giw.data", "w") do fout
    for i in eachindex(giw)
        z = giw[i]
        @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iωₙ[i], real(z), imag(z), err[i])
    end
end

# Write green's function (real frequency axis)
open("gre.data", "w") do fout
    for i in eachindex(gre)
        z = gre[i]
        @printf(fout, "%20.16f %20.16f %20.16f\n", ω[i], real(z), imag(z))
    end
end

function test_curve_fit(iωₙ::Vector{F64}, giw::Vector{C64}, η, ω)
    function fitfun(x, p)
        G = zeros(C64, niw)
        for i = 1:niw
            for j = 1:5
                G[i] = G[i] + p[2*j-1] / (x[i] * im - p[2*j])
            end
        end
        return real.(G) - real.(giw) + imag.(G) - imag.(giw) 
    end

    function fitfun1(x, p)
        G1 = zeros(C64, niw)
        for i = 1:niw
            for j = 1:5
                G1[i] = G1[i] + p[2*j-1] / (x[i] * im - p[2*j])
            end
        end
        return G1
    end

    function fitfun2(x, p)
        G2 = zeros(C64, nmesh)
        for i = 1:nmesh
            for j = 1:5
                G2[i] = G2[i] + p[2*j-1] / (x[i]  + η * im - p[2*j])
            end
        end
        return G2
    end

    fit = curve_fit(fitfun, iωₙ, zeros(F64,niw), [0.4, 2.2, 0.5, -2.3, 0.1, 0.2, 0.0, 0.0, 0.0, 0.0])
    G = fitfun1(iωₙ, fit.param)
    Gout = fitfun2(ω, fit.param)

    @show fit.param, sum(fit.resid)

    open("repr.data", "w") do fout
        for i in eachindex(G)
            z = G[i]
            @printf(fout, "%20.16f %20.16f %20.16f %20.16f\n", iωₙ[i], real(z), imag(z), err[i])
        end
    end

    open("Gout.data", "w") do fout
        for i in eachindex(Gout)
            z = Gout[i]
            @printf(fout, "%20.16f %20.16f %20.16f\n", ω[i], real(z), imag(z))
        end
    end
end

@show typeof(iωₙ), typeof(giw)
test_curve_fit(iωₙ, giw, η, ω)
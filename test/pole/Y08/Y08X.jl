#!/usr/bin/env julia

using Distributed
addprocs(8)

@everywhere using ACFlow
@everywhere using Printf

welcome()

# Setup key parameters
# This parameters should be consistent with those in gendata.jl
nkpt = 151
niw = 20
nmesh = 501

# Allocate memories at advance
chiw = zeros(C64, nkpt, niw)
Akw = zeros(F64, nkpt, nmesh)
grid = zeros(F64, niw)
mesh = zeros(F64, nmesh)

# Read momentum-resolved Lindhard function
open("chiw.data", "r") do fin
    for k = 1:nkpt
        for i = 1:niw
            arr = line_to_array(fin)
            _re, _im = parse.(F64, arr[3:4])
            grid[i] = parse(F64, arr[2])
            chiw[k,i] = _re + _im * im
        end
    end
end

# Analytically continuation k by k
for k = 1:nkpt

    # For StochPX solver

    # Setup parameters
    #
    # For [BASE] block
    # See types.jl/_PBASE for default setup
    B = Dict{String,Any}(
        "finput" => "chiw.data",
        "solver" => "StochPX",
        "ktype"  => "bsymm",
        "mtype"  => "flat",
        "grid"   => "bfreq",
        "mesh"   => "linear",
        "ngrid"  => 20,
        "nmesh"  => 501,
        "wmax"   => 2.5,
        "wmin"   => 0.0,
        "beta"   => 50.0,
    )
    #
    # For [StochPX] block
    # See types.jl/_PStochPX for default setup
    S = Dict{String,Any}(
        "method" => "mean",
        "nfine"  => 100000,
        "npole"  => 2000,
        "ntry"   => 2000,
        "nstep"  => 100,
        "theta"  => 1e+8,
        "eta"    => 1e-3,
    )
    #
    setup_param(B, S)

    # Call the solver
    mesh[:], Aout, Gout = solve(grid, chiw[k,:])

    # Store spectral density
    Akw[k,:] = Aout

    # Backup calculated results
    cp("Aout.data", "Aout.data.$k", force = true)
    cp("Gout.data", "Gout.data.$k", force = true)
    cp("repr.data", "repr.data.$k", force = true)

    println("k -> $k [finished]\n")
end

# Write Akw
open("Akw.data", "w") do fout
    for k = 1:nkpt
        for i = 1:nmesh
            @printf(fout, "%5i %20.16f %20.16f\n", k, mesh[i], Akw[k,i])
        end
    end
end

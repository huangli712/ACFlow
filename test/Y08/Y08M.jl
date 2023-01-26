#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACFlow

welcome()

# Setup key parameters
nkpt = 151
niw = 10

# Allocate memories at advance
chiw = zeros(C64, nkpt, niw)
grid = zeros(F64, niw)

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

    println("k -> ", k, " [finished]")

    # For MaxEnt solver

    # Setup parameters
    #
    # For [BASE] block
    # See types.jl/_PBASE for default setup
    B = Dict{String,Any}(
        "finput" => "chiw.data",
        "solver" => "MaxEnt",
        "ktype"  => "bsymm",
        "mtype"  => "flat",
        "grid"   => "bfreq",
        "mesh"   => "linear",
        "ngrid"  => 10,
        "nmesh"  => 501,
        "wmax"   => 2.5,
        "wmin"   => 0.0,
        "beta"   => 50.0,
    )
    #
    # For [MaxEnt] block
    # See types.jl/_PMaxEnt for default setup
    S = Dict{String,Any}(
        "method" => "chi2kink",
        "nalph"  => 14,
        "alpha"  => 1e12,
    )
    #
    setup_param(B, S)

    # Call the solver
    mesh, Aout, Gout = solve(grid, chiw[k,:])

    # Backup calculated results
    cp("Aout.data", "Aout.mem.$k", force = true)
    cp("Gout.data", "Gout.mem.$k", force = true)
    cp("repr.data", "repr.mem.$k", force = true)

end

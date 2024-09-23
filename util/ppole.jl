#!/usr/bin/env julia

#
# This script is used to regenerate the retarded Green's function and the
# corresponding spectral function by using the pole expansion.
#
# The η parameter is quite essential for the StochPX solver. It controls
# the broadening of the spectral function. Perhaps you would like to try
# a different η, but you might not want to perform a regular StochPX job
# again. At this time, this script can help. What you have to do is just
# to modify the η parameter (`eta`) in the `case.toml` file, and then
# make sure the existences of the `pole.data` file and `passed.data` file.
#
# Note that in addition to the `eta` parameter, the `method` parameter in
# the `case.toml` file can be changed as well. However, the `wmax` and
# `wmin` parameters should not be changed.
#
# Dr. Jia-Ming Wang (Renmin University) also contributes to this script.
#
# Usage:
#
#     $ ppole.jl ac.toml
#

push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src/")

using Printf
using ACFlow

"""
    calc_green_function(
        spe::StochPXElement,
        mesh::AbstractMesh,
        fmesh::AbstractMesh,
        Gᵥ::Vector{F64}
    )

Reconstruct Green's function at real axis by using the pole expansion. It
just calculates the contribution of the current solution (`spe`) to the
final Green's function.

Please see `calc_green()` in `src/spx.jl` for more details.

### Arguments
* spe -> A StochPXElement struct.
* mesh -> Mesh for output spectrum.
* fmesh -> Very dense mesh for the poles.
* Gᵥ -> Input data for correlator.

### Returns
* G -> Reconstructed Green's function, G(ω).
"""
function calc_green_function(
    spe::StochPXElement,
    mesh::AbstractMesh,
    fmesh::AbstractMesh,
    Gᵥ::Vector{F64}
    )
    ktype = get_b("ktype")

    χ₀ = -Gᵥ[1]
    @cswitch ktype begin
        @case "fermi"
            G = calc_green(spe.P, spe.A, spe.𝕊, mesh, fmesh)
            break

        @case "boson"
            G = calc_green(spe.P, spe.A, spe.𝕊, mesh, fmesh, χ₀, false)
            break

        @case "bsymm"
            G = calc_green(spe.P, spe.A, spe.𝕊, mesh, fmesh, χ₀, true)
            break
    end

    return G
end

"""
    parse_pole_data()

Try to parse the `pole.data` file, and return all the poles and χ².

### Arguments
N/A

### Returns
* χ²ᵥ -> All the χ².
* SPE -> A vector of StochPXElement struct. It contains all the poles.
"""
function parse_pole_data()
    # Extract key parameters
    ntry = get_x("ntry")
    npole = get_x("npole")

    # Prepare arrays
    SPE = StochPXElement[]
    P = zeros(I64, npole)
    A = zeros(F64, npole)
    𝕊 = zeros(F64, npole)
    χ²ᵥ = zeros(F64, ntry)

    # Check whether the `pole.data` file is available
    fn = "pole.data"
    @assert isfile(fn)

    open(fn, "r") do fin
        # There are `ntry` blocks in `pole.data` file.
        for i = 1:ntry
            # Extract χ²
            ldata = line_to_array(fin)
            χ²ᵥ[i] = parse(F64, ldata[5])
            #
            # Extract information about poles
            # For each block, there are `npole` lines
            for j = 1:npole
                ldata = line_to_array(fin)
                ind = parse(I64, ldata[1])
                @assert ind == j
                P[j] = parse(I64, ldata[2])
                A[j] = parse(F64, ldata[4])
                𝕊[j] = parse(F64, ldata[5])
            end
            #
            # Store the poles
            push!(SPE, StochPXElement(copy(P), copy(A), copy(𝕊)))
            #
            # Skip two blank lines between blocks
            readline(fin)
            readline(fin)
        end
    end

    return χ²ᵥ, SPE
end

"""
    filter_pole_data()

Parse the `passed.data` file, and tell us which solutions can be used to
reconstruct the Green's function and spectral function.

### Arguments
N/A

### Returns
* passed -> It contains the indices for selected solutions.
"""
function filter_pole_data()
    # Check whether the `passed.data` file is available.
    fn = "passed.data"
    @assert isfile(fn)

    # Prepare the array
    passed = I64[]

    # Read the data
    open(fn, "r") do fin
        ldata = line_to_array(fin)
        npass = parse(I64, ldata[3])
        for i = 1:npass
            ldata = line_to_array(fin)
            passed[i] = parse(I64, ldata[2])
        end
    end

    return passed
end

function pole_to_green()
    solver = get_b("solver")
    nmesh = get_b("nmesh")
    @assert solver == "StochPX"
    method = get_x("method")

    Gout = zeros(C64, nmesh)

    S = StochPXSolver()
    mesh = make_mesh()
    fmesh = calc_fmesh(S)
    Gᵥ, _ = init_iodata(S, read_data())

    χ²ᵥ, SPE = parse_pole_data()

    if method == "best"
        # The χ² of the best solution should be the smallest.
        p = argmin(χ²ᵥ)
        @printf("Best solution: try = %6i -> [χ² = %9.4e]\n", p, χ²ᵥ[p])
        #
        # Calculate G(ω)
        Gout = calc_green_function(SPE[p], mesh, fmesh, Gᵥ)
    else
        passed = filter_pole_data()
        for p in passed
            # Calculate and accumulate G(ω)
            G = calc_green_function(SPE[p], mesh, fmesh, Gᵥ)
            @. Gout = Gout + G
        end
        npass = length(passed)
        @. Gout = Gout / npass
        println("Accumulate $npass solutions to get the spectral density")
    end

    write_complete(mesh, Gout)
end

welcome()
overview()
read_param()
pole_to_green()
goodbye()

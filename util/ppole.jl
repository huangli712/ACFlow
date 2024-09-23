#!/usr/bin/env julia

#
# Usage:
#
#     $ ppole.jl ac.toml
#

push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src/")

using Printf
using ACFlow

"""
    calc_green(
        t::I64,
        SPE::Vector{StochPXElement},
        mesh::AbstractMesh,
        fmesh::AbstractMesh,
        Gᵥ::Vector{F64}
    )

Reconstruct Green's function at real axis by using the pole expansion. It
is a driver function.

### Arguments
* t -> Index of the current attempt.
* SPE -> A vector of StochPXElement. It contains all the poles.
* mesh -> Mesh for output spectrum.
* fmesh -> Very dense mesh for the poles.
* Gᵥ -> Input data for correlator.

### Returns
* G -> Reconstructed Green's function, G(ω).
"""
function calc_green(
    t::I64,
    SPE::Vector{StochPXElement},
    mesh::AbstractMesh,
    fmesh::AbstractMesh,
    Gᵥ::Vector{F64}
    )
    ktype = get_b("ktype")
    ntry = get_x("ntry")
    @assert t ≤ ntry

    # Calculate G(ω)
    χ₀ = -Gᵥ[1]
    @cswitch ktype begin
        @case "fermi"
            G = ACFlow.calc_green(SPE[t].P, SPE[t].A, SPE[t].𝕊, mesh, fmesh)
            break

        @case "boson"
            G = ACFlow.calc_green(SPE[t].P, SPE[t].A, SPE[t].𝕊, mesh, fmesh, χ₀, false)
            break

        @case "bsymm"
            G = ACFlow.calc_green(SPE[t].P, SPE[t].A, SPE[t].𝕊, mesh, fmesh, χ₀, true)
            break
    end

    return G
end

function parse_pole_data()
    ntry = get_x("ntry")
    npole = get_x("npole")

    SPE = StochPXElement[]
    P = zeros(I64, npole)
    A = zeros(F64, npole)
    𝕊 = zeros(F64, npole)
    χ²ᵥ = zeros(F64, ntry)

    fn = "pole.data"
    @assert isfile(fn)

    open(fn, "r") do fin
        for i = 1:ntry
            ldata = line_to_array(fin)
            χ²ᵥ[i] = parse(F64, ldata[5])
            for j = 1:npole
                ldata = line_to_array(fin)
                ind = parse(I64, ldata[1])
                @assert ind == j
                P[j] = parse(I64, ldata[2])
                A[j] = parse(F64, ldata[4])
                𝕊[j] = parse(F64, ldata[5])
            end
            push!(SPE, StochPXElement(copy(P), copy(A), copy(𝕊)))
            readline(fin)
            readline(fin)
        end
    end

    return χ²ᵥ, SPE
end

function filter_pole_data()
    fn = "passed.data"
    @assert isfile(fn)

    passed = I64[]

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
    @assert solver == "StochPX"
    method = get_x("method")

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
        Gout = calc_green(p, SPE, mesh, fmesh, Gᵥ)
    else
    end
end

welcome()
overview()
read_param()
pole_to_green()
goodbye()

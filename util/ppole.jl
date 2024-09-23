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
    calc_green(t::I64, SC::StochPXContext, real_axis::Bool)

Reconstruct Green's function at imaginary axis or real axis by using the
pole expansion. It is a driver function. If `real_axis = true`, it will
returns G(ω), or else G(iωₙ).

### Arguments
* t -> Index of the current attempt.
* SC -> A StochPXContext struct.
* real_axis -> Working at real axis (true) or imaginary axis (false)?

### Returns
* G -> Reconstructed Green's function, G(ω) or G(iωₙ).
"""
function calc_green(t::I64, SC::StochPXContext, real_axis::Bool)
    ktype = get_b("ktype")
    ntry = get_x("ntry")
    @assert t ≤ ntry

    # Calculate G(iωₙ)
    if real_axis == false
        return calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.Λ)
    end

    # Calculate G(ω). Now we don't need SC.Λ.
    χ₀ = -SC.Gᵥ[1]
    @cswitch ktype begin
        @case "fermi"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh)
            break

        @case "boson"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh, χ₀, false)
            break

        @case "bsymm"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh, χ₀, true)
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

    χ²ᵥ, SPE = parse_pole_data()

    if method == "best"
        # The χ² of the best solution should be the smallest.
        p = argmin(χ²ᵥ)
        @printf("Best solution: try = %6i -> [χ² = %9.4e]\n", p, χ²ᵥ[p])
        #
        # Calculate G(ω)
        Gout = calc_green(p, SC, true)
    else
    end
end

welcome()
overview()
read_param()
pole_to_green()
goodbye()

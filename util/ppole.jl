#!/usr/bin/env julia

#
# Usage:
#
#     $ ppole.jl ac.toml
#

push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src/")

using Printf
using ACFlow

function calc_green_function()
    η = get_x("eta")
end

function parse_green_data()
end

function parse_pole_data()
    ntry = get_x("ntry")
    npole = get_x("npole")

    SPE = StochPXElement[]
    P = zeros(I64, npole)
    A = zeros(F64, npole)
    𝕊 = zeros(F64, npole)
    χ² = zeros(F64, ntry)

    fn = "pole.data"
    @assert isfile(fn)

    open(fn, "r") do fin
        for i = 1:ntry
            ldata = line_to_array(fin)
            χ²[i] = parse(F64, ldata[5])
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

    return χ², SPE
end

function filter_pole_data()
end

function pole_to_green()
    solver = get_b("solver")
    ktype = get_b("ktype")
    @assert solver == "StochPX"

    parse_pole_data()
end

welcome()
overview()
read_param()
pole_to_green()
goodbye()

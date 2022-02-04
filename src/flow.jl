#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

abstract type AbstractData end

struct RawData{T} <: AbstractData
    _grid :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

mutable struct GreenData <: AbstractData
    value :: Vector{F64}
    error :: Vector{F64}
    covar :: Vector{F64}
end

function solve(rd::RawData)
    solver = get_c("solver")
    
    @cswitch solver begin
        @case "MaxEnt"
            MaxEnt.solve(rd)
            break

        @case "StochOM"
            break

        @case "StochAC"
            StochAC.solve(rd)
            break

        @default
            sorry()
            break
    end
end

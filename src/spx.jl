#
# Project : Gardenia
# Source  : spx.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/03
#

#=
### *Customized Structs* : *StochPX Solver*
=#

mutable struct StochPXElement
    P :: Vector{I64}
    A :: Vector{F64}
end

mutable struct StochPXContext
    Gᵥ    :: Vector{F64}
    Gᵧ    :: Vector{F64}
    σ¹    :: Vector{F64}
    allow :: Vector{I64}
    grid  :: AbstractGrid
    mesh  :: AbstractMesh
    Aout  :: Vector{F64}
end

function solve(S::StochPXSolver, rd::RawData)
end

function init(S::StochPXSolver, rd::RawData)
end

function run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
end

function prun(S::StochPXSolver,
            p1::Dict{String,Vector{Any}},
            p2::Dict{String,Vector{Any}},
            MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
end

function average(step::F64, SC::StochPXContext)
end

function last(SC::StochPXContext, Aout::Array{F64,2})
end

function warmup()
end

function sample()
end

function measure()
end

#=
### *Service Functions*
=#

function init_mc(S::StochPXSolver)
end

function init_element(S::StochPXSolver)
end

function init_iodata(S::StochPXSolver, rd::RawData)
end

function constraints(S::StochPXSolver)
end


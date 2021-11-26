#
# Project : Pansy
# Source  : ZenCore.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/19
#

abstract type AbstractData end

struct GreenData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

function GreenData()
    return GreenData(Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct SigmaData <: AbstractData
    
end

function SigmaData()
end

struct ChiData <: AbstractData end

function ChiData()
end

struct MomentsData <: AbstractData
    𝑀₀ :: F64
    𝑀₁ :: F64
    𝑀₂ :: F64
    𝑀₃ :: F64
end

function MomentsData()
    return MomentsData(zero(F64), zero(F64), zero(F64), zero(F64))
end

abstract type AbstractGrid end

struct RealFrequencyGrid <: AbstractGrid end

function RealFrequencyGrid()
end

struct ImaginaryTimeGrid <: AbstractGrid
    grid :: Vector{F64}
end

function ImaginaryTimeGrid()
end

struct FermionicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

function FermionicMatsubaraGrid()
    return FermionicMatsubaraGrid(Vector{F64}[])
end

struct BosonicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

function BosonicMatsubaraGrid()
end

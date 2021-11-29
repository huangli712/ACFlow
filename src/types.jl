#
# Project : Gardenia
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/29
#

abstract type AbstractData end
abstract type AbstractGrid end

struct GreenData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

function GreenData()
    return GreenData(Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct SigmaData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}    
end

function SigmaData()
    return SigmaData(Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct ChiData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

function ChiData()
    return ChiData(Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct MomentsData <: AbstractData
    𝑀₀ :: N64
    𝑀₁ :: N64
    𝑀₂ :: N64
    𝑀₃ :: N64
end

function MomentsData(::T) where {T <: N64}
    return MomentsData(zero(T), zero(T), zero(T), zero(T))
end

struct VectorMomentsData <: AbstractData
    V𝑀₀ :: Vector{N64}
    V𝑀₁ :: Vector{N64}
    V𝑀₂ :: Vector{N64}
    V𝑀₃ :: Vector{N64}
end

function VectorMomentsData(::T) where {T <: N64}
    return VectorMomentsData(Vector{N64}[], Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct ImaginaryTimeGrid <: AbstractGrid
    grid :: Vector{F64}
end

function ImaginaryTimeGrid()
    return ImaginaryTimeGrid(Vector{F64}[])
end

struct RealFrequencyGrid <: AbstractGrid
    grid :: Vector{F64}
end

function RealFrequencyGrid()
    return RealFrequencyGrid(Vector{F64}[])
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
    return BosonicMatsubaraGrid(Vector{F64}[])
end

abstract type AbstractGrid end

struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

struct BosonicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

abstract type AbstractMesh end

struct UniformMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

struct NonUniformMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

abstract type AbstractData end

struct RawData{T} <: AbstractData
    mesh  :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

mutable struct GreenData <: AbstractData
    value :: Vector{F64}
    error :: Vector{F64}
    var   :: Vector{F64}
end
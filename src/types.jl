abstract type AbstractMesh end

struct UniformMesh <: AbstractMesh
    nmesh :: I64
    wmax  :: F64
    wmin  :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

struct NonUniformMesh <: AbstractMesh end

abstract type AbstractGrid end

struct FermionicImaginaryTimeGrid <: AbstractGrid
    ngrid :: I64
    β :: F64
    τ :: Vector{F64}
end

struct FermionicMatsubaraGrid <: AbstractGrid
    ngrid :: I64
    β :: F64
    ω :: Vector{F64}
end

struct BosonicImaginaryTimeGrid <: AbstractGrid end

struct BosonicMatsubaraGrid <: AbstractGrid end

abstract type AbstractData end

struct RawData{T} <: AbstractData
    mesh  :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

mutable struct GreenData <: AbstractData
    value :: Vector{C64}
    error :: Vector{F64}
    var   :: Vector{F64}
end
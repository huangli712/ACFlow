#
# Project : Gardenia
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/16
#

abstract type AbstractData end
abstract type AbstractGrid end

struct GreenData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

struct SigmaData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}    
end

struct ChiData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

struct MomentsData <: AbstractData
    𝑀₀ :: N64
    𝑀₁ :: N64
    𝑀₂ :: N64
    𝑀₃ :: N64
    𝐶𝑀 :: Matrix{N64}
end

#=
struct VectorMomentsData <: AbstractData
    V𝑀₀ :: Vector{N64}
    V𝑀₁ :: Vector{N64}
    V𝑀₂ :: Vector{N64}
    V𝑀₃ :: Vector{N64}
end
=#

struct KernelData <: AbstractData
    𝐾 :: Matrix{N64}
end

struct KernelMomentsData <: AbstractData
    𝐾𝑀₀ :: Matrix{N64}
    𝐾𝑀₁ :: Matrix{N64}
    𝐾𝑀₂ :: Matrix{N64}
    𝐾𝑀₃ :: Matrix{N64}
end

struct ImaginaryTimeGrid <: AbstractGrid
    grid :: Vector{F64}
end

struct RealFrequencyGrid <: AbstractGrid
    nul :: I64
    nur :: I64
    w0l :: F64
    wl  :: F64
    w0r :: F64
    wr  :: F64
    dw  :: F64
    grid :: Vector{F64}
end

struct FermionicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

struct BosonicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

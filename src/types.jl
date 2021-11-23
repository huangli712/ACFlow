
abstract type AbstractData end
struct GreenData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

function GreenData()
    return GreenData(Vector{N64}[], Vector{N64}[], Vector{N64}[])
end

struct SigmaData <: AbstractData end
struct ChiData <: AbstractData end

abstract type AbstractGrid end
struct RealFrequencyGrid <: AbstractGrid end

struct ImaginaryTimeGrid <: AbstractGrid
    grid :: Vector{F64}
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


abstract type AbstractData end
struct GreenData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
end

struct SigmaData <: AbstractData end
struct ChiData <: AbstractData end




abstract type AbstractGrid end
struct RealFrequencyGrid <: AbstractGrid end

struct ImaginaryTimeGrid <: AbstractGrid
    β :: F64
    grid :: Vector{F64}
end

struct FermionicMatsubaraGrid <: AbstractGrid
    β :: F64
    grid :: Vector{F64}
end

struct BosonicMatsubaraGrid <: AbstractGrid
    β :: F64
    grid :: Vector{F64}
end

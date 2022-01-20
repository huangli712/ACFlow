abstract type AbstractMesh end
struct UniformMesh <: AbstractMesh end
struct NonUniformMesh <: AbstractMesh end

abstract type AbstractGrid end
struct FermionicImaginaryTimeGrid <: AbstractGrid end

struct FermionicMatsubaraGrid <: AbstractGrid
    ngrid :: I64
    β :: F64
    ω :: Vector{F64}
end

function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ngrid
    return fg.ω[ind]
end

struct BosonicImaginaryTimeGrid <: AbstractGrid end

struct BosonicMatsubaraGrid <: AbstractGrid
    ngrid :: I64
    β :: F64
    ω :: Vector{F64}
end

function Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ngrid
    return bg.ω[ind]
end
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

function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ngrid
    return fg.ω[ind]
end

function make_grid(fg::FermionicMatsubaraGrid)
    for n in eachindex(fg.ω)
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end